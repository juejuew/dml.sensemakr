# Compatibility is established from the saved panel and a cell ATT/IF replay.
# Read the saved estimation sample; never rerun global preprocessing.
.did_adapter_context <- function(fit) {
  # An installed version is diagnostic metadata, not the producer of a saved fit.
  version <- .adapter_package_version("did")
  drdid.version <- .drdid_require()
  if (!is.list(fit) || !inherits(fit, "MP"))
    .drdid_stop("fit must be an MP object from did::att_gt(); aggte() is not supported.")
  dp <- fit$DIDparams
  if (!is.list(dp) || !is.data.frame(dp$data))
    .drdid_stop("The source object must retain DIDparams$data.")
  required <- c("yname", "tname", "idname", "gname", "xformla", "panel",
                "est_method", "faster_mode", "control_group", "base_period", "anticipation")
  missing <- required[vapply(required, function(key) is.null(dp[[key]]), logical(1))]
  if (length(missing))
    .drdid_stop("The source object lacks required DIDparams fields: ", paste(missing, collapse = ", "), ".")
  is.string <- function(x) is.character(x) && length(x) == 1L && !is.na(x) && nzchar(x)
  column.fields <- c("idname", "tname", "gname", "yname")
  if (!all(vapply(dp[column.fields], is.string, logical(1))))
    .drdid_stop("Saved ID, period, cohort and outcome column names must be nonempty strings.")
  if (!isTRUE(dp$panel) || isTRUE(dp$true_repeated_cross_sections))
    .drdid_stop("Only the balanced panel estimator path is supported, not repeated cross-sections or unbalanced estimation.")
  if (identical(dp$fix_weights, "varying"))
    .drdid_stop("fix_weights = 'varying' uses the repeated-cross-section estimator, even for balanced panels.")
  if (!is.null(dp$fix_weights) && !dp$fix_weights %in% c("base_period", "first_period"))
    .drdid_stop("Unknown fix_weights setting.")
  if (!identical(dp$est_method, "dr"))
    .drdid_stop("Only est_method = 'dr' is supported; custom, reg and ipw estimators are excluded.")
  if (length(setdiff(dp$clustervars, dp$idname)))
    .drdid_stop("Clustering above the individual level requires clustered joint sensitivity inference.")
  if (!is.logical(dp$faster_mode) || length(dp$faster_mode) != 1L || is.na(dp$faster_mode) ||
      !identical(dp$control_group, "nevertreated") && !identical(dp$control_group, "notyettreated") ||
      !identical(dp$base_period, "varying") && !identical(dp$base_period, "universal"))
    .drdid_stop("Missing or unsupported saved estimation settings.")
  a <- dp$anticipation
  if (!is.numeric(a) || length(a) != 1L || !is.finite(a) || a < 0 || a != floor(a))
    .drdid_stop("The saved anticipation must be a nonnegative integer.")
  dat <- as.data.frame(dp$data)
  columns <- c(dp$idname, dp$tname, dp$gname, dp$yname)
  if (length(columns) != 4L || anyNA(columns) || any(!columns %in% names(dat)) || !nrow(dat))
    .drdid_stop("The saved estimation data lack required columns or observations.")
  id <- dat[[dp$idname]]; tt <- dat[[dp$tname]]; gg <- dat[[dp$gname]]
  if (!is.numeric(id) || any(!is.finite(id)) || !is.numeric(tt) || any(!is.finite(tt)) ||
      !is.numeric(gg) || anyNA(gg) || any(gg == -Inf) ||
      !is.numeric(dat[[dp$yname]]) || any(!is.finite(dat[[dp$yname]])))
    .drdid_stop("The saved panel contains invalid IDs, periods, cohorts or outcomes.")
  ids <- sort(unique(id)); times <- sort(unique(tt)); N <- length(ids)
  if (length(times) < 2L || anyDuplicated(dat[c(dp$idname, dp$tname)]) ||
      nrow(dat) != N * length(times) || any(tabulate(match(id, ids), N) != length(times)))
    .drdid_stop("DIDparams$data must be a balanced panel with unique ID-time pairs.")
  # Canonical ID order, independent of the fast/slow source row order.
  rows <- lapply(times, function(t) { r <- which(tt == t); r[match(ids, id[r])] })
  cohorts <- gg[rows[[1L]]]
  if (any(vapply(rows, function(r) !identical(gg[r], cohorts), logical(1))))
    .drdid_stop("Cohort membership changes within a saved panel unit.")
  if (!is.null(dp$weightsname)) {
    if (!dp$weightsname %in% names(dat)) .drdid_stop("The saved sampling-weight column is missing.")
    .drdid_check_weights(dat[[dp$weightsname]], nrow(dat))
  }
  if (!".w" %in% names(dat)) .drdid_stop("The saved normalized sampling weights are missing.")
  .drdid_check_weights(dat$.w, nrow(dat))
  if (!inherits(dp$xformla, "formula") || length(dp$xformla) != 2L)
    .drdid_stop("The saved xformla must be a one-sided formula.")
  saved.times <- if (dp$faster_mode) dp$time_periods else dp$tlist
  if (!identical(as.numeric(saved.times), as.numeric(times)))
    .drdid_stop("Saved periods do not match DIDparams$data.")
  if (!is.numeric(fit$n) || length(fit$n) != 1L || !is.finite(fit$n) || fit$n != N)
    .drdid_stop("The source sample size does not match the saved panel units.")
  m <- length(fit$att)
  if (!is.numeric(fit$att) || !m || !is.numeric(fit$group) || !is.numeric(fit$t) ||
      length(fit$group) != m || length(fit$t) != m || any(!is.finite(c(fit$group, fit$t))) ||
      anyDuplicated(data.frame(group = fit$group, time = fit$t)))
    .drdid_stop("The source ATT cells must have unique, finite group-time labels.")
  inf <- fit$inffunc
  if (is.null(inf)) .drdid_stop("The source has no IF. Refit did::att_gt() with compute_inffunc = TRUE.")
  if (length(dim(inf)) != 2L || !identical(as.integer(dim(inf)), c(as.integer(N), as.integer(m))))
    .drdid_stop("The source IF dimensions do not match its units and ATT cells.")
  inf.ids <- rownames(inf)
  if (is.null(inf.ids) || anyNA(inf.ids) || anyDuplicated(inf.ids) ||
      !setequal(inf.ids, as.character(ids)))
    .drdid_stop("The source IF must have unique row names matching the saved unit IDs.")
  # The fast path saves the exact global design basis. The slow path builds
  # model.matrix once on ALL retained periods before selecting a cell.
  if (dp$faster_mode) {
    ti <- as.data.frame(dp$time_invariant_data)
    source.ids <- ti[[dp$idname]]
    if (length(source.ids) != N || anyNA(source.ids) || anyDuplicated(source.ids) ||
        !setequal(source.ids, ids)) .drdid_stop("Saved tensor unit IDs cannot be aligned.")
    ix <- match(ids, source.ids)
    if (length(dp$covariates_tensor) != length(times) ||
        length(dp$outcomes_tensor) != length(times) || length(dp$weights_tensor) != length(times))
      .drdid_stop("The saved panel tensors have incompatible dimensions.")
    X <- lapply(seq_along(times), function(k) {
      r <- rows[[k]]
      if (!isTRUE(all.equal(as.numeric(dp$outcomes_tensor[[k]][ix]), dat[[dp$yname]][r], tolerance = 1e-12)) ||
          !isTRUE(all.equal(as.numeric(dp$weights_tensor[[k]][ix]), dat$.w[r], tolerance = 1e-12)))
        .drdid_stop("The saved outcomes/weights tensors disagree with DIDparams$data.")
      mat <- as.matrix(dp$covariates_tensor[[k]])
      if (nrow(mat) != N) .drdid_stop("The saved covariate tensor has an incompatible unit count.")
      mat[ix, , drop = FALSE]
    })
  } else {
    for (v in all.vars(dp$xformla)) {
      if (!v %in% names(dat)) .drdid_stop("A formula variable is absent from the saved data.")
      if (is.factor(dat[[v]])) dat[[v]] <- droplevels(dat[[v]])
    }
    mm <- stats::model.matrix(dp$xformla, data = dat, na.action = stats::na.pass)
    if (nrow(mm) != nrow(dat)) .drdid_stop("Evaluating the saved formula changed the estimation sample.")
    X <- lapply(rows, function(r) mm[r, , drop = FALSE])
  }
  for (k in seq_along(X)) {
    mat <- X[[k]]
    if (!is.numeric(mat) || !ncol(mat) || any(!is.finite(mat)))
      .drdid_stop("The saved design must be finite and have at least one column.")
    if (is.null(colnames(mat))) colnames(mat) <- if (ncol(mat) == 1L && all(mat == 1))
      "(Intercept)" else paste0("x", seq_len(ncol(mat)))
    if (anyNA(colnames(mat)) || any(!nzchar(colnames(mat))) || anyDuplicated(colnames(mat)))
      .drdid_stop("The saved design needs unique, nonempty column names.")
    rownames(mat) <- NULL
    X[[k]] <- mat
  }
  list(fit = fit, dp = dp, data = dat, ids = ids, times = times, rows = rows,
       cohorts = cohorts, X = X, if.rows = match(as.character(ids), inf.ids),
       N = N, version = version, drdid.version = drdid.version,
       layout = if (dp$faster_mode) "panel_tensors" else "panel_formula")
}

.did_cell_sample <- function(ctx, group, time) {
  scalar <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x)
  if (!scalar(group) || !scalar(time)) .drdid_stop("group and time must be finite numeric scalars.")
  if (time < group) .drdid_stop("Only post-treatment cells (time >= group) are supported.")
  j <- which(ctx$fit$group == group & ctx$fit$t == time)
  if (length(j) != 1L) .drdid_stop("The requested (group, time) cell is absent from the MP object.")
  theta <- ctx$fit$att[j]
  if (!is.finite(theta)) .drdid_stop("The requested cell has no finite source ATT.")
  pre <- which(ctx$times + ctx$dp$anticipation < group)
  if (!length(pre)) .drdid_stop("The requested cohort has no valid pre-treatment base period.")
  pre <- utils::tail(pre, 1L); post <- match(time, ctx$times)
  if (is.na(post)) .drdid_stop("The requested time is absent from the saved data.")
  treated <- ctx$cohorts == group
  never <- ctx$cohorts == 0 | ctx$cohorts == Inf
  control <- if (ctx$dp$control_group == "nevertreated") never else
    never | (ctx$cohorts > time + ctx$dp$anticipation & !treated)
  take <- which(treated | control)
  D <- as.numeric(treated[take]); n <- length(take)
  if (n < 4L || length(unique(D)) != 2L) .drdid_stop("The cell needs both treated and control units.")
  phi.full <- as.numeric(ctx$fit$inffunc[, j])
  if (length(phi.full) != ctx$N || any(!is.finite(phi.full)))
    .drdid_stop("The selected source IF column must be finite.")
  if.rows <- ctx$if.rows[take]
  outside <- setdiff(seq_len(ctx$N), if.rows)
  if (any(abs(phi.full[outside]) > 1e-10))
    .drdid_stop("The source IF has nonzero contributions outside the reconstructed cell.")
  # did stores (N / n_cell) * phi_cell. Undo this before cell-level DML.
  phi <- (n / ctx$N) * phi.full[if.rows]
  w.period <- if (identical(ctx$dp$fix_weights, "first_period")) 1L else pre
  y <- ctx$data[[ctx$dp$yname]]
  sample <- list(y1 = y[ctx$rows[[post]][take]], y0 = y[ctx$rows[[pre]][take]],
                 D = D, X = ctx$X[[pre]][take, , drop = FALSE],
                 i.weights = ctx$data$.w[ctx$rows[[w.period]][take]],
                 unit_id = ctx$ids[take], id.supplied = TRUE,
                 input.info = list(format = "did saved panel", formula = ctx$dp$xformla,
                   times = ctx$times[c(pre, post)], source_rows = list(
                     pre = ctx$rows[[pre]][take], post = ctx$rows[[post]][take])))
  sample$i.weights <- sample$i.weights / mean(sample$i.weights)
  spec <- list(method = "trad", trim.level = 0.995, phi = phi)
  diagnostics <- .panel_check_source(theta, sample, spec)
  list(sample = sample, theta = theta, phi = phi, diagnostics = diagnostics,
       cell = list(group = group, time = time, base_period = ctx$times[pre],
                   column = j, n_cell = n, n_global = ctx$N, if_rows = if.rows,
                   if_scale = n / ctx$N))
}
