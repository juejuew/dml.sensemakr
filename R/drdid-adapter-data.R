# Input and source-fit checks for the panel DRDID adapter.
.drdid_stop <- function(...) stop(..., call. = FALSE)

.adapter_package_version <- function(package) {
  tryCatch(as.character(utils::packageVersion(package)), error = function(e) NA_character_)
}

.drdid_require <- function() {
  if (!requireNamespace("DRDID", quietly = TRUE))
    .drdid_stop("Install DRDID to reconstruct and check the source panel fit.")
  .adapter_package_version("DRDID")
}

# Check only the interfaces used by the selected path; never infer compatibility
# from a version number or assume that an absent field has a new default value.
.drdid_function <- function(name, required) {
  fun <- tryCatch(utils::getFromNamespace(name, "DRDID"), error = function(e) NULL)
  if (!is.function(fun)) .drdid_stop("DRDID lacks the required function: ", name, ".")
  missing <- setdiff(required, names(formals(fun)))
  if (length(missing))
    .drdid_stop("Incompatible DRDID interface ", name, ": missing arguments ",
                paste(missing, collapse = ", "), ".")
  fun
}

.drdid_spec <- function(fit) {
  if (!inherits(fit, "drdid") || !identical(fit$argu$type, "dr"))
    .drdid_stop("fit must be a DRDID doubly robust result (argu$type = 'dr').")
  if (!isTRUE(fit$argu$panel))
    .drdid_stop("Only two-period panel data are supported; repeated cross-sections require different sensitivity moments.")
  method <- fit$argu$estMethod[1L]
  if (length(method) != 1L || is.na(method) || !method %in% c("imp", "trad"))
    .drdid_stop("Cannot establish the source estimation method; use 'imp' or 'trad'.")
  if (!is.numeric(fit$ATT) || length(fit$ATT) != 1L || !is.finite(fit$ATT))
    .drdid_stop("fit$ATT must be a finite scalar.")
  phi <- fit[["att.inf.func"]]
  if (is.null(phi))
    .drdid_stop("The source fit has no influence function. Refit DRDID with inffunc = TRUE.")
  if (!is.numeric(phi) || any(!is.finite(phi)) ||
      (!is.null(dim(phi)) && (length(dim(phi)) != 2L || ncol(phi) != 1L)))
    .drdid_stop("fit$att.inf.func must be a finite vector or single-column matrix.")
  trim <- fit$argu$trim.level
  if (!is.numeric(trim) || length(trim) != 1L || !is.finite(trim) || trim <= 0 || trim > 1)
    .drdid_stop("The source fit must record a trim.level in (0, 1].")
  if (method == "imp" &&
      (length(fit$ps.flag) != 1L || is.na(fit$ps.flag) || !fit$ps.flag %in% c(0, 1)))
    .drdid_stop("The improved source fit must have successful IPT estimation (ps.flag 0 or 1); MLE fallback is not supported.")
  list(method = method, trim.level = trim, phi = as.numeric(phi))
}

.drdid_column <- function(value, call, name) {
  if (is.null(value) && !is.null(call)) {
    candidate <- call[[name]]
    if (is.character(candidate) && length(candidate) == 1L) value <- candidate
  }
  if (!is.character(value) || length(value) != 1L || is.na(value) || !nzchar(value))
    .drdid_stop("Supply ", name, " explicitly; a call expression is not a resolved column name.")
  value
}

.drdid_check_weights <- function(weights, n) {
  if (!is.numeric(weights) || length(weights) != n || any(!is.finite(weights)) ||
      any(weights <= 0) || any(weights != weights[1L]))
    .drdid_stop("Only constant positive sampling weights are supported; weighted sensitivity moments require a separate extension.")
  invisible(NULL)
}

.drdid_panel_data <- function(fit, spec, data, yname, tname, idname, dname,
                              xformla, weightsname, y1, y0, D, covariates,
                              i.weights, unit_id) {
  source.call <- fit[["call.param"]]
  if (!is.null(data)) {
    if (any(!vapply(list(y1, y0, D, covariates, i.weights, unit_id), is.null, logical(1))))
      .drdid_stop("Supply either long data with column names, or wide panel vectors, not both.")
    if (!is.data.frame(data)) .drdid_stop("data must be a data.frame.")
    yname <- .drdid_column(yname, source.call, "yname")
    tname <- .drdid_column(tname, source.call, "tname")
    idname <- .drdid_column(idname, source.call, "idname")
    dname <- .drdid_column(dname, source.call, "dname")
    if (is.null(weightsname) && !is.null(source.call[["weightsname"]]))
      weightsname <- .drdid_column(NULL, source.call, "weightsname")
    if (!is.null(weightsname))
      weightsname <- .drdid_column(weightsname, source.call, "weightsname")
    required <- c(yname, tname, idname, dname, weightsname)
    if (any(!required %in% names(data))) .drdid_stop("Some specified columns are absent from data.")
    # Check the original values across both periods: DRDID normalizes weights
    # before returning them, which can otherwise turn all-negative weights positive.
    if (!is.null(weightsname))
      .drdid_check_weights(data[[weightsname]], nrow(data))
    if (is.null(xformla)) {
      candidate <- source.call[["xformla"]]
      if (is.null(candidate)) xformla <- ~ 1
      else if (inherits(candidate, "formula") ||
               (is.call(candidate) && identical(candidate[[1L]], as.name("~"))))
        xformla <- stats::as.formula(candidate)
      else .drdid_stop("Supply xformla explicitly; the original formula was saved as an expression.")
    }
    if (!inherits(xformla, "formula") || length(xformla) != 2L)
      .drdid_stop("xformla must be a one-sided formula.")
    id <- data[[idname]]; time <- data[[tname]]
    if (!is.numeric(id) || !is.numeric(time) || any(!is.finite(id)) || any(!is.finite(time)))
      .drdid_stop("Long panel IDs and times must be finite numeric values, as required by DRDID.")
    times <- sort(unique(time))
    if (length(times) != 2L || anyDuplicated(data[c(idname, tname)]) || any(table(id) != 2L))
      .drdid_stop("Provide the exact balanced two-period estimation sample with unique ID-time pairs; the adapter does not drop units.")
    X <- stats::model.matrix(xformla, stats::model.frame(xformla, data,
                                                       na.action = stats::na.pass))
    if (any(!is.finite(X)) || anyNA(data[required]))
      .drdid_stop("Missing values in estimation variables are not supported. Prepare the common sample and refit DRDID first.")
    prep <- .drdid_function("pre_process_drdid", c("yname", "tname", "idname", "dname",
      "xformla", "data", "panel", "estMethod", "weightsname", "boot", "inffunc"))
    dp <- prep(yname = yname, tname = tname, idname = idname, dname = dname,
               xformla = xformla, data = data, panel = TRUE, estMethod = spec$method,
               weightsname = weightsname, boot = FALSE, inffunc = TRUE)
    y1 <- dp$y1; y0 <- dp$y0; D <- dp$D
    covariates <- as.matrix(dp$covariates); i.weights <- dp$i.weights
    unit_id <- sort(unique(id))
    if (length(y1) != length(unit_id) || length(y0) != length(unit_id))
      .drdid_stop("DRDID preprocessing changed the supplied sample; cannot establish the IF mapping.")
    input.info <- list(format = "long", columns = list(yname = yname, tname = tname,
                        idname = idname, dname = dname, weightsname = weightsname),
                       formula = xformla, times = times,
                       source_rows = list(
                         pre = which(time == times[1L])[match(unit_id, id[time == times[1L]])],
                         post = which(time == times[2L])[match(unit_id, id[time == times[2L]])]))
  } else {
    if (any(!vapply(list(yname, tname, idname, dname, xformla, weightsname), is.null, logical(1))))
      .drdid_stop("Column names and xformla require long data; use covariates for wide input.")
    if (is.null(y1) || is.null(y0) || is.null(D))
      .drdid_stop("Supply y1, y0 and D, or supply long data with column names.")
    input.info <- list(format = "wide")
  }
  n <- length(D)
  vectors <- list(y1 = y1, y0 = y0, D = D)
  if (n < 4L || any(vapply(vectors, function(v)
      !is.numeric(v) || length(v) != n || any(!is.finite(v)), logical(1))))
    .drdid_stop("y1, y0 and D must be finite numeric vectors of the same length.")
  if (!all(D %in% c(0, 1)) || length(unique(D)) != 2L)
    .drdid_stop("D must contain both treatment groups, coded 0 and 1.")
  if (is.null(covariates)) covariates <- matrix(1, n, 1L, dimnames = list(NULL, "(Intercept)"))
  X <- as.matrix(covariates)
  if (!is.numeric(X) || nrow(X) != n || !ncol(X) || any(!is.finite(X)))
    .drdid_stop("covariates must be a finite numeric matrix with one row per unit.")
  if (is.null(colnames(X))) colnames(X) <- paste0("x", seq_len(ncol(X)))
  if (anyNA(colnames(X)) || any(!nzchar(colnames(X))) || anyDuplicated(colnames(X)))
    .drdid_stop("Covariate columns need unique, nonempty names.")
  if (is.null(i.weights)) i.weights <- rep(1, n)
  .drdid_check_weights(i.weights, n)
  id.supplied <- !is.null(unit_id)
  if (!id.supplied) unit_id <- seq_len(n)
  if (length(unit_id) != n || anyNA(unit_id) || anyDuplicated(unit_id))
    .drdid_stop("unit_id must identify each supplied unit exactly once, in source IF order.")
  if (length(spec$phi) != n) .drdid_stop("The source IF length does not match the supplied panel units.")
  list(y1 = as.numeric(y1), y0 = as.numeric(y0), D = as.numeric(D), X = X,
       i.weights = i.weights / mean(i.weights), unit_id = unit_id,
       input.info = input.info, id.supplied = id.supplied)
}

.drdid_refit_source <- function(sample, spec) {
  name <- if (spec$method == "imp") "drdid_imp_panel" else "drdid_panel"
  fun <- .drdid_function(name, c("y1", "y0", "D", "covariates", "i.weights",
                                 "boot", "inffunc", "trim.level"))
  out <- fun(y1 = sample$y1, y0 = sample$y0, D = sample$D, covariates = sample$X,
      i.weights = sample$i.weights, boot = FALSE, inffunc = TRUE,
      trim.level = spec$trim.level)
  if (!is.list(out) || !is.numeric(out$ATT) || length(out$ATT) != 1L ||
      !is.finite(out$ATT) || !is.numeric(out$att.inf.func) ||
      length(out$att.inf.func) != length(sample$D) || any(!is.finite(out$att.inf.func)))
    .drdid_stop("DRDID reconstruction must return a finite ATT and one finite IF value per unit.")
  out
}

.drdid_check_source <- function(fit, sample, spec) {
  .panel_check_source(fit$ATT, sample, spec, fit$ps.flag)
}

.panel_check_source <- function(theta, sample, spec, ps.flag = NULL) {
  replay <- .drdid_refit_source(sample, spec)
  close <- function(x, y) isTRUE(all.equal(as.numeric(x), as.numeric(y),
                                          tolerance = 1e-7, check.attributes = FALSE))
  if (!close(theta, replay$ATT) || !close(spec$phi, replay$att.inf.func))
    .drdid_stop("The supplied data/settings do not reproduce the source ATT and ordered IF. Check sample, weights, covariates and row order.")
  if (spec$method == "imp") {
    ps.fun <- .drdid_function("pscore.cal", c("D", "int.cov", "i.weights", "n"))
    ps <- ps.fun(
      sample$D, sample$X, i.weights = sample$i.weights, n = length(sample$D))
    if (!is.list(ps) || length(ps$flag) != 1L || !ps$flag %in% c(0, 1) ||
        !identical(as.numeric(ps$flag), as.numeric(ps.flag)))
      .drdid_stop("The source IPT convergence branch could not be reproduced.")
    pscore <- ps$pscore
  } else {
    ps.fun <- .drdid_function("fastglm_fit", c("x", "y", "family", "weights", "method"))
    ps <- ps.fun(x = sample$X, y = sample$D, family = stats::binomial(),
                 weights = sample$i.weights, method = 3)
    if (!is.list(ps) || !isTRUE(ps$converged)) .drdid_stop("Source propensity estimation did not converge.")
    pscore <- ps$fitted.values
  }
  if (!is.numeric(pscore) || length(pscore) != length(sample$D) || any(!is.finite(pscore)))
    .drdid_stop("DRDID propensity reconstruction must return one finite probability per unit.")
  pscore <- as.numeric(pscore)
  if (any(!is.finite(pscore)) || any(pscore <= 0 | pscore >= 1 - 1e-6) ||
      any(pscore[sample$D == 0] >= spec$trim.level))
    .drdid_stop("Source propensity trimming/clipping is active. This adapter currently supports the untrimmed panel estimand only.")
  list(reproduced = TRUE, propensity = pscore, trim.active = FALSE,
       analytic.se = replay$se)
}
