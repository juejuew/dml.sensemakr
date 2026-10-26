#' Bounds benchmarked against a pre-treatment placebo comparison
#'
#' @description
#' Compare a post-treatment panel ATT with a pre-treatment placebo ATT on the
#' same units. Propagate the joint estimation uncertainty of both comparisons.
#' @param model A post-treatment, equal-weight panel adapter object returned by
#'   [dml_from_did()] or [dml_from_drdid()].
#' @param pre_model A panel adapter object for a genuinely untreated comparison.
#'   It must use the same unit IDs, treatment-group indicator, estimation method,
#'   and (for cross-fitting) partitions and repetitions as `model`. Both objects
#'   must have explicitly supplied IDs; row numbers generated for wide input are
#'   not sufficient. Pre-period outcomes and covariates may differ.
#' @param k A nonempty vector of finite, nonnegative extrapolation multipliers.
#' @param type `"bias"` transports the bias magnitude: `theta +/- k * abs(theta.pre)`.
#'   `"scaled"` transports the bias factors: multiply this magnitude by `S / S.pre`.
#' @param level Confidence level for each one-sided endpoint, strictly between
#'   0.5 and 1. Default 0.95. Use 0.975 for a region covering both population bound
#'   endpoints with asymptotic probability at least 0.95 by the union bound.
#' @param combine.method How repeated cross-fitting estimates are combined:
#'   `"median"` (default) or `"mean"`, using the existing package rule. Failed
#'   repetitions are never silently dropped.
#' @details
#' Implements the two restrictions in Section 4.2 and the influence functions in
#' Appendix C.2 of Wang et al. (2026). The pre-period causal effect must be zero,
#' including no anticipation. Transportability across periods is a user-imposed
#' restriction, not an implication of a pre-trend test. Independent panel units,
#' overlap and the nuisance conditions of each adapter fit are required.
#'
#' All IFs are aligned by explicit unit ID before forming cross-period covariances.
#' The current implementation requires identical units and group membership; it
#' does not silently intersect samples or treat the two estimates as independent.
#' Inference is cellwise, for each fixed `k`, not simultaneous over cells or `k`.
#'
#' For positive `k`, the usual delta method requires the population placebo ATT
#' to be nonzero. This is not a testable guarantee and the approximation is not
#' uniform near zero. An exactly zero estimated placebo yields point bounds but
#' `NA` standard errors/confidence bounds and a warning. `k = 0` always recovers
#' the original post-treatment analytic IF interval. The scaled method also
#' requires strictly positive finite pre/post sensitivity scales.
#'
#' This low-level function cannot verify dates in arbitrary wide-panel objects.
#' The caller is responsible for selecting genuinely pre-treatment outcomes.
#' Use [did_pretrend_bounds()] to construct and check the comparison from a did fit.
#' @return A `dml_pretrend_bounds` object with `bounds` (one row per `k`),
#'   `results` (per-repetition bounds and unit-level IFs), `joint.covariance`
#'   (per-repetition covariance of the six pre/post estimates), aligned `unit_id`,
#'   and the two adapter objects `model` and `pre_model`. The covariance includes
#'   the `1/n` factor for estimates, not just the covariance of unit IFs.
#' @references
#' Wang, J., Sant'Anna, P. H. C., Chernozhukov, V., and Cinelli, C. (2026).
#' *Omitted Variable Bias in Difference-in-Differences Designs*, Sections 4.2,
#' 5.2 and Appendix C.2.
#' @seealso [did_pretrend_bounds()], [benchmark_bounds()]
#' @md
#' @export
pretrend_bounds <- function(model, pre_model, k = 1,
                            type = c("bias", "scaled"), level = 0.95,
                            combine.method = c("median", "mean")) {
  type <- match.arg(type)
  combine.method <- match.arg(combine.method)
  .pretrend_options(k, level)
  pair <- .pretrend_pair(model, pre_model)
  combine <- if (combine.method == "mean") combine.mean else combine.median
  average <- if (combine.method == "mean") mean else stats::median
  covariance <- lapply(pair$joint, function(x) {
    centered <- sweep(x, 2L, colMeans(x), "-")
    crossprod(centered) / nrow(x)^2
  })
  results <- lapply(k, function(multiplier) {
    lapply(seq_along(pair$joint), function(r)
      .pretrend_one(pair$estimates[[r]], pair$joint[[r]], multiplier, type))
  })
  zero <- FALSE
  rows <- lapply(seq_along(k), function(j) {
    rr <- results[[j]]
    get <- function(key) vapply(rr, function(r) r[[key]], numeric(1))
    lo <- combine(get("theta.minus"), get("se.minus"), na.rm = FALSE)
    hi <- combine(get("theta.plus"), get("se.plus"), na.rm = FALSE)
    regular <- all(vapply(rr, function(r) r$regular, logical(1)))
    if (!regular) zero <<- TRUE
    data.frame(k = k[j], theta.s = average(get("theta.s")),
      theta.pre = average(get("theta.pre")),
      scale.ratio = average(get("scale.ratio")),
      bias.bound = average(get("bias.bound")),
      theta.minus = unname(lo["estimate"]), theta.plus = unname(hi["estimate"]),
      se.minus = unname(lo["se"]), se.plus = unname(hi["se"]),
      lwr = unname(lo["estimate"] - stats::qnorm(level) * lo["se"]),
      upr = unname(hi["estimate"] + stats::qnorm(level) * hi["se"]),
      status = if (regular) "ok" else "zero_pretrend", stringsAsFactors = FALSE)
  })
  if (zero) warning("The estimated pre-trend is exactly zero. For positive k the absolute value is non-differentiable; returning point bounds but NA delta-method SEs and confidence bounds.", call. = FALSE)
  out <- list(bounds = do.call(rbind, rows), type = type, level = level,
    combine.method = combine.method, method = pair$method,
    unit_id = pair$ids, results = results, joint.covariance = covariance,
    model = model, pre_model = pre_model, call = match.call(),
    assumptions = c("zero pre-period causal effect including no anticipation",
      "user-imposed transport restriction for each fixed k",
      "independent equal-weight panel units and valid nuisance fits in both periods",
      "nonzero population placebo ATT for delta inference at positive k"))
  class(out) <- "dml_pretrend_bounds"
  out
}

.pretrend_options <- function(k, level) {
  if (!is.numeric(k) || !length(k) || any(!is.finite(k)) || any(k < 0))
    .drdid_stop("k must be a nonempty vector of finite nonnegative numbers.")
  if (!is.numeric(level) || length(level) != 1L || !is.finite(level) ||
      level <= 0.5 || level >= 1)
    .drdid_stop("level must be a single number strictly between 0.5 and 1 (one-sided endpoint coverage).")
}

.pretrend_pair <- function(model, pre_model) {
  fits <- list(model, pre_model)
  for (m in fits) {
    if (!(inherits(m, "dml_did") || inherits(m, "dml_drdid")) ||
        !identical(m$info$target, "att") || !isTRUE(m$info$conditional) ||
        !identical(m$info$centered, FALSE))
      .drdid_stop("Both models must be equal-weight, uncentered conditional ATT panel adapter objects.")
    s <- m$adapter_data
    n <- length(s$D)
    if (!isTRUE(s$id.supplied) || length(s$unit_id) != n || n < 4L ||
        anyNA(s$unit_id) || anyDuplicated(s$unit_id))
      .drdid_stop("Both fits need explicit unique unit IDs; supply unit_id to the wide DRDID adapter.")
    .drdid_check_weights(s$i.weights, n)
  }
  ids <- model$adapter_data$unit_id
  pre.ids <- pre_model$adapter_data$unit_id
  if (length(ids) != length(pre.ids) || !setequal(ids, pre.ids))
    .drdid_stop("Pre/post fits must use exactly the same units; refit on a common sample explicitly.")
  ix <- match(ids, pre.ids)
  if (!identical(as.numeric(model$adapter_data$D),
                 as.numeric(pre_model$adapter_data$D[ix])))
    .drdid_stop("Treatment-group membership differs between the pre/post fits.")
  method <- function(m) if (is.null(m$refit_spec$method)) "crossfit" else m$refit_spec$method
  modes <- vapply(fits, method, character(1))
  if (!all(modes %in% c("crossfit", "parametric")) || modes[1L] != modes[2L])
    .drdid_stop("Pre/post fits must use the same method (crossfit or parametric).")
  post <- model$results$main$treat
  pre <- pre_model$results$main$treat
  R <- length(post)
  if (!R || length(pre) != R || (modes[1L] == "parametric" && R != 1L))
    .drdid_stop("Pre/post fits must contain the same number of estimation repetitions.")
  if (modes[1L] == "crossfit") {
    folds <- function(m) {
      f <- m$adapter_info$fold_id
      if (is.null(dim(f))) f <- matrix(f, ncol = 1L)
      if (!identical(dim(f), c(length(ids), R)) || anyNA(f))
        .drdid_stop("Cross-fitting partitions are missing or malformed.")
      f
    }
    if (!isTRUE(all.equal(unname(folds(model)), unname(folds(pre_model)[ix, , drop = FALSE]),
                          check.attributes = FALSE)))
      .drdid_stop("Pre/post cross-fitting partitions differ; use matched folds and repetition order.")
  }
  keys <- c("theta.s", "sigma2.s", "nu2.s")
  labels <- c("theta.post", "sigma2.post", "nu2.post", "theta.pre", "sigma2.pre", "nu2.pre")
  estimates <- joint <- vector("list", R)
  for (r in seq_len(R)) {
    es <- lapply(list(post[[r]], pre[[r]]), function(z) vapply(keys, function(key) {
      v <- z$estimates[[key]]
      if (!is.numeric(v) || length(v) != 1L || !is.finite(v))
        .drdid_stop("Both fits must contain finite scalar ATT and sensitivity moments.")
      as.numeric(v)
    }, numeric(1)))
    ps <- lapply(list(post[[r]], pre[[r]]), function(z) vapply(keys, function(key) {
      v <- z$psis[[paste0("psi.", key)]]
      if (!is.numeric(v) || length(v) != length(ids) || any(!is.finite(v)))
        .drdid_stop("Both fits must contain finite aligned unit-level influence functions.")
      as.numeric(v)
    }, numeric(length(ids))))
    estimates[[r]] <- stats::setNames(c(es[[1L]], es[[2L]]), labels)
    joint[[r]] <- cbind(ps[[1L]], ps[[2L]][ix, , drop = FALSE])
    colnames(joint[[r]]) <- labels
  }
  list(ids = ids, method = modes[1L], estimates = estimates, joint = joint)
}

# Algebraic kernel; the IFs supplied here already include nuisance estimation.
.pretrend_one <- function(est, joint, k, type) {
  theta <- unname(est[1L]); placebo <- unname(est[4L])
  ratio <- 1
  psi.ratio <- rep(0, nrow(joint))
  if (type == "scaled") {
    if (any(est[c(2L, 3L, 5L, 6L)] <= 0))
      .drdid_stop("Scaled pre-trend bounds require strictly positive pre/post sensitivity moments.")
    S <- sqrt(est[2L]) * sqrt(est[3L])
    Spre <- sqrt(est[5L]) * sqrt(est[6L])
    ratio <- unname(S / Spre)
    # Relative IFs for S avoid unnecessary products of large scale components.
    psi.ratio <- ratio / 2 * (joint[, 2L] / est[2L] + joint[, 3L] / est[3L] -
                              joint[, 5L] / est[5L] - joint[, 6L] / est[6L])
    if (!is.finite(ratio) || ratio <= 0 || any(!is.finite(psi.ratio)))
      .drdid_stop("The pre/post sensitivity scale ratio or its IF is not finite and positive.")
  }
  width <- k * ratio * abs(placebo)
  if (!is.finite(width)) .drdid_stop("The extrapolated bias bound is not finite.")
  regular <- k == 0 || placebo != 0
  psi.width <- if (k == 0) rep(0, nrow(joint)) else if (!regular)
    rep(NA_real_, nrow(joint)) else
    k * (ratio * sign(placebo) * joint[, 4L] + abs(placebo) * psi.ratio)
  minus <- joint[, 1L] - psi.width
  plus <- joint[, 1L] + psi.width
  se <- function(x) sqrt(sum((x - mean(x))^2)) / length(x)
  list(theta.s = theta, theta.pre = placebo, scale.ratio = ratio,
    bias.bound = width, theta.minus = theta - width, theta.plus = theta + width,
    se.minus = se(minus), se.plus = se(plus), regular = regular,
    psis = list(psi.bias.bound = psi.width, psi.theta.minus = minus,
                psi.theta.plus = plus))
}

#' Construct a pre-trend benchmark for an existing did cell
#'
#' @description
#' Fit a pre-treatment placebo comparison on the frozen post-treatment cell
#' sample, then call [pretrend_bounds()] with joint uncertainty propagation.
#' @param model An equal-weight `dml_did` object from [dml_from_did()], retaining
#'   its original `MP` object and saved panel data. Benchmark-refit objects whose
#'   external fit is no longer an `MP` object are not accepted.
#' @param pre_period Numeric vector `c(base, end)` selecting two distinct saved
#'   periods in chronological order. Both must be strictly before treatment minus
#'   anticipation, and the end cannot exceed the post cell's baseline.
#' @inheritParams pretrend_bounds
#' @details
#' The original post-treatment ATT and IF are preserved. The placebo is a new
#' DRDID fit using exactly the post cell's units and group membership, with the
#' saved covariate design at the pre-period baseline. This is deliberate: with
#' not-yet-treated controls, a pre-treatment cell stored by `att_gt()` can have a
#' different control group. With a universal baseline it can be normalized to
#' zero or have the opposite time orientation. No such stored cell is substituted
#' for the explicit forward comparison requested here.
#'
#' The method and DML settings are inherited from `model`; repeated cross-fitting
#' uses matching partitions. Period lengths are not automatically rescaled: `k`
#' concerns the changes over the selected intervals. This function handles one
#' pre-period contrast at a time, not maxima over several pre-trends or event-study
#' aggregation. It can be used as `FUN` in [did_cell_apply()].
#' @return A `dml_pretrend_bounds` object, additionally recording `pre_period`,
#'   `post_cell` and `sample_policy = "post_cell"`. Its `pre_model` can be reused
#'   in [pretrend_bounds()] to change the type or multipliers without refitting.
#' @references
#' Wang, J., Sant'Anna, P. H. C., Chernozhukov, V., and Cinelli, C. (2026).
#' *Omitted Variable Bias in Difference-in-Differences Designs*, Section 4.2 and
#' Appendix C.2.
#' @examples
#' \dontrun{
#' m <- dml_from_did(fit, group = 2007, time = 2007, method = "parametric")
#' b <- did_pretrend_bounds(m, pre_period = c(2005, 2006), k = c(0, 1, 2))
#' b$bounds
#' pretrend_bounds(m, b$pre_model, k = c(0, 1, 2), type = "scaled")
#' }
#' @md
#' @export
did_pretrend_bounds <- function(model, pre_period, k = 1,
                                type = c("bias", "scaled"), level = 0.95,
                                combine.method = c("median", "mean")) {
  type <- match.arg(type)
  combine.method <- match.arg(combine.method)
  .pretrend_options(k, level)
  if (!inherits(model, "dml_did") || !inherits(model$external_fit, "MP"))
    .drdid_stop("model must be a did cell retaining its original MP object; use pretrend_bounds() for an explicit pair of panel fits.")
  if (!is.numeric(pre_period) || length(pre_period) != 2L ||
      any(!is.finite(pre_period)) || pre_period[1L] >= pre_period[2L])
    .drdid_stop("pre_period must be c(base, end) with two finite, strictly increasing periods.")
  ctx <- .did_adapter_context(model$external_fit)
  if (ctx$unweighted.refit)
    .drdid_stop("Pre-trend benchmarking currently requires equal sampling weights.")
  cell <- model$adapter_info$cell
  if (any(pre_period >= cell$group - ctx$dp$anticipation) ||
      pre_period[2L] > cell$base_period)
    .drdid_stop("Both pre_period dates must precede treatment minus anticipation, and end no later than the post baseline.")
  periods <- match(pre_period, ctx$times)
  if (anyNA(periods)) .drdid_stop("A requested pre_period is absent from the saved panel.")
  # Replay the original post cell, but never refit its DML auxiliary estimator.
  src <- .did_cell_sample(ctx, cell$group, cell$time)
  close <- function(a, b) isTRUE(all.equal(a, b, tolerance = 1e-10, check.attributes = FALSE))
  fields <- c("y1", "y0", "D", "X", "i.weights", "unit_id")
  if (!all(vapply(fields, function(key) close(src$sample[[key]], model$adapter_data[[key]]), logical(1))) ||
      !close(as.numeric(src$theta), model$results$main$treat[[1L]]$estimates$theta.s) ||
      !close(src$phi, model$results$main$treat[[1L]]$psis$psi.theta.s))
    .drdid_stop("The saved post cell no longer matches its source data or ATT/IF.")
  sample <- src$sample
  take <- match(sample$unit_id, ctx$ids)
  y <- ctx$data[[ctx$dp$yname]]
  sample$y0 <- y[ctx$rows[[periods[1L]]][take]]
  sample$y1 <- y[ctx$rows[[periods[2L]]][take]]
  sample$X <- ctx$X[[periods[1L]]][take, , drop = FALSE]
  source <- .drdid_refit_source(sample, model$refit_spec$source)
  method <- model$refit_spec$method
  if (is.null(method)) method <- "crossfit"
  pre_model <- dml_from_drdid(source, y1 = sample$y1, y0 = sample$y0,
    D = sample$D, covariates = sample$X, i.weights = sample$i.weights,
    unit_id = sample$unit_id, dml_args = model$refit_spec$dml_args, method = method)
  pre_model$adapter_info$pre_period <- pre_period
  pre_model$adapter_info$sample_policy <- "post_cell"
  out <- pretrend_bounds(model, pre_model, k = k, type = type, level = level,
                         combine.method = combine.method)
  out$pre_period <- pre_period
  out$post_cell <- cell[c("group", "time", "base_period")]
  out$sample_policy <- "post_cell"
  out$call <- match.call()
  out
}

#' @param x A `dml_pretrend_bounds` object.
#' @param digits Number of significant digits for printing.
#' @param ... Further arguments passed to `print.data.frame()`.
#' @rdname pretrend_bounds
#' @export
print.dml_pretrend_bounds <- function(x, digits = 4, ...) {
  cat("Pre-trend benchmark:", x$type, "extrapolation;", x$method, "fits\n")
  if (!is.null(x$pre_period))
    cat("Pre-period comparison:", x$pre_period[1L], "to", x$pre_period[2L],
        "(post-cell units and groups)\n")
  cat("One-sided endpoint level:", format(x$level),
      "; pointwise in k and cell; benchmark uncertainty propagated.\n\n")
  print(x$bounds, digits = digits, row.names = FALSE, ...)
  invisible(x)
}
