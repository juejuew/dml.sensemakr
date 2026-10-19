#' Adapt a did group-time ATT for sensitivity analysis
#'
#' Uses the estimation data saved by `did::att_gt()`. With equal sampling weights,
#' preserves the selected ATT and IF and estimates sensitivity moments using
#' cross-fitting (default) or full-sample parametric fitting. With unequal weights,
#' the default method refits the entire cell using unweighted DML.
#'
#' @param fit An `MP` object from `did::att_gt()`, estimated with `est_method = "dr"`.
#'   The equal-weight hybrid path requires saved influence functions
#'   (`compute_inffunc = TRUE`) and compatible DRDID panel interfaces.
#' @param group,time Numeric cohort and outcome-period labels of one estimated
#'   post-treatment cell (`time >= group`).
#' @param dml_args Named list of DML settings; see [dml_from_drdid()].
#'   Must be empty for `method = "parametric"`.
#' @param method `"crossfit"` (default) retains the existing cross-fitted method.
#'   `"parametric"` uses the source cell's full-sample logistic MLE and control-arm
#'   linear regression, with no sample splitting. Requires equal weights and no
#'   active source trimming/clipping. It never switches to an unweighted refit.
#' @return With equal weights, a `dml_did` object inheriting from `dml`. With
#'   unequal weights, an ordinary `dml` object with native DML refitting behavior.
#'   Both contain `adapter_data` (the frozen two-period panel), `external_fit`
#'   (the original MP object), and `adapter_info$cell` (cell labels, base period
#'   and sample sizes). `adapter_info$mode` is `"hybrid"`, `"parametric"` or
#'   `"unweighted_dml"`. Source IF row indices and scale are recorded for the
#'   two source-preserving modes. Parametric results include `parametric`
#'   (coefficients, gradients, assumptions) and each result's `nuisance.correction`.
#' @details Supports the balanced panel estimator path, either control group,
#'   either base-period setting, anticipation, and both `faster_mode` settings.
#'   Hybrid compatibility is checked from the saved data, settings, design and labeled
#'   IF, and by reproducing the cell ATT and ordered IF. Installed package versions
#'   are diagnostic metadata, not acceptance gates or evidence of the producer
#'   version. Tested with did 2.5.1 and 2.5.1.902 and DRDID 1.3.0. Missing required
#'   information or incompatible runtime interfaces still causes an error.
#'   Sampling weights must be finite and positive. With the default method, if unequal weights occur in
#'   the retained data, the adapter warns and treats all selected units equally.
#'   It calls [dml()] to estimate the ATT, its IF and all sensitivity components
#'   anew, without using or requiring the original ATT/IF. Numerical results
#'   match direct unweighted `dml()` on the same ordered cell data with the same
#'   resolved DML settings and seed; the settings are in `adapter_info$dml_args`.
#'   Original cell weights are in `adapter_info$source_weights`. This changes
#'   the weighting of the estimand and is not weighted sensitivity analysis.
#'   The refit uses the processed sample saved by did and cannot restore rows
#'   removed by its original cleaning.
#'
#'   Active source trimming or DML clipping is rejected on the hybrid path.
#'   The unweighted refit follows native `dml()` behavior, including `ps.trim`.
#'   Coarser clustering, repeated cross-sections, `fix_weights = "varying"`,
#'   custom estimators, pre-treatment cells and `aggte()` remain unsupported.
#'
#'   For the hybrid path, the source stores `(N / n_cell) * IF_cell`. Row names align units; multiplying
#'   by `n_cell / N` recovers the cell IF. A DRDID reconstruction checks both the
#'   original ATT and this ordered IF. Fast mode uses the saved design tensors;
#'   slow mode reevaluates the formula on all saved periods before subsetting.
#'   The slow formula's environment and contrast options must remain compatible
#'   with the source fit. Do not modify the source object's data or tensors.
#'
#'   Inference uses analytic joint influence functions for this cell; original
#'   bootstrap SEs and simultaneous critical values are not reused. Validity
#'   requires the identifying and nuisance conditions of the estimators. These
#'   checks establish computational consistency, not model correctness.
#'
#'   The parametric method uses DRDID's actual logistic-MLE/OLS solvers and the
#'   saved cell design. It retains the source ATT/IF and computes new sensitivity
#'   moments: these moments are not quantities originally returned by did.
#'   Joint inference includes the full-sample nuisance estimation effects via
#'   estimating-equation derivatives. Its usual sensitivity interpretation
#'   requires correct nuisance functions, regular fixed-dimensional models,
#'   overlap, finite moments and the original identifying conditions. Sandwich
#'   corrections do not repair a misspecified sensitivity target. Each reduced
#'   benchmark model needs its own nuisance validity conditions.
#'   Parametric fits are deterministic and store `info$cf.folds = 0L`, one result,
#'   empty `fold_id`/`rep_seeds`, and no auxiliary DML estimates. Benchmarking
#'   refits the same logistic/OLS procedure on the frozen cell and rejects learner
#'   overrides. See `inst/notes/parametric-inference.md` in the source package.
#'
#'   [dml_benchmark()] refits DRDID and DML for hybrid fits, and DML alone for
#'   unweighted refits, after removing design columns on the same cell and folds.
#'   Parametric fits refit the source estimator and full-sample sensitivity moments.
#'   [did_cell_apply()] applies sensitivity or benchmark
#'   functions across a batch; its intervals remain cellwise, not simultaneous.
#'   `dml_args$cf.reps` enables repeated cross-fitting (default 1). For hybrid fits, each repetition
#'   uses the same source ATT/IF and its own sensitivity moments and joint IF
#'   covariance; downstream methods retain their mean/median aggregation.
#'   A failure in any hybrid repetition rejects that cell. Its `adapter_info$fold_id` is a
#'   vector for one repetition or a unit-by-repetition matrix otherwise;
#'   `rep_seeds` and `dml.att` contain one value per repetition.
#' @seealso [dml_from_did_cells()], [did_cell_apply()], [dml_from_drdid()]
#' @md
#' @export
dml_from_did <- function(fit, group, time, dml_args = list(),
                          method = c("crossfit", "parametric")) {
  method <- match.arg(method)
  if (method == "parametric") .panel_parametric_args(dml_args)
  ctx <- .did_adapter_context(fit)
  out <- .did_fit_cell(ctx, group, time, dml_args, method)
  if (inherits(out, "dml_did")) out$call <- match.call()
  else out$adapter_call <- match.call()
  out
}

.did_fit_cell <- function(ctx, group, time, dml_args, method = "crossfit") {
  if (method == "parametric" && ctx$unweighted.refit)
    .drdid_stop("method = 'parametric' requires equal sampling weights; it does not silently switch to an unweighted DML refit.")
  src <- .did_cell_sample(ctx, group, time)
  info <- list(source = "did", mode = "hybrid", reconstruction.version = ctx$version,
                DRDID.version = ctx$drdid.version, estMethod = "dr (DRDID trad)",
                compatibility = list(basis = "saved structure and ATT/IF replay", layout = ctx$layout),
                inference = "analytic joint IF", cell = src$cell,
                control_group = ctx$dp$control_group,
                anticipation = ctx$dp$anticipation, base_period = ctx$dp$base_period,
                faster_mode = ctx$dp$faster_mode, fix_weights = ctx$dp$fix_weights,
                original.se = ctx$fit$se[src$cell$column],
                original.boot = ctx$dp$bstrap, original.cband = ctx$dp$cband)
  if (ctx$unweighted.refit) {
    warning("Unequal sampling weights detected: refitting this cell with equal weights using DML; the original ATT and IF are not used.",
            call. = FALSE)
    sample <- src$sample
    args <- .drdid_dml_args(dml_args, length(sample$D))
    # Keep dml's own estimates, IFs, class, call and call environment. Native
    # benchmarking will therefore refit DML rather than the source estimator.
    out <- do.call(dml, c(list(y = sample$y1 - sample$y0, d = sample$D, x = sample$X), args))
    info$mode <- "unweighted_dml"
    info$estMethod <- "DML (equal-weight refit)"
    info$compatibility$basis <- "saved data and design"
    info$inference <- "analytic DML IF"
    info$source_weights <- sample$i.weights
    info$dml_args <- args
    sample$i.weights <- rep(1, length(sample$D))
    out$adapter_data <- sample
    out$adapter_info <- info
    out$external_fit <- ctx$fit
    return(out)
  }
  .did_assemble(src$sample, src$theta, src$phi, src$diagnostics, info,
                ctx$fit, dml_args, method)
}

.did_assemble <- function(sample, theta, phi, diagnostics, info, external_fit,
                           dml_args, method = "crossfit") {
  source <- list(method = "trad", trim.level = 0.995)
  hybrid <- if (method == "parametric")
    .panel_parametric_aux(sample, theta, phi, source, diagnostics) else
    .panel_hybrid_aux(sample, theta, phi, dml_args)
  out <- hybrid$fit
  out$external_fit <- external_fit
  out$adapter_data <- sample
  out$adapter_info <- c(info, list(unit_id = sample$unit_id, fold_id = hybrid$fold,
    rep_seeds = hybrid$rep_seeds,
    diagnostics = diagnostics))
  if (method == "parametric") {
    out$adapter_info$mode <- "parametric"
    out$adapter_info$inference <- out$parametric$inference
  } else {
    out$adapter_info$dml.att <- vapply(out$auxiliary_dml_results$main$treat,
                                      function(r) r$estimates$theta.s, numeric(1))
  }
  out$refit_spec <- list(source = source, dml_args = hybrid$args, method = method)
  out$call <- NULL
  out$call.env <- NULL
  class(out) <- c("dml_did", "dml")
  out
}

.did_benchmark_refit <- function(model, x, dreg = NULL, yreg = NULL) {
  if (!ncol(x)) .drdid_stop("A benchmark must retain at least one design column; retain the intercept for an unconditional refit.")
  method <- model$refit_spec$method
  if (is.null(method)) method <- "crossfit"
  if (method == "parametric" && (!is.null(dreg) || !is.null(yreg)))
    .drdid_stop("Parametric benchmarks use the source logistic/OLS models; dreg/yreg overrides are not supported.")
  sample <- model$adapter_data
  sample$X <- x
  source <- .drdid_refit_source(sample, model$refit_spec$source)
  spec <- c(model$refit_spec$source, list(phi = as.numeric(source$att.inf.func)))
  diagnostics <- .panel_check_source(source$ATT, sample, spec)
  args <- model$refit_spec$dml_args
  if (!is.null(dreg)) args$dreg <- dreg
  if (!is.null(yreg)) args$yreg <- yreg
  info <- model$adapter_info
  info[c("unit_id", "fold_id", "rep_seeds", "diagnostics", "dml.att")] <- NULL
  info$original.se <- source$se
  info$original.boot <- FALSE
  info$original.cband <- FALSE
  info$benchmark.refit <- TRUE
  # The reduced estimator is a fresh DRDID fit on the frozen cell, not an MP.
  out <- .did_assemble(sample, source$ATT, spec$phi, diagnostics, info, source, args, method)
  .panel_check_benchmark_alignment(model, out)
  out
}

#' Convert several did ATT cells without dropping failures
#'
#' @inheritParams dml_from_did
#' @param cells A data.frame with numeric `group` and `time` columns and unique
#'   rows. NULL selects all post-treatment cells in source order, including cells
#'   whose estimation failed. Pre-treatment cells are listed in `excluded`.
#' @param on_error `"record"` keeps an error message and a NULL entry for each
#'   failed cell; `"stop"` stops at the first failure. Invalid source objects or
#'   invalid selection syntax always stop before fitting.
#' @return A `dml_did_cells` object with `models` (one entry per selected cell),
#'   `status` (cell labels, success/error status, messages, warnings and sample
#'   sizes), and `excluded` (unselected source cells). Names remain in source or
#'   explicit selection order. No ATT aggregation is performed.
#' @md
#' @export
dml_from_did_cells <- function(fit, cells = NULL, dml_args = list(),
                                on_error = c("record", "stop"),
                                method = c("crossfit", "parametric")) {
  on_error <- match.arg(on_error)
  method <- match.arg(method)
  if (method == "parametric") .panel_parametric_args(dml_args)
  ctx <- .did_adapter_context(fit)
  all.cells <- data.frame(group = fit$group, time = fit$t)
  default <- is.null(cells)
  if (default) cells <- all.cells[all.cells$time >= all.cells$group, , drop = FALSE]
  if (!is.data.frame(cells) || !all(c("group", "time") %in% names(cells)))
    .drdid_stop("cells must be a data.frame with group and time columns.")
  cells <- as.data.frame(cells[c("group", "time")])
  if (!nrow(cells) || !is.numeric(cells$group) || !is.numeric(cells$time) ||
      any(!is.finite(c(cells$group, cells$time))) || anyDuplicated(cells))
    .drdid_stop("Select at least one cell with unique, finite numeric group-time labels.")
  # Catch malformed DML options once, independently of per-cell estimation.
  if (method == "crossfit") .drdid_dml_args(dml_args, ctx$N)
  labels <- .did_cell_labels(cells)
  status <- .did_status_table(cells, labels)
  models <- stats::setNames(vector("list", nrow(cells)), labels)
  for (i in seq_len(nrow(cells))) {
    ans <- .did_capture(function() .did_fit_cell(ctx, cells$group[i], cells$time[i], dml_args, method))
    if (!is.null(ans$error) && on_error == "stop")
      .drdid_stop(labels[i], ": ", ans$error)
    models[i] <- list(ans$value)
    status <- .did_record_status(status, i, ans)
    if (is.null(ans$error)) {
      status$n_cell[i] <- ans$value$adapter_info$cell$n_cell
      status$base_period[i] <- ans$value$adapter_info$cell$base_period
    }
  }
  selected <- vapply(seq_len(nrow(all.cells)), function(i)
    any(cells$group == all.cells$group[i] & cells$time == all.cells$time[i]), logical(1))
  excluded <- all.cells[!selected, , drop = FALSE]
  excluded$reason <- ifelse(default & excluded$time < excluded$group,
                            "pre-treatment or normalized reference cell", "not selected")
  structure(list(models = models, status = status, excluded = excluded,
                  call = match.call()), class = "dml_did_cells")
}

#' Apply a function to each successfully adapted did cell
#'
#' @param object A result from [dml_from_did_cells()].
#' @param FUN Function taking a `dml` fit as its first argument, for example
#'   [sensemakr()], [confidence_bounds()] or [dml_benchmark()].
#' @param ... Additional arguments passed to `FUN`.
#' @inheritParams dml_from_did_cells
#' @return A `did_cell_results` object with `values`, `status` and `excluded`.
#'   Every selected cell retains an entry, including upstream and new failures.
#'   A function returning NULL successfully is distinguished by status `"ok"`.
#' @md
#' @export
did_cell_apply <- function(object, FUN, ..., on_error = c("record", "stop")) {
  if (!inherits(object, "dml_did_cells")) .drdid_stop("object must come from dml_from_did_cells().")
  FUN <- match.fun(FUN)
  on_error <- match.arg(on_error)
  status <- object$status
  values <- stats::setNames(vector("list", nrow(status)), status$cell)
  for (i in seq_len(nrow(status))) {
    if (status$status[i] != "ok") {
      if (on_error == "stop") .drdid_stop(status$cell[i], ": ", status$message[i])
      next
    }
    ans <- .did_capture(function() FUN(object$models[[i]], ...))
    if (!is.null(ans$error) && on_error == "stop") .drdid_stop(status$cell[i], ": ", ans$error)
    values[i] <- list(ans$value)
    status <- .did_record_status(status, i, ans)
  }
  structure(list(values = values, status = status, excluded = object$excluded),
             class = "did_cell_results")
}

.did_cell_labels <- function(cells) {
  paste0("g=", sprintf("%.17g", cells$group), ",t=", sprintf("%.17g", cells$time))
}

.did_status_table <- function(cells, labels) {
  data.frame(cell = labels, group = cells$group, time = cells$time,
             status = "pending", message = "", warnings = "",
             n_cell = NA_integer_, base_period = NA_real_, stringsAsFactors = FALSE)
}

.did_capture <- function(fun) {
  warnings <- character(); error <- NULL
  value <- tryCatch(withCallingHandlers(fun(), warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }), error = function(e) { error <<- conditionMessage(e); NULL })
  list(value = value, error = error, warnings = unique(warnings))
}

.did_record_status <- function(status, i, ans) {
  status$status[i] <- if (is.null(ans$error)) "ok" else "error"
  status$message[i] <- if (is.null(ans$error)) "" else ans$error
  status$warnings[i] <- paste(Filter(nzchar, c(status$warnings[i], ans$warnings)), collapse = " | ")
  status
}

#' @export
summary.dml_did <- function(object, combine.method = "median", ...) {
  out <- summary.dml(object, combine.method = combine.method, ...)
  out$adapter_info <- object$adapter_info
  class(out) <- c("summary_dml_did", "summary_dml")
  out
}

#' @export
print.summary_dml_did <- function(x, ...) {
  cell <- x$adapter_info$cell
  parametric <- identical(x$adapter_info$mode, "parametric")
  cat(if (parametric) "Parametric did ATT(" else "Hybrid did ATT(",
      cell$group, ", ", cell$time, ") sensitivity fit\n", sep = "")
  if (parametric) cat("Full-sample logistic/OLS; no cross-fitting\n")
  cat("Base period: ", cell$base_period, "; cell units: ", cell$n_cell,
      "; inference: analytic joint IF\n\n", sep = "")
  print(x$main)
  invisible(x)
}

#' @export
print.dml_did <- function(x, ...) {
  print(summary(x, ...))
  invisible(x)
}

#' @export
print.dml_did_cells <- function(x, ...) {
  print(x$status, row.names = FALSE)
  if (nrow(x$excluded)) cat(nrow(x$excluded), "source cells excluded; see $excluded.\n")
  invisible(x)
}

#' @export
print.did_cell_results <- function(x, ...) {
  print(x$status, row.names = FALSE)
  invisible(x)
}
