#' Adapt a panel DRDID estimate for conditional ATT sensitivity analysis
#'
#' Retains the external ATT and influence function and estimates the residual
#' variance and conditional imbalance using cross-fitted DML. The source fit is
#' replayed to check its consistency with the supplied, ordered estimation data.
#'
#' @param fit A DRDID panel DR result estimated with `inffunc = TRUE`.
#' @param data Optional balanced, complete long panel data.frame.
#' @param yname,tname,idname,dname Column names for long input. Literal names can
#'   be recovered from the matched source call; unresolved expressions cannot.
#' @param xformla One-sided covariate formula for long input.
#' @param weightsname Optional sampling-weight column in long input. Raw values
#'   must be constant and positive across units and periods, before normalization.
#' @param y1,y0,D Post/pre outcomes and treatment indicator for wide input,
#'   in the exact source influence-function order.
#' @param covariates Numeric wide-input design matrix, including any intercept
#'   used by the source estimator. NULL means intercept only.
#' @param i.weights Optional constant positive sampling weights for wide input.
#' @param unit_id Optional unique unit IDs in wide-input order. They label the
#'   supplied order and do not independently establish the original IF mapping.
#' @param dml_args Named list of DML options, such as learners, `cf.folds`,
#'   `cf.reps`, `cf.seed` and `ps.trim`. Repetitions default to 1. The estimand
#'   and clean tuning are fixed. Every repetition must pass the overlap checks;
#'   active probability clipping is not supported.
#' @return A `dml_drdid` object inheriting from `dml`, with the external fit,
#'   ordered panel data, source diagnostics and joint estimator covariance.
#' @details Supports panel `imp` and `trad`, equal weights, successful source
#'   estimation, and no active source or DML propensity trimming/clipping.
#'   Acceptance uses the input structure, the required DRDID function interfaces,
#'   and a numerical ATT/IF replay, not an installed-version equality check.
#'   DRDID 1.3.0 is tested; other versions must satisfy these same checks.
#'   Valid inference additionally requires compatible asymptotic linearity,
#'   nuisance convergence, overlap and moment conditions. Reconstruction checks
#'   cannot prove model correctness or certify the provenance of historical data.
#'   Benchmark covariates are design-matrix columns (or named groups of columns).
#'   Both DRDID and DML are refit on the frozen sample for every benchmark.
#'   With repeated cross-fitting, the external ATT and IF are retained in every
#'   repetition and combined with that repetition's DML moment IFs. Downstream
#'   methods use their existing mean/median aggregation, without dividing the
#'   external ATT variance by the number of repetitions. One failed repetition
#'   rejects the fit. `adapter_info$fold_id` is a vector for one repetition or a
#'   unit-by-repetition matrix otherwise; `rep_seeds` records the DML seeds and
#'   `dml.att` contains the auxiliary ATT for every repetition.
#' @export
dml_from_drdid <- function(fit, data = NULL, yname = NULL, tname = NULL,
                           idname = NULL, dname = NULL, xformla = NULL,
                           weightsname = NULL, y1 = NULL, y0 = NULL, D = NULL,
                           covariates = NULL, i.weights = NULL, unit_id = NULL,
                           dml_args = list()) {
  version <- .drdid_require()
  spec <- .drdid_spec(fit)
  sample <- .drdid_panel_data(fit, spec, data, yname, tname, idname, dname,
                              xformla, weightsname, y1, y0, D, covariates,
                              i.weights, unit_id)
  diagnostics <- .drdid_check_source(fit, sample, spec)
  hybrid <- .panel_hybrid_aux(sample, fit$ATT, spec$phi, dml_args)
  aux <- hybrid$fit
  args <- hybrid$args
  fold <- hybrid$fold
  aux$external_fit <- fit
  aux$adapter_data <- sample
  aux$adapter_info <- list(source = "DRDID", reconstruction.version = version,
                           compatibility = list(basis = "input structure and ATT/IF replay"),
                           estMethod = spec$method, inference = "analytic joint IF",
                           unit_id = sample$unit_id, fold_id = fold,
                           rep_seeds = hybrid$rep_seeds,
                           diagnostics = diagnostics, original.se = fit$se,
                           original.boot = fit$argu$boot,
                           dml.att = vapply(aux$auxiliary_dml_results$main$treat,
                                            function(r) r$estimates$theta.s, numeric(1)))
  aux$refit_spec <- list(source = list(method = spec$method, trim.level = spec$trim.level),
                         dml_args = args)
  aux$call <- match.call()
  aux$call.env <- NULL
  class(aux) <- c("dml_drdid", "dml")
  aux
}

.drdid_constant_design_learner <- function(reg) {
  if (is.null(reg)) return(NULL)
  if (is.character(reg)) reg <- list(method = reg)
  if (is.list(reg) && identical(tolower(names(reg)), c("yreg0", "yreg1")))
    return(lapply(reg, .drdid_constant_design_learner))
  if (is.null(reg$preProcess)) reg$preProcess <- "ignore"
  reg
}

.drdid_dml_args <- function(args, n) {
  if (!is.list(args) || (length(args) && (is.null(names(args)) ||
      anyNA(names(args)) || any(!nzchar(names(args))) || anyDuplicated(names(args)))))
    .drdid_stop("dml_args must be a uniquely named list.")
  locked <- list(model = "npm", target = "att", conditional = TRUE, centered = FALSE,
                 dirty.tuning = FALSE, y.class = FALSE)
  for (key in intersect(names(args), names(locked))) {
    if (!isTRUE(all.equal(args[[key]], locked[[key]])))
      .drdid_stop("dml_args$", key, " is fixed by the panel adapter.")
  }
  allowed <- setdiff(names(formals(dml)), c("y", "d", "x", "groups"))
  if (any(!names(args) %in% allowed)) .drdid_stop("Unsupported dml_args: ", paste(setdiff(names(args), allowed), collapse = ", "))
  out <- utils::modifyList(list(cf.folds = 5L, cf.reps = 1L, cf.seed = 1L, ps.trim = 0.01,
                                verbose = FALSE), args, keep.null = TRUE)
  out <- utils::modifyList(out, locked)
  k <- out$cf.folds; seed <- out$cf.seed
  if (!is.numeric(k) || length(k) != 1L || !is.finite(k) || k != floor(k) || k < 2 || k > n)
    .drdid_stop("cf.folds must be an integer between 2 and the number of units.")
  # dml() draws repetition seeds without replacement from seq_len(1e6).
  reps <- out$cf.reps
  if (!is.numeric(reps) || length(reps) != 1L || !is.finite(reps) ||
      reps != floor(reps) || reps < 1 || reps > 1e6)
    .drdid_stop("cf.reps must be an integer between 1 and 1000000 (DML's seed-pool size).")
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) || seed != floor(seed) || seed < 0 || seed > .Machine$integer.max)
    .drdid_stop("Provide a finite, nonnegative integer cf.seed for reproducible benchmarking.")
  trim <- out$ps.trim
  if (is.numeric(trim) && length(trim) == 1L) limits <- c(trim, 1 - trim)
  else if (is.list(trim) && all(c("lower", "upper") %in% names(trim)))
    limits <- c(trim$lower, trim$upper)
  else .drdid_stop("ps.trim must be a scalar or list(lower = ..., upper = ...).")
  if (!is.numeric(limits) || length(limits) != 2L || any(!is.finite(limits)) ||
      limits[1L] < 0 || limits[2L] > 1 || limits[1L] >= limits[2L])
    .drdid_stop("ps.trim must define an ordered interval within [0, 1].")
  out
}

.drdid_combine_moments <- function(result, theta, phi, y, D) {
  n <- length(phi)
  a <- result$estimates$sigma2.s; b <- result$estimates$nu2.s
  pa <- result$psis$psi.sigma2.s; pb <- result$psis$psi.nu2.s
  if (length(pa) != n || length(pb) != n || any(!is.finite(c(pa, pb))))
    .drdid_stop("DML moment influence functions are not aligned finite unit-level vectors.")
  joint <- cbind(theta.s = phi, sigma2.s = pa, nu2.s = pb)
  centered <- sweep(joint, 2L, colMeans(joint), "-")
  covariance <- crossprod(centered) / n^2
  pS <- b * pa + a * pb
  se.if <- function(v) sqrt(sum((v - mean(v))^2)) / n
  # Define the mean components consistently with the external ATT, rather than
  # leaving a DML decomposition that no longer sums to the reported estimate.
  theta1 <- mean(y[D == 1])
  phi1 <- D / mean(D) * (y - theta1)
  result$psis$psi.theta.s <- phi
  result$psis$psi.theta.1s <- phi1
  result$psis$psi.theta.0s <- phi1 - phi
  result$psis$psi.S2 <- pS
  values <- list(theta.s = as.numeric(theta), theta.1s = theta1,
                  theta.0s = theta1 - theta, se.theta.s = se.if(phi),
                  se.theta.1s = se.if(phi1), se.theta.0s = se.if(phi1 - phi),
                  se.sigma2.s = se.if(pa), se.nu2.s = se.if(pb),
                  S2 = a * b, se.S2 = se.if(pS),
                  cov.theta.S2 = sum((phi - mean(phi)) * (pS - mean(pS))) / n^2)
  result$estimates[names(values)] <- values
  result$joint.covariance <- covariance
  result
}

.drdid_benchmark_refit <- function(model, x, dreg = NULL, yreg = NULL) {
  sample <- model$adapter_data
  if (!ncol(x)) .drdid_stop("A benchmark must retain at least one design column; retain the source intercept for an unconditional refit.")
  sample$X <- x
  spec <- model$refit_spec$source
  fit <- .drdid_refit_source(sample, spec)
  args <- model$refit_spec$dml_args
  if (!is.null(dreg)) args$dreg <- dreg
  if (!is.null(yreg)) args$yreg <- yreg
  out <- dml_from_drdid(fit, y1 = sample$y1, y0 = sample$y0, D = sample$D,
                 covariates = sample$X, i.weights = sample$i.weights,
                 unit_id = sample$unit_id, dml_args = args)
  .panel_check_benchmark_alignment(model, out)
  out
}

.panel_check_benchmark_alignment <- function(model, refit) {
  if (!identical(model$adapter_info$unit_id, refit$adapter_info$unit_id) ||
      !identical(model$adapter_info$fold_id, refit$adapter_info$fold_id) ||
      length(model$results$main$treat) != length(refit$results$main$treat))
    .drdid_stop("Benchmark units or cross-fitting partitions differ from the saved fit. Restore the original RNG settings and refit specifications.")
  invisible(NULL)
}

#' @export
summary.dml_drdid <- function(object, combine.method = "median", ...) {
  out <- summary.dml(object, combine.method = combine.method, ...)
  out$adapter_info <- object$adapter_info
  class(out) <- c("summary_dml_drdid", "summary_dml")
  out
}

#' @export
print.summary_dml_drdid <- function(x, ...) {
  cat("Hybrid panel DiD sensitivity fit\n")
  cat("ATT: DRDID (", x$adapter_info$estMethod, "); sensitivity moments: cross-fitted DML\n", sep = "")
  cat("Inference: analytic joint influence functions\n\n")
  print(x$main)
  invisible(x)
}

#' @export
print.dml_drdid <- function(x, ...) {
  print(summary(x, ...))
  invisible(x)
}

# Shared auxiliary fit: frontends supply an ordered panel ATT and its IF.
.panel_hybrid_aux <- function(sample, theta, phi, dml_args) {
  args <- .drdid_dml_args(dml_args, length(sample$D))
  # caret's range preprocessing has no scale to estimate on a constant design.
  # Preserve explicit preprocessing choices; omit the default in this case.
  if (all(apply(sample$X, 2L, function(v) length(unique(v)) == 1L))) {
    # dml's yreg/dreg defaults are promises referring to reg, not literal
    # learner specifications. Resolve the shared learner before either override.
    default.reg <- args$reg
    if (is.null(default.reg)) default.reg <- eval(formals(dml)$reg, environment(dml))
    for (key in c("dreg", "yreg")) {
      learner <- args[[key]]
      if (is.null(learner)) learner <- default.reg
      args[[key]] <- .drdid_constant_design_learner(learner)
    }
  }
  # Preflight all exact partitions before fitting, including control-arm scaling.
  # Use the same seed draw and fold assignment as dml()/cross.fitting().
  set.seed(args$cf.seed)
  rep.seeds <- sample.int(1e6, args$cf.reps)
  n <- length(sample$D)
  folds <- vapply(rep.seeds, function(seed) {
    set.seed(seed)
    rep.int(seq_len(args$cf.folds), ceiling(n / args$cf.folds))[sample.int(n)]
  }, integer(n))
  for (i in seq_len(args$cf.reps)) {
    for (k in seq_len(args$cf.folds)) {
      train <- folds[, i] != k
      control <- train & sample$D == 0
      if (length(unique(sample$D[train])) != 2L || sum(control) < 2L ||
          length(unique((sample$y1 - sample$y0)[control])) < 2L)
        .drdid_stop("Repetition ", i, ", training fold ", k,
          " lacks both groups or control-outcome variation; change cf.folds/cf.seed or provide more data.")
    }
  }
  aux <- do.call(dml, c(list(y = sample$y1 - sample$y0, d = sample$D, x = sample$X), args))
  if (length(aux$fits) != args$cf.reps || length(aux$results$main$treat) != args$cf.reps)
    .drdid_stop("DML did not return every requested cross-fitting repetition.")
  aux$auxiliary_dml_results <- aux$results
  results <- lapply(seq_len(args$cf.reps), function(i) {
    result <- aux$auxiliary_dml_results$main$treat[[i]]
    pscore <- aux$fits[[i]]$preds$dhat
    if (any(!is.finite(pscore)) || any(pscore <= 0 | pscore >= 1) ||
        result$trim.summary$trimmed_num$all > 0)
      .drdid_stop("Repetition ", i,
        ": DML propensity clipping is active or predictions are outside (0, 1). Use a suitable learner/overlap specification; no repetitions were silently removed.")
    if (!is.finite(result$estimates$S2) || result$estimates$S2 <= 0)
      .drdid_stop("Repetition ", i,
        ": DML produced non-positive or non-finite S2; regular sensitivity inference is unavailable.")
    tryCatch(.drdid_combine_moments(result, theta, phi, sample$y1 - sample$y0, sample$D),
      error = function(e) .drdid_stop("Repetition ", i, ": ", conditionMessage(e)))
  })
  aux$results$main <- list(treat = results)
  aux$coefs$main <- list(treat = combine.cross.fits(results))
  colnames(folds) <- paste0("rep", seq_len(args$cf.reps))
  # Keep the single-repetition metadata shape compatible with earlier releases.
  fold <- if (args$cf.reps == 1L) unname(folds[, 1L]) else folds
  list(fit = aux, args = args, fold = fold, rep_seeds = rep.seeds)
}
