##' Benchmarks for the strength of latent variables using observed covariates
##' @description
##' Compute benchmarks for the strength of latent variables, under the assumption that the gains in explanatory power due to latent variables is proportional to the gains of observed covariates.
##' @param model an object of class \code{\link{dml}}, fit with a single target (\code{"ate"}, \code{"att"}, or \code{"atu"}).
##' @param benchmark_covariates the observed covariates to benchmark. Either a character vector of column names in the model's \code{x} (each benchmarked on its own), or a named list where each element is a character vector of column names to drop \emph{together} in the leave-one-out refit (e.g. the dummy columns of a factor: \code{list(region = c("region3", "region4"))}). List element names become the row labels; unnamed elements are labelled by the column name (singletons) or the columns joined by \code{"+"}.
##' @returns An object of class \code{dml_benchmark} containing benchmark results. \code{gain.Y} and \code{gain.D} are population non-negative quantities, but their DML estimates are reported as-is (not floored at 0): a small negative value is expected sampling noise when the covariate's true gain is near zero, and reflects that rather than an impossible negative population value.
##' @export
dml_benchmark <- function(model, benchmark_covariates){
  bench <- bench_fun(model = model, benchmark_covariates = benchmark_covariates)
  class(bench) <- "dml_benchmark"
  return(bench)
}

##' Print and summary methods for DML benchmarks
##' @description Print and summary methods for benchmark results.
##' @param x an object of class \code{\link{dml_benchmark}} or \code{summary_dml_benchmark}.
##' @param digits minimal number of significant digits.
##' @rdname summary.dml_benchmark
##' @export
print.dml_benchmark <- function(x, digits = max(3L, getOption("digits") - 3L),
                                 combine.method = "median", ...){
  print(summary(x, combine.method = combine.method), digits = digits, ...)
}

##' @param object an object of class \code{\link{dml_benchmark}}.
##' @param combine.method method to combine cross-fitting repetitions. Either \code{"median"} (default, matching \code{summary.dml}/\code{summary.dml.bounds} elsewhere in this package) or \code{"mean"}.
##' @param na.rm logical. Should NA values be removed? Default is \code{TRUE}.
##' @param ... arguments passed to other methods.
##' @returns For \code{print}: the object, printed to console. For \code{summary}: an object of class \code{summary_dml_benchmark} holding a table of the benchmark components -- the leave-one-out gains \code{gain.Y} and \code{gain.D}, the alignment \code{rho}, and the bias contribution \code{delta} -- each with a standard error derived from its influence function and combined across cross-fitting repetitions. \code{gain.Y}/\code{gain.D} are not floored at 0 (see \code{\link{dml_benchmark}}), so a small negative combined value is possible and should be read as "statistically indistinguishable from zero," not as an error.
##' @rdname summary.dml_benchmark
##' @export
summary.dml_benchmark <- function(object, combine.method = c("median", "mean"),
                                  na.rm = TRUE, ...){
  combine.method <- match.arg(combine.method)
  combine <- if (combine.method == "mean") combine.mean else combine.median
  covars  <- names(object$benchmarks)

  rows <- lapply(covars, function(v) {
    est  <- object$benchmarks[[v]]        # per-rep point estimates (one row/rep)
    psis <- object$benchmarks_psis[[v]]   # per-rep influence functions
    reps <- seq_len(nrow(est))
    se_of <- function(psi.list) sapply(reps, function(i) psi.sd(psi.list[[i]]))
    # delta = theta.sj - theta.s (see bench_fun()), so its influence function
    # is psi.theta.s.wo - psi.theta.s; psi.sd() is sign-invariant, so the
    # order doesn't affect the resulting SE.
    se.delta <- sapply(reps, function(i)
      psi.sd(psis$psi.theta.s.wo[[i]] - psis$psi.theta.s[[i]]))
    cGY <- combine(est$gain.Y, se_of(psis$psi.GY))
    cGD <- combine(est$gain.D, se_of(psis$psi.GD))
    cRH <- combine(est$rho,    se_of(psis$psi.rho))
    cDL <- combine(est$delta,  se.delta)
    warn_if_negative_gain(v, unname(cGY["estimate"]), unname(cGD["estimate"]),
                          unname(cGY["se"]), unname(cGD["se"]))
    c(gain.Y = unname(cGY["estimate"]), se.gain.Y = unname(cGY["se"]),
      gain.D = unname(cGD["estimate"]), se.gain.D = unname(cGD["se"]),
      rho    = unname(cRH["estimate"]), se.rho    = unname(cRH["se"]),
      delta  = unname(cDL["estimate"]), se.delta  = unname(cDL["se"]))
  })

  tab <- do.call(rbind, rows)
  rownames(tab) <- covars

  # keep the influence functions on the summary object too (unlike a plain
  # aggregation, which would silently drop them) -- downstream code may want
  # to recompute SEs at a different combine.method without refitting anything
  out <- list(benchmarks = tab,
             benchmarks_psis = object$benchmarks_psis,
             combine.method = combine.method)
  class(out) <- "summary_dml_benchmark"
  out
}

##' @rdname summary.dml_benchmark
##' @export
print.summary_dml_benchmark <- function(x, digits = max(3L, getOption("digits") - 3L), ...){
  cat("\nDebiased Machine Learning: Covariate Benchmarks\n\n")
  print(round(x$benchmarks, digits), ...)
  cat("\nNote: benchmarks combined using the", x$combine.method, "method.\n")
  invisible(x)
}

# bench_plm <- function(model, benchmark_covariates) {
#   # if (is.null(model$results$main$all)) stop("Benchmarks implemented for ATE only. ATT/ATU coming soon.")
#   x <- model$data$x
#   which.not <- which(!benchmark_covariates %in% colnames(x))
#   if (any(which.not)){
#     stop("Covariates not found: ", paste(benchmark_covariates[which.not], collapse = ", "), ".")
#   }
#
#   resY.D   <- sapply(model$fits,
#                      function(x) lm(model$data$y - x$preds$yhat ~ model$data$d - x$preds$dhat)$res)
#   resD   <- sapply(model$fits, function(x) model$data$d - x$preds$dhat)
#
#   R2.Y <- (apply(resY.D, 2, function(x) max(1-mean(x^2)/var(model$data$y),0)))
#   R2.D <- (apply(resD, 2, function(x) max(1-mean(x^2)/var(model$data$d),0)))
#
#   theta.short <- extract_estimate(model$results$main[[1]], "theta.s")
#   benchmarks <- list()
#   for (i in seq_along(benchmark_covariates)){
#     covar <- benchmark_covariates[i]
#     cat("\n=== Computing benchmarks using covariate:", covar, " ===\n\n")
#     index.o <- which(colnames(x) == covar)
#     xo <- x[,-index.o]
#     model.call <- model$call
#     model.call["x"] <- call("xo")
#     model.wo <- eval(model.call)
#
#     resY.D.wo   <- sapply(model.wo$fits,
#                        function(x) lm(model.wo$data$y - x$preds$yhat ~ model.wo$data$d - x$preds$dhat)$res)
#     resD.wo   <- sapply(model.wo$fits, function(x) model.wo$data$d - x$preds$dhat)
#
#     R2.Ywo <- (apply(resY.D.wo, 2, function(x) max(1-mean(x^2)/var(model.wo$data$y),0)))
#     R2.Dwo <- (apply(resD.wo, 2, function(x) max(1-mean(x^2)/var(model.wo$data$d),0)))
#
#     ## Bias Decomposition
#     theta.short.wo <- extract_estimate(model.wo$results$main[[1]], "theta.s")
#     Bias <-  theta.short.wo - theta.short
#     V.g <- apply(resY.D.wo, 2, function(x) mean(x^2)) -
#       apply(resY.D,2, function(x) mean(x^2)) # var( g - g_s)
#     V.a <- apply(resD, 2, function(x) mean((x/mean(x^2))^2))-
#       apply(resD.wo, 2, function(x) mean((x/mean(x^2))^2)) # Var (a-a_s)
#     valid <- V.g > 0 & V.a > 0
#     Cor <- 0
#     Cor[valid] <- (abs(Bias[valid])/sqrt(V.g[valid]*V.a[valid]))
#     Cor <- pmin(1, Cor)
#     Cor <- Cor*sign(Bias)
#
#     #Gain metrics:
#     Gain.Y <- pmax(0, (R2.Y-R2.Ywo)/(1-R2.Y))
#     Gain.D <- pmax(0, (R2.D-R2.Dwo)/(1-R2.D))
#
#     bench <- data.frame(gain.Y = Gain.Y,
#                         gain.D =  Gain.D,
#                         rho = Cor,
#                         theta.s  = theta.short,
#                         theta.sj = theta.short.wo,
#                         delta = Bias)
#
#     benchmarks[[covar]] = bench
#   }
#   return(benchmarks)
# }

# gain in outcome-model fit attributable to the benchmark covariate: the
# relative increase in outcome residual variance (sigma2.s) from dropping
# it (V.g = sigma.sq.wo - sigma.sq). In population V.g >= 0 always (dropping
# a covariate can only weakly increase the true residual variance), but the
# *estimated* V.g can come out negative purely from DML sampling noise --
# sigma.sq and sigma.sq.wo come from two separate, independently cross-fitted
# dml() fits, so their difference doesn't benefit from any shared-fold noise
# cancellation. We deliberately do NOT floor this at 0: doing so would bias
# the estimator upward whenever the true gain is small (you'd keep the whole
# positive tail of the sampling distribution but clip the negative tail),
# and it would silently discard the information that a small negative value
# conveys (statistically indistinguishable from zero) in favor of an
# uninformative exact 0. This also keeps the point estimate consistent with
# its own influence function/SE (see psi.GY in bench_fun()), which is not
# floored either.
benchmark_gain_y <- function(sigma.sq, V.g) {
  V.g / sigma.sq
}

# gain in treatment/Riesz-representer precision attributable to the
# benchmark covariate: the relative increase in nu2.s from including it
# (V.a = nu.sq - nu.sq.wo). Same reasoning as benchmark_gain_y() above: V.a
# is >= 0 in population, but its DML estimate can be negative due to sampling
# noise, and is deliberately left unfloored for the same bias/consistency
# reasons.
benchmark_gain_d <- function(nu.sq.wo, V.a) {
  V.a / nu.sq.wo
}

# warns whenever gain.Y and/or gain.D (both non-negative in the population by
# construction -- see benchmark_gain_y()/benchmark_gain_d()) come out
# negative for a covariate. Checks both quantities together and issues a
# single warning covering whichever of the two are negative, so the
# accompanying guidance isn't printed twice when both happen to be negative
# at once. This is deliberately just a warning, not an error or a silent
# floor: it does not tell the researcher what to conclude or do about it,
# only reports the plain facts (which quantity, the value(s), and the SE if
# available), since a negative estimate here reflects sampling noise around
# a true value near zero, not an impossible population value or necessarily
# a mistake.
warn_if_negative_gain <- function(covar, gain.Y, gain.D, se.gain.Y = NULL, se.gain.D = NULL,
                                  suggest_solutions = TRUE) {
  describe <- function(quantity, value, se) {
    neg <- which(value < 0)
    if (length(neg) == 0L) return(NULL)
    detail <- if (is.null(se)) {
      rep_word <- if (length(value) == 1L) "repetition" else "repetitions"
      val_word <- if (length(neg) == 1L) "value" else "values"
      paste0(length(neg), " of ", length(value), " cross-fitting ", rep_word, ", ", val_word, ": ",
            paste(sprintf("%.4f", value[neg]), collapse = ", "))
    } else {
      paste0("estimate = ", sprintf("%.4f", value[neg]), ", SE = ", sprintf("%.4f", se[neg]))
    }
    paste0(quantity, " (", detail, ")")
  }
  parts <- c(describe("gain.Y", gain.Y, se.gain.Y), describe("gain.D", gain.D, se.gain.D))
  if (length(parts) == 0L) return(invisible(NULL))
  verb <- if (length(parts) == 1L) "is" else "are"
  # Deliberately calibrated: the true value cannot be negative (certain,
  # proven), but "sampling noise" is the most plausible explanation for why
  # the estimate is, not a verified diagnosis -- this warning only observes
  # the symptom, it cannot rule out other sources of estimation error.
  msg <- paste0("For covariate '", covar, "', ", paste(parts, collapse = " and "), " ", verb,
               " negative. gain.Y and gain.D cannot be negative in the population (see ",
               "?dml_benchmark); the true value here is most plausibly near zero, with ",
               "sampling noise producing a negative estimate -- though this warning cannot ",
               "rule out other sources of estimation error.")
  # suggest_solutions is FALSE at the per-repetition call site (bench_fun())
  # so that a normal fit -> dml_benchmark() -> summary()/print() workflow --
  # which triggers this warning at both the per-repetition and combined
  # level -- only prints this guidance once, attached to the combined
  # estimate (the number the user will actually read or act on), instead of
  # repeating it at both stages.
  if (suggest_solutions) {
    msg <- paste0(msg, " Two ways to address this: (1) reconsider whether this covariate is a ",
                 "substantively strong benchmark -- try a different one, or group several ",
                 "related covariates together (see the 'benchmark_covariates' argument of ",
                 "dml_benchmark()); or (2) reduce estimation noise by re-fitting the underlying ",
                 "dml() model with more cross-fitting repetitions (cf.reps) or a different ",
                 "learner (yreg/dreg), then re-running dml_benchmark().")
  }
  warning(msg, call. = FALSE)
}

# alignment (rho) between the covariate's effect on the outcome and on
# selection: |Bias|/sqrt(V.g*V.a), signed by sign(Bias) and capped at 1 in
# magnitude. Forced to 0 wherever `valid` is FALSE (i.e., V.g or V.a is
# non-positive -- dropping the covariate did not actually reduce outcome
# fit/treatment precision, so there is no meaningful bias to attribute to
# it).
benchmark_rho <- function(Bias, V.g, V.a, valid) {
  Cor <- rep(0, length(Bias))
  Cor[valid] <- abs(Bias[valid]) / sqrt(V.g[valid] * V.a[valid])
  Cor <- pmin(1, Cor)
  Cor * sign(Bias)
}

# maps a dml target to the name of its slot in model$results$main (only
# meaningful for npm models -- see resolve_benchmark_slot())
.target_to_slot <- c(ate = "all", att = "treat", atu = "untr")

# resolves which slot of model$results$main to read from. plm models always
# store their (single) results under "all", regardless of the declared
# target -- dml() accepts target = "att"/"atu" for a plm model without
# erroring, but only npm's ate.att.atu.npm() actually computes distinct
# target-specific slots, so model type is checked explicitly here rather
# than trusting target alone. Errors if the model was fit with more than one
# target (benchmarking a multi-target fit is ambiguous: which slot's
# theta.s/sigma2.s/nu2.s should the leave-one-out comparison use?).
#
# Note: the resulting rho/psi.rho have only been verified (against the
# manuscript's published numbers) for att-slot models; ate/atu are computed
# with the same, unmodified formula but that has not been independently
# verified for those targets specifically.
resolve_benchmark_slot <- function(model) {
  if (identical(model$info$model, "plm")) {
    return("all")
  }
  slot <- unname(.target_to_slot[model$info$target])
  if (length(slot) != 1L || is.na(slot)) {
    stop("dml_benchmark() requires a model fit with a single target ",
        "('ate', 'att', or 'atu').")
  }
  slot
}

# normalizes benchmark_covariates into a named list of column-name groups. A
# character vector benchmarks each column on its own; a named list lets a
# group of columns (e.g. the dummies of a factor) be dropped together in the
# leave-one-out refit. Element names become the row labels; unnamed elements
# are labelled by the column name (singletons) or the columns joined by "+".
normalize_benchmark_groups <- function(benchmark_covariates, x) {
  groups <- if (is.list(benchmark_covariates)) benchmark_covariates
           else as.list(benchmark_covariates)

  group.names <- names(groups)
  if (is.null(group.names)) group.names <- rep("", length(groups))

  for (i in seq_along(groups)) {
    cols <- groups[[i]]
    if (!is.character(cols) || length(cols) < 1L) {
      stop("Each entry of 'benchmark_covariates' must be a non-empty ",
          "character vector of column names.")
    }
    if (is.na(group.names[i]) || group.names[i] == "") {
      group.names[i] <- if (length(cols) == 1L) cols else paste(cols, collapse = "+")
    }
  }
  names(groups) <- group.names

  all.cols  <- unlist(groups, use.names = FALSE)
  which.not <- which(!all.cols %in% colnames(x))
  if (length(which.not) > 0) {
    stop("Covariates not found: ", paste(all.cols[which.not], collapse = ", "), ".")
  }

  groups
}

bench_fun <- function(model, benchmark_covariates){

  slot <- resolve_benchmark_slot(model)
  x <- model$data$x
  covariate.groups <- normalize_benchmark_groups(benchmark_covariates, x)

  nu.sq <- extract_estimate(model$results$main[[slot]], param = "nu2.s")
  sigma.sq <- extract_estimate(model$results$main[[slot]], param = "sigma2.s")

  # resY  <- sapply(model$fits, function(x)model$data$y-x$preds$yhat)
  # R2.Y  <- apply(resY, 2, function(x) max(1-var(x)/var(model$data$y),0))

  theta.short <- extract_estimate(model$results$main[[slot]], "theta.s")

  # short IFs
  psi.theta.s  <- lapply(model$results$main[[slot]], function(x) x$psis$psi.theta.s)
  psi.sigma2.s <- lapply(model$results$main[[slot]], function(x) x$psis$psi.sigma2.s)
  psi.nu2.s    <- lapply(model$results$main[[slot]], function(x) x$psis$psi.nu2.s)

  benchmarks <- list()
  benchmarks_psis <- list()
  for (covar in names(covariate.groups)) {
    cols <- covariate.groups[[covar]]
    cat("\n=== Computing benchmarks using covariate:", covar, " ===\n\n")
    index.o <- which(colnames(x) %in% cols)   # drop all columns in the group
    xo <- x[, -index.o, drop = FALSE]
    model.call <- model$call
    model.call["x"] <- call("xo")
    model.wo <- eval(model.call)

    nu.sq.wo <- extract_estimate(model.wo$results$main[[slot]], param = "nu2.s")
    sigma.sq.wo <- extract_estimate(model.wo$results$main[[slot]], param = "sigma2.s")

    # resY.wo  <- sapply(model.wo$fits, function(x) model.wo$data$y - x$preds$yhat)
    # R2.Y.wo  <- apply(resY.wo, 2, function(x) max(1-var(x)/var(model.wo$data$y),0))

    ## (Debiased) Bias Decomposition
    theta.short.wo <- extract_estimate(model.wo$results$main[[slot]], "theta.s")

    # benchmark IFs
    psi.theta.s.wo  <- lapply(model.wo$results$main[[slot]], function(x) x$psis$psi.theta.s)
    psi.sigma2.s.wo <- lapply(model.wo$results$main[[slot]], function(x) x$psis$psi.sigma2.s)
    psi.nu2.s.wo    <- lapply(model.wo$results$main[[slot]], function(x) x$psis$psi.nu2.s)

    # Bias = theta.short.wo - theta.short (without minus with). This
    # specific direction is load-bearing: benchmark_rho()/psi.rho below
    # implement rho = Bias/sqrt(V.g*V.a) with NO leading minus sign, which
    # only reproduces the correct (true-correlation) sign of rho when Bias
    # is defined this way. It also matches the manuscript's own definition
    # of the benchmarking "change in estimate" (Appendix E.1, Chernozhukov
    # et al., 2026): Delta_{s,j} := Em(W,g_{s,-j}) - Em(W,g_s), i.e.
    # theta_{s,j} (without the covariate) minus theta_s (with it) -- the
    # same without-minus-with direction as Bias. So `delta` is just `Bias`,
    # not its negation; verified against Table 3's reported values (e.g.
    # income: Delta_hat = 3,349, positive, matching theta.sj > theta.s).
    Bias  <- theta.short.wo - theta.short
    Delta <- Bias  # theta.short.wo - theta.short, matches manuscript's Delta_{s,j}

    # V.g <- apply(resY.wo,2,var) - apply(resY,2,var)
    V.g <- sigma.sq.wo - sigma.sq
    V.a <- nu.sq - nu.sq.wo
    valid <- V.g > 0 & V.a > 0

    Cor <- benchmark_rho(Bias, V.g, V.a, valid)

    #(1- R^2_{a~a_s}) =  (Ea^2 - Ea_s^2)/ E a^2
    Gain.Y <- benchmark_gain_y(sigma.sq, V.g)
    Gain.D <- benchmark_gain_d(nu.sq.wo, V.a)
    warn_if_negative_gain(covar, Gain.Y, Gain.D, suggest_solutions = FALSE)

    bench <- data.frame(gain.Y = Gain.Y,
                        gain.D = Gain.D,
                        rho = Cor,
                        theta.s  = theta.short,
                        theta.sj = theta.short.wo,
                        delta = Delta)

    benchmarks[[covar]] = bench

    ## bounds IFs
    psi.GY <- Map(function(psi.wo, psi, s, s.wo) {
      (s * psi.wo - s.wo * psi) / (s^2)
    }, psi.sigma2.s.wo, psi.sigma2.s, sigma.sq, sigma.sq.wo)

    psi.GD <- Map(function(psi.wo, psi, nu, nu.wo) {
      (nu.wo * psi - nu * psi.wo) / (nu.wo^2)
    }, psi.nu2.s.wo, psi.nu2.s, nu.sq, nu.sq.wo)

    psi.rho <- lapply(seq_len(length(valid)), function(k) {
      if (!valid[k]) {
        rep(0, nrow(x))
      } else {
        ((psi.theta.s.wo[[k]] - psi.theta.s[[k]]) /
           (sqrt(V.g[k] * V.a[k]))) -
          ((Bias[k] * (psi.sigma2.s.wo[[k]] - psi.sigma2.s[[k]])) /
             (2 * (V.g[k]^(3/2)) * sqrt(V.a[k]))) -
          ((Bias[k] * (psi.nu2.s[[k]] - psi.nu2.s.wo[[k]])) /
             (2 * (V.a[k]^(3/2)) * sqrt(V.g[k])))
      }
    })

    benchmarks_psis[[covar]] = list(
      psi.theta.s  = psi.theta.s,
      psi.sigma2.s = psi.sigma2.s,
      psi.nu2.s    = psi.nu2.s,

      psi.theta.s.wo  = psi.theta.s.wo,
      psi.sigma2.s.wo = psi.sigma2.s.wo,
      psi.nu2.s.wo    = psi.nu2.s.wo,

      psi.GY  = psi.GY,
      psi.GD  = psi.GD,
      psi.rho = psi.rho
    )
  }
  return(list(
    benchmarks      = benchmarks,
    benchmarks_psis = benchmarks_psis))
}
