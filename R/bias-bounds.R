
# bounds functions --------------------------------------------------------
bias.factor <- function(cf.y = 0.03, cf.d = 0.04, rho2 = 1){
  sqrt(rho2*cf.y*(cf.d/(1 - cf.d)))
}


# Extracts and precomputes the quantities that bounds() needs but that do
# NOT depend on the hypothesized cf.y/cf.d/rho2 scenario -- theta.s, nu2.s
# (for the validity check), the influence functions, S = sqrt(sigma2.s*nu2.s),
# psi.S2, and se.theta.s. None of these formulas involve cf.y/cf.d/rho2, so
# their value is identical no matter what confounding scenario is later
# evaluated. Computing this once and reusing it across repeated calls to
# bounds() with different cf.y/cf.d -- e.g. inside an optim()-based search,
# as in robustness_value()/extreme_robustness_value() -- avoids redundant
# recomputation of psi.S2 and se.theta.s (both O(n) in the sample size) on
# every one of the ~20-50 evaluations such a search typically needs.
prep_bounds <- function(short.results) {
  theta.s  <- short.results$estimates$theta.s
  sigma2.s <- short.results$estimates$sigma2.s
  nu2.s    <- short.results$estimates$nu2.s

  # Checked here (once per prep_bounds() call) rather than in bounds()
  # itself: nu2.s does not depend on cf.y/cf.d/rho2, so when a `fixed`
  # built from this call is reused across many bounds() calls -- e.g. the
  # ~20-50 optim() evaluations per row in robustness_value()/
  # extreme_robustness_value() -- the underlying fact doesn't change from
  # one call to the next, and shouldn't be re-warned every time.
  if (nu2.s < 0) {
    warning("Unable to compute bounds because nu^2 is negative,
            this is an indication of a very poor fitting of the treatment model.
            Consider choosing a different learner for the treatment.")
  }

  psi.theta.s  <- short.results$psis$psi.theta.s
  psi.sigma2.s <- short.results$psis$psi.sigma2.s
  psi.nu2.s    <- short.results$psis$psi.nu2.s

  S      <- sqrt(sigma2.s*nu2.s)
  psi.S2 <- (sigma2.s*(psi.nu2.s) + nu2.s*psi.sigma2.s)

  list(theta.s = theta.s, nu2.s = nu2.s,
      psi.theta.s = psi.theta.s, psi.sigma2.s = psi.sigma2.s, psi.nu2.s = psi.nu2.s,
      S = S, psi.S2 = psi.S2, se.theta.s = psi.sd(psi.theta.s))
}

# `fixed`, if supplied, must be the output of prep_bounds(short.results) --
# passing it in skips recomputing the cf.y/cf.d/rho2-invariant quantities
# (see prep_bounds() above). When fixed = NULL (the default, and the
# behavior of every pre-existing call site), this computes them inline
# exactly as before, so nothing changes for any caller that doesn't opt in.
bounds <- function(short.results, cf.y = 0.04, cf.d = 0.03, rho2 = 1, fixed = NULL){

  if (is.null(fixed)) fixed <- prep_bounds(short.results)

  theta.s <- fixed$theta.s
  nu2.s   <- fixed$nu2.s

  # See prep_bounds(): the nu2.s < 0 check (and its warning) happens there,
  # once per `fixed`, rather than being repeated on every bounds() call
  # that reuses it.

  # bias bound
  bf <- bias.factor(cf.y = cf.y, cf.d = cf.d, rho2 = rho2)
  S  <- fixed$S
  bias.bound <- S*bf

  # plug-in bounds
  theta.m <- theta.s - bias.bound
  theta.p <- theta.s + bias.bound

  # short IFs
  psi.theta.s  <- fixed$psi.theta.s
  psi.sigma2.s <- fixed$psi.sigma2.s
  psi.nu2.s    <- fixed$psi.nu2.s

  # bounds IFs
  psi.S2         <- fixed$psi.S2
  psi.bias.bound <- (bf/2) * (1/S) * psi.S2
  psi.theta.m    <- psi.theta.s - psi.bias.bound
  psi.theta.p    <- psi.theta.s + psi.bias.bound

  out <- list(
    psis      = list(psi.theta.s    = psi.theta.s,
                     psi.sigma2.s   = psi.sigma2.s,
                     psi.nu2.s      = psi.nu2.s,
                     psi.bias.bound = psi.bias.bound,
                     psi.theta.m    = psi.theta.m,
                     psi.theta.p    = psi.theta.p),

    estimates = list(theta.s = theta.s,
                     se.theta.s = fixed$se.theta.s,
                     bias.bound = bias.bound,
                     se.bias.bound = psi.sd(psi.bias.bound),
                     theta.m = theta.m,
                     se.theta.m = psi.sd(psi.theta.m),
                     theta.p = theta.p,
                     se.theta.p = psi.sd(psi.theta.p))
  )
  return(out)
}





extract_coefs <- function(bounds.results){
  coefs.names <- c("theta.s", "bias.bound", "theta.m", "theta.p")
  names(coefs.names) <- coefs.names
  lapply(coefs.names,
         function(param) combine.cross.fits(bounds.results, param = param))
}

get_bounds <- function(bounds, combine.method = "mean"){
  f <- function(coefs) t(sapply(coefs, function(x) x[combine.method, ]))
  lapply(bounds$coefs, f)
}

##' Bounds on the omitted variable bias for causal machine learning
##'
##' @param model an object of class \code{\link{dml}} or \code{\link{dml_bounds}}.
##' @param cf.y (nonparametric) partial R2 of the omitted variables with the outcome. Must be a number between (0, 1).
##' @param cf.d how much variation latent variables create in the Riesz Representer of the target parameters. Must be a number between (0, 1). When the target of interest is the ATE in a partially linear model, this corresponds to the partial R2 of omitted variables with the treatment. When the target of interest is the ATE in a non-parametric model with a binary treatment, this corresponds to the gains in precision (i.e, 1/variance) when predicting who is assigned to treatment.
##' @param rho2 degree of adversity. Default is \code{rho=1}, which assumes the maximum degree of adversity of confounding.
##' @param ... arguments passed to other methods.
##' @param combine.method method to combine the results of each repetition. Options are \code{mean} and \code{median}. Default is \code{median}.
##' @param level confidence level. Default is \code{0.95}.
##' @param only restricts computation to a single main slot and/or a single group,
##'   instead of every slot/group in \code{model$results}. Must be a list with
##'   elements \code{main} (a slot name, e.g. \code{"all"}, or \code{NULL} to skip
##'   all main slots) and \code{groups} (a group name, or \code{NULL} to skip all
##'   groups). Default \code{NULL} computes everything, exactly as before. Intended
##'   for internal use by \code{\link{robustness_value}}/\code{\link{extreme_robustness_value}},
##'   whose \code{optim()} search only ever needs one row's result per evaluation.
##' @returns For \code{dml_bounds}: an object of class \code{dml.bounds}. For \code{confidence_bounds}: a matrix or numeric vector of confidence bounds.
##' @export
dml_bounds <- function(model, cf.y, cf.d, rho2 = 1, fixed = NULL, only = NULL){

  # bounds for
  out <- list()
  out$info <- list(cf.y = cf.y,
                     cf.d = cf.d,
                     rho2 = rho2)

  out$dml.fit <- model

  # main <- model$results$main
  # bounds.results   <- lapply(main, bounds, cf.y = cf.y, cf.d = cf.d, rho2 = rho2)
  # out$results$main <- bounds.results
  #
  # main.coefs <- extract_coefs(bounds.results)
  # out$coefs$main <- main.coefs
  main <- model$results$main
  fixed.main <- fixed$main
  if (!is.null(only)) {
    if (is.null(only$main)) {
      main <- NULL
    } else {
      # Name-based subsetting (not positional): keeps `main` and `fixed.main`
      # aligned to the SAME slot even though `only$main` need not be the
      # first slot in model$results$main -- positional Map() below would
      # otherwise silently pair the wrong slot's cached `fixed` quantities
      # with this slot's bounds() computation.
      main <- main[only$main]
      if (!is.null(fixed.main)) fixed.main <- fixed.main[only$main]
    }
  }
  if (!is.null(main)) {
    main.bounds <- if (is.null(fixed.main)) {
      lapply(main,
            function(x) lapply(x,
                               bounds, cf.y = cf.y, cf.d = cf.d, rho2 = rho2))
    } else {
      Map(function(x, fx) Map(function(xx, fxx) bounds(xx, cf.y = cf.y, cf.d = cf.d, rho2 = rho2, fixed = fxx), x, fx),
         main, fixed.main[names(main)])
    }
    out$results$main <- main.bounds
    out$coefs$main <- lapply(main.bounds, extract_coefs)
  }

  groups <- model$results$groups
  fixed.groups <- fixed$groups
  if (!is.null(only)) {
    if (is.null(only$groups)) {
      groups <- NULL
    } else {
      groups <- groups[only$groups]
      if (!is.null(fixed.groups)) fixed.groups <- fixed.groups[only$groups]
    }
  }
  if (!is.null(groups)) {
    groups.bounds <- if (is.null(fixed.groups)) {
      lapply(groups,
            function(x) lapply(x,
                               bounds, cf.y = cf.y, cf.d = cf.d, rho2 = rho2))
    } else {
      Map(function(x, fx) Map(function(xx, fxx) bounds(xx, cf.y = cf.y, cf.d = cf.d, rho2 = rho2, fixed = fxx), x, fx),
         groups, fixed.groups[names(groups)])
    }
    out$results$groups <- groups.bounds
    out$coefs$groups <- lapply(groups.bounds, extract_coefs)
    }

  class(out) <- "dml.bounds"
  return(out)
}

# Precomputes, once, the prep_bounds() quantities for every target slot and
# every cross-fitting repetition in a fitted dml() model (and any groupwise
# results), mirroring the structure of model$results$main / $groups. Passing
# the result into dml_bounds() (or confidence_bounds.dml(), rv_fun(),
# xrv_fun()) lets repeated calls -- e.g. from an optim()-based search over
# cf.y/cf.d, as in robustness_value()/extreme_robustness_value() -- skip
# recomputing these fixed quantities on every evaluation.
prep_dml_bounds <- function(model) {
  out <- list()
  main <- model$results$main
  if (!is.null(main)) {
    out$main <- lapply(main, function(x) lapply(x, prep_bounds))
  }
  groups <- model$results$groups
  if (!is.null(groups)) {
    out$groups <- lapply(groups, function(x) lapply(x, prep_bounds))
  }
  out
}

# For each row of confint(model)/coef(model) -- main "ate.*" rows first, then
# "gate.*" group rows, in that order (see se.dml()'s c(ate = , gate = )
# construction and ate.att.atu.npm()/group.ate.*()'s preserved ordering) --
# returns the `only` value that restricts dml_bounds() to just that row's
# slot/group. Used by robustness_value()/extreme_robustness_value() so their
# optim() search, which only ever needs one row's result per evaluation,
# doesn't force dml_bounds() to compute every other slot and group too.
row_only_list <- function(model) {
  main.slots <- if (identical(model$info$model, "plm")) {
    rep("all", length(model$info$target))
  } else {
    unname(.target_to_slot[model$info$target])
  }
  group.names <- names(model$results$groups)
  c(
    lapply(main.slots, function(s) list(main = s, groups = NULL)),
    lapply(group.names, function(g) list(main = NULL, groups = g))
  )
}

# TRUE if the row identified by `only` (see row_only_list()) has a negative
# nu2.s in any cross-fitting repetition. S = sqrt(sigma2.s * nu2.s) is then
# NaN for that repetition regardless of cf.y/cf.d (bias.bound = S * bf is
# NaN for ANY bf), and combine.mean()/combine.median() have no na.rm, so a
# single bad repetition makes the row's combined bound NA for every
# cf.y/cf.d -- the confidence_bounds() objective robustness_value()/
# extreme_robustness_value() feed into optim() is then a CONSTANT NA over
# the whole search domain. optim(method = "Brent") doesn't error or return
# NA on a constant-NA objective: it deterministically converges to its
# `upper` bound (verified empirically), which looks like an ordinary RV/XRV
# value (e.g. 1) but carries no information about the actual data -- an
# honest NA is preferable to a plausible-looking number that isn't real.
row_is_undefined <- function(fixed, only) {
  reps <- if (!is.null(only$main)) fixed$main[[only$main]] else fixed$groups[[only$groups]]
  any(vapply(reps, function(r) r$nu2.s < 0, logical(1)))
}


#' Compute confidence bounds
#' @description Computes confidence bounds on the target parameter of interest accounting for omitted variable biases.
#'
#' @rdname dml_bounds
#' @export
confidence_bounds <- function(...){
  UseMethod("confidence_bounds")
}



##' @param theta.s short estimate of the target parameter (used as the first argument for the numeric method).
##' @param S2 estimated variance product (sigma2 * nu2).
##' @param se.theta.s standard error of the short estimate.
##' @param se.S2 standard error of S2.
##' @param cov.theta.S2 covariance of theta.s and S2.
##' @param return character vector specifying which bounds to return. Options are \code{"lwr"} and \code{"upr"}.
#' @export
#' @rdname dml_bounds
confidence_bounds.numeric <- function(theta.s, S2,
                                      se.theta.s, se.S2,
                                      cov.theta.S2,
                                      cf.y, cf.d,
                                      rho2 = 1,
                                      combine.method = "median",
                                      level = 0.95, ...){
  k = bias.factor(cf.y = cf.y, cf.d = cf.d, rho2 = rho2)
  se.m <- sqrt((se.theta.s)^2 + (k^2/(4*S2))*se.S2^2 - (k/sqrt(S2))*cov.theta.S2)
  se.p <- sqrt((se.theta.s)^2 + (k^2/(4*S2))*se.S2^2 + (k/sqrt(S2))*cov.theta.S2)
  theta.m <- theta.s - k*sqrt(S2)
  theta.p <- theta.s + k*sqrt(S2)
  level[level < 0.5] <- 0.5
  t_crit <- qnorm(level)
  lwr <- combine.median(theta.m, se.m)
  upr <- combine.median(theta.p, se.p)
  lwr <- unname(lwr["estimate"] - t_crit*lwr["se"])
  upr <- unname(upr["estimate"] + t_crit*upr["se"])
  c(lwr = lwr, upr = upr)
}

#' @export
#' @rdname dml_bounds
confidence_bounds.dml <- function(model,
                                  cf.y,
                                  cf.d,
                                  rho2 = 1,
                                  level = 0.95,
                                  combine.method = "median",
                                  fixed = NULL,
                                  only = NULL, ...){
  object <- dml_bounds(model, cf.y = cf.y, cf.d = cf.d, rho2 = rho2, fixed = fixed, only = only)
  confidence_bounds(object, level = level,combine.method = combine.method, ...)
}



#' @export
#' @rdname dml_bounds
confidence_bounds.dml.bounds <- function(model,
                                         cf.y = NULL,
                                         cf.d = NULL,
                                         rho2 = NULL,
                                         level = 0.95,
                                         combine.method = "median",
                                         return = c("lwr", "upr"),
                                         ...){

  if (!is.null(cf.y) | !is.null(cf.d) | !is.null(rho2)) {

    if (is.null(cf.y)) {
      cf.y <- model$info$cf.y
    }

    if (is.null(cf.d)) {
      cf.d <- model$info$cf.d
    }

    if (is.null(rho2)) {
      rho2 <- model$info$rho2
    }

    new_bounds <- dml_bounds(model$dml.fit, cf.y = cf.y, cf.d = cf.d, rho2 = rho2)
    return(confidence_bounds(new_bounds, combine.method = combine.method, return = return))

  }

  level2 = max(0, 1 - (1 - level)*2)

  confs <- confint(model, level = level2, combine.method = combine.method)

  if (is.list(confs)) {
    out <- t(sapply(confs, function(x) c(lwr = x["theta.m",1], upr = x["theta.p",2])))
  } else {
    out <- rbind(ate = c(lwr = confs["theta.m",1], upr = confs["theta.p",2]))
  }
  out <- out[, return, drop = FALSE]
  attr(out, "conf.levels") <- c(point = level, region = level2)
  attr(out, "sens.param")  <- model$info
  class(out) <- c("confidence.bounds", "matrix")
  out
}


rv_fun <- function(dml.fit, rv, par, side = "lwr", theta = 0, alpha = 0.05, fixed = NULL, only = NULL){
  (confidence_bounds(dml.fit,  cf.y = rv,cf.d = rv, level = 1 - alpha, fixed = fixed, only = only)[par,side] - theta)^2
}

# objective for the extreme robustness value: fixes cf.y = 1 (the maximum
# possible outcome-side confounding strength) and searches only over cf.d,
# collapsing the two-dimensional RV search into one dimension. This is the
# single-parameter simplification in Chernozhukov et al. (2026), Section 2.4:
# setting both rho2 and cf.y to their trivial bound of 1 reduces the bound to
# |theta - theta.s| <= sqrt(cf.d/(1-cf.d)) * S.
xrv_fun <- function(dml.fit, xrv, par, side = "lwr", theta = 0, alpha = 0.05, rho2 = 1, fixed = NULL, only = NULL){
  (confidence_bounds(dml.fit, rho2 = rho2, cf.y = 1, cf.d = xrv, level = 1 - alpha, fixed = fixed, only = only)[par,side] - theta)^2
}

##' Computes Extreme Robustness Values for Debiased Machine Learning
##'
##' @description
##' This function computes the extreme robustness value of a target parameter estimated via debiased machine learning.
##'
##' The extreme robustness value describes the minimum strength of association (parameterized in terms of partial R2) that omitted variables would need to have with the Riesz Representer alone -- i.e. fixing \code{cf.y = 1}, the maximum possible outcome-side confounding strength -- so that the confidence bounds for the target parameter includes zero (or another threshold of interest). This is the single-parameter simplification discussed in Chernozhukov et al. (2026), Section 2.4, obtained by setting both \code{rho2} and \code{cf.y} to their trivial bound of 1.
##'
##' @returns A named numeric vector of extreme robustness values.
##' @export
extreme_robustness_value <- sensemakr::extreme_robustness_value

##' @rdname extreme_robustness_value
##' @param model an object of class \code{\link{dml}}.
##' @param theta the null hypothesis of interest for the target parameter theta. Default is \code{theta = 0} (zero null hypothesis).
##' @param alpha significance level. Default is \code{alpha = 0.05}. Setting \code{alpha = 1} uses a closed-form point XRV (no confidence-level adjustment) instead of the numerical search used for \code{alpha < 1}.
##' @param rho2 degree of adversity. Default is \code{rho2 = 1}, which assumes the maximum degree of adversity of confounding.
##' @inheritParams summary.dml
##' @exportS3Method sensemakr::extreme_robustness_value dml
extreme_robustness_value.dml <- function(model, theta = 0, alpha = 0.05, rho2 = 1, ...){
  conf <- confint(model, level = 1 - alpha, ...)
  out <- setNames(rep(NA, nrow(conf)), rownames(conf))
  # Precomputed once, outside the optim() search below, so the cf.y/cf.d-
  # invariant quantities (see prep_bounds()) aren't recomputed on every one
  # of the ~20-50 evaluations Brent's method needs per row. Unused (and
  # harmless to compute) for rows resolved by the alpha == 1 closed form.
  fixed <- prep_dml_bounds(model)

  # `only` for each row of conf, in the same order: main rows (one per
  # model$info$target slot) first, then group rows -- see row_only_list().
  # NOTE: rownames(conf)[i] is NOT simply the target name (e.g. "ate") --
  # se.dml() builds it as c(ate = <named vector by slot>), and R's naming
  # rule for c(prefix = named_vector) concatenates names with a dot, so a
  # PLM ATE row is actually named "ate.all", not "ate". row_only_list()
  # derives the slot from model$info$target/model$info$model instead of
  # parsing that compound string. Passing `only` into the general search
  # below also restricts dml_bounds() to just this row's slot/group,
  # instead of recomputing every slot and every group on every optim()
  # evaluation.
  row.only <- row_only_list(model)

  for (i in 1:nrow(conf)) {
    if (conf[i,1] <= theta & theta <= conf[i,2]) {
      out[i] <- 0
      next
    }
    # See row_is_undefined(): a negative nu2.s in any cf.rep makes this
    # row's bound NA for every cf.y/cf.d, so neither the closed form nor
    # optim()'s general search can return a meaningful value here.
    if (row_is_undefined(fixed, row.only[[i]])) {
      out[i] <- NA_real_
      next
    }
    side <- ifelse(theta < conf[i,1], "lwr", "upr")
    slot.i <- row.only[[i]]$main  # NULL for group rows -- closed form only covers main rows
    if (alpha == 1 && !is.null(slot.i)) {
      # closed-form point XRV for row i's target: f0 = |theta - theta.s|/sqrt(S2),
      # XRV = f0^2/(1+f0^2) (see Section 2.4's single-parameter bound, solved
      # for the critical cf.d at the theta boundary). Combine S2 across
      # cf.reps via the median, matching the package's default combine.method.
      S2.i <- stats::median(extract_estimate(model$results$main[[slot.i]], "S2"))
      f0   <- unname(abs(theta - coef(model)[rownames(conf)[i]]) / sqrt(S2.i))
      out[i] <- f0^2 / (1 + f0^2)
    } else {
      fn <- function(xrv) xrv_fun(xrv, dml.fit = model, par = names(out)[i],
                                  side = side, theta = theta, alpha = alpha,
                                  rho2 = rho2, fixed = fixed, only = row.only[[i]])
      out[i] <- optim(par = c(0.01), fn, lower = 0, upper = 1, method = "Brent")$par
    }
  }
  return(out)
}


##' Computes Robustness Values for Debiased Machine Learning
##'
##' @description
##' This function computes the robustness value of a target parameter estimated via debiased machine learning.
##'
##' The robustness value describes the minimum strength of association (parameterized in terms of  R2) that omitted variables would need to have both with the outcome and with the Riesz Representer so that the confidence bounds for the target parameter includes zero (or another threshold of interest).
##'
##'
##' @returns A named numeric vector of robustness values.
##' @export
robustness_value <- sensemakr::robustness_value

##' @rdname robustness_value
##' @param model an object of class \code{\link{dml}} or \code{\link{dml_bounds}}.
##' @param theta the null hypothesis of interest for the target parameter theta. Default is \code{theta =0} (zero null hypothesis).
##' @param alpha significance level. Default is \code{alpha = 0.05}.
##' @inheritParams summary.dml
##' @exportS3Method sensemakr::robustness_value dml
robustness_value.dml <- function(model, theta = 0, alpha = 0.05, ...){
  conf <- confint(model, level = 1 - alpha,...)
  out <- setNames(rep(NA,nrow(conf)), rownames(conf))
  # See prep_dml_bounds()/prep_bounds(): precomputed once so the cf.y/cf.d-
  # invariant quantities aren't redone on every optim() evaluation below.
  fixed <- prep_dml_bounds(model)
  # See row_only_list(): restricts dml_bounds() to just row i's slot/group,
  # instead of recomputing every slot and every group on every evaluation.
  row.only <- row_only_list(model)
  for (i in 1:nrow(conf)) {
    if (conf[i,1] <= theta & theta <= conf[i,2]) {
      out[i] <- 0
      next
    }
    # See row_is_undefined(): a negative nu2.s in any cf.rep makes this
    # row's bound NA for every cf.y/cf.d, so optim()'s search would be
    # over a constant-NA objective -- Brent doesn't error on that, it
    # deterministically returns its `upper` bound as a value that looks
    # like a real RV but isn't one.
    if (row_is_undefined(fixed, row.only[[i]])) {
      out[i] <- NA_real_
      next
    }
    side <- ifelse(theta < conf[i,1], "lwr", "upr")
    fn <- function(rv) rv_fun(rv, dml.fit = model, par = names(out)[i],
                              side = side, theta = theta, alpha = alpha, fixed = fixed,
                              only = row.only[[i]])
    out[i] <- optim(par = c(0.01), fn, lower = 0, upper = 1, method = "Brent")$par
  }
  return(out)
  # grid <- seq(0, 0.99,by = 0.001)
  # values <- mapply(function(x,y) confidence_bounds(dml.fit, cf.y= x, cf.d = y), x = grid, y = grid)
  # rv.idx <- which(values[1,] <= theta & theta <= values[2,])[1]
  # grid[rv.idx]
}

##' @rdname robustness_value
##' @exportS3Method sensemakr::robustness_value dml.bounds
robustness_value.dml.bounds <- function(model, theta = 0, alpha = 0.05, ...){
  model <- model$dml.fit
  conf <- confint(model, level = 1 - alpha,...)
  out <- setNames(rep(NA, nrow(conf)), rownames(conf))
  # See prep_dml_bounds()/prep_bounds(): precomputed once so the cf.y/cf.d-
  # invariant quantities aren't redone on every optim() evaluation below.
  fixed <- prep_dml_bounds(model)
  # See row_only_list(): restricts dml_bounds() to just row i's slot/group,
  # instead of recomputing every slot and every group on every evaluation.
  row.only <- row_only_list(model)
  for (i in 1:nrow(conf)) {
    if (conf[i,1] <= theta & theta <= conf[i,2]) {
      out[i] <- 0
      next
    }
    # See row_is_undefined(): a negative nu2.s in any cf.rep makes this
    # row's bound NA for every cf.y/cf.d, so optim()'s search would be
    # over a constant-NA objective -- Brent doesn't error on that, it
    # deterministically returns its `upper` bound as a value that looks
    # like a real RV but isn't one.
    if (row_is_undefined(fixed, row.only[[i]])) {
      out[i] <- NA_real_
      next
    }
    side <- ifelse(theta < conf[i,1], "lwr", "upr")
    fn <- function(rv) rv_fun(rv, dml.fit = model, par = names(out)[i],
                              side = side, theta = theta, alpha = alpha, fixed = fixed,
                              only = row.only[[i]])
    out[i] <- optim(par = c(0.01), fn, lower = 0, upper = 1, method = "Brent")$par
  }
  return(out)
  # grid <- seq(0, 0.99,by = 0.001)
  # values <- mapply(function(x,y) confidence_bounds(dml.fit, cf.y= x, cf.d = y), x = grid, y = grid)
  # rv.idx <- which(values[1,] <= theta & theta <= values[2,])[1]
  # grid[rv.idx]
}

# ##' Robustness Value DML
# ##'
# ##'@exportS3Method sensemakr::robustness_value dml.bounds
# ##'@exportS3Method dml.sensemakr::robustness_value dml.bounds
# robustness_value.dml.bounds <- function(object, theta = 0, alpha = 0.05, ...){
#   conf <- confint(object, parm = "theta.s", level = 1-alpha)
#   if(conf[,1] <= theta & theta <= conf[,2]){
#     return(0)
#   }
#   side <- ifelse(theta < conf[,1], "lwr", "upr")
#   fn <- function(rv) rv_fun(rv, dml.fit = object, side = side, theta = theta)
#   out <- optim(par = c(0.01), fn, lower=0, upper = 1, method = "Brent")$par
#   return(out)
#   # grid <- seq(0, 0.99,by = 0.001)
#   # values <- mapply(function(x,y) confidence_bounds(dml.fit, cf.y= x, cf.d = y), x = grid, y = grid)
#   # rv.idx <- which(values[1,] <= theta & theta <= values[2,])[1]
#   # grid[rv.idx]
# }
