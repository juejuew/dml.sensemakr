# Bridges DRDID::drdid_panel() output into the short.results structure that
# prep_bounds()/bounds() expect, given the raw (D, deltaY, X, w) sample that
# produced it.
#
# Unlike the did::att_gt() adapter (R/did-adapter.R), DRDID's own return
# object does not carry its data along with it -- ret$call.param is only an
# unevaluated match.call(), and ret$argu is settings, not data. So the
# caller must supply the same (D, deltaY, X, w) sample they gave to
# drdid_panel() themselves.
#
# theta.s/psi.theta.s are read directly from fit$ATT/fit$att.inf.func (only
# available if drdid_panel() was called with inffunc = TRUE). sigma2.s/
# nu2.s are computed by reusing did_cell_nuisances()/did_cell_sigma2_nu2()
# unchanged (R/did-adapter.R) -- those functions only ever operate on a
# plain (D, deltaY, X, w) sample and were already generic, not
# did::att_gt()-specific, so nothing new needed deriving for that part.
#
# Because there is no way to independently reconstruct the sample here (the
# way did_cell_sample() does from DIDparams$data), we can't verify on our
# own that the caller's (D, deltaY, X, w) actually corresponds to `fit`.
# Instead, drdid_cell_short_results() recomputes theta.s from a fresh
# nuisance refit on the supplied sample, using DRDID::drdid_panel()'s own
# dr.att formula (see its source), and errors if that doesn't reproduce
# fit$ATT closely -- this is the same three-way check (att_gt() vs. direct
# drdid_panel() vs. our own refit) that was verified by hand to match to
# machine precision on a real mpdta cell before this was written.
#
# Scope: panel data only (drdid_panel(), not drdid_rc()'s repeated
# cross-section, which has a materially different estimator internally and
# has not been derived/verified here).

# Checks that a DRDID fit is within the scope this adapter supports.
validate_drdid_fit <- function(fit) {
  if (!inherits(fit, "drdid")) {
    stop("Expected an object of class 'drdid' (the output of DRDID::drdid_panel()).")
  }
  if (!isTRUE(fit$argu$panel)) {
    stop("This adapter currently only supports panel data (argu$panel = TRUE). ",
        "Repeated-cross-section fits (e.g. from DRDID::drdid_rc()) are not yet supported.")
  }
  if (is.null(fit$att.inf.func)) {
    stop("This adapter requires the fit to have been produced with inffunc = TRUE ",
        "(e.g. DRDID::drdid_panel(..., inffunc = TRUE)).")
  }
  invisible(TRUE)
}

# Recomputes theta.s from a fresh (D, deltaY, X, w) sample and already-fit
# nuisances, replicating DRDID::drdid_panel()'s own dr.att formula exactly
# (see its source: w.treat/w.cont/eta.treat/eta.cont). Used only as a
# consistency check against fit$ATT -- not used elsewhere.
recompute_drdid_att <- function(sample, nuis) {
  D <- sample$D; deltaY <- sample$deltaY; X <- sample$X; w <- sample$w
  p_hat <- nuis$p_hat; trim_ps <- nuis$trim_ps; beta0 <- nuis$beta0

  out.delta <- as.vector(X %*% beta0)
  w.treat <- trim_ps * w * D
  w.cont  <- trim_ps * w * p_hat * (1 - D) / (1 - p_hat)
  eta.treat <- mean(w.treat * (deltaY - out.delta)) / mean(w.treat)
  eta.cont  <- mean(w.cont  * (deltaY - out.delta)) / mean(w.cont)
  eta.treat - eta.cont
}

# Assembles a short.results-shaped object (matching what ate.plm()/
# ate.npm()/did_cell_short_results() return) for a single DRDID::
# drdid_panel() fit, given the raw sample that produced it.
drdid_cell_short_results <- function(fit, D, deltaY, X, w = NULL, tol = 1e-6) {
  validate_drdid_fit(fit)

  D <- as.integer(D)
  deltaY <- as.numeric(deltaY)
  X <- as.matrix(X)
  n <- length(D)

  if (length(deltaY) != n || nrow(X) != n) {
    stop("D, deltaY, and X must all have the same number of observations.")
  }
  if (length(fit$att.inf.func) != n) {
    stop("fit$att.inf.func has length ", length(fit$att.inf.func),
        ", but the supplied sample has ", n, " observations -- the supplied ",
        "(D, deltaY, X, w) likely doesn't match the sample that produced `fit`.")
  }

  if (is.null(w)) w <- rep(1, n)
  w <- as.numeric(w) / mean(as.numeric(w))  # drdid_panel() normalizes weights the same way

  sample <- list(D = D, deltaY = deltaY, X = X, w = w)
  nuis <- did_cell_nuisances(sample)

  # Consistency check: does refitting nuisances on this exact sample and
  # recomputing theta.s reproduce fit$ATT? If not, the supplied sample does
  # not correspond to what produced `fit` -- see file header.
  att.check <- recompute_drdid_att(sample, nuis)
  if (abs(att.check - fit$ATT) > tol) {
    stop("The supplied (D, deltaY, X, w) sample does not reproduce fit$ATT when ",
        "refit (recomputed ", signif(att.check, 6), ", expected ", signif(fit$ATT, 6),
        "). Check that D, deltaY, X, and w are exactly what was passed to ",
        "drdid_panel() (including weight normalization).")
  }

  comp <- did_cell_sigma2_nu2(sample, nuis)

  theta.s <- fit$ATT
  psi.theta.s <- as.numeric(fit$att.inf.func)
  sigma2.s <- comp$sigma2.s
  nu2.s <- comp$nu2.s
  psi.S2 <- sigma2.s * comp$psi.nu2.s.cell + nu2.s * comp$psi.sigma2.s.cell

  list(
    estimates = list(theta.s = theta.s,
                     se.theta.s = psi.sd(psi.theta.s),
                     sigma2.s = sigma2.s,
                     nu2.s = nu2.s,
                     S2 = sigma2.s * nu2.s,
                     se.S2 = psi.sd(psi.S2),
                     cov.theta.S2 = mean(psi.theta.s * psi.S2) / n),
    psis = list(psi.theta.s = psi.theta.s,
               psi.sigma2.s = comp$psi.sigma2.s.cell,
               psi.nu2.s = comp$psi.nu2.s.cell,
               psi.S2 = psi.S2)
  )
}

# Assembles one or more DRDID::drdid_panel() fits into a single pseudo-
# "dml" object. `fits` must be a named list, each element itself a list
# with elements fit/D/deltaY/X/w (w optional) -- the name becomes that
# comparison's group label.
#
# Mirrors did_to_dml()'s architecture exactly: every comparison goes into
# results$groups, results$main is left NULL, so the existing
# coef.dml()/se.dml()/confint.dml()/row_only_list()/dml_bounds()/
# robustness_value()/extreme_robustness_value() machinery works unchanged.
drdid_to_dml <- function(fits) {
  if (!is.list(fits) || length(fits) == 0L || is.null(names(fits)) || any(names(fits) == "")) {
    stop("`fits` must be a non-empty named list, e.g. ",
        "list(my_comparison = list(fit = ..., D = ..., deltaY = ..., X = ...)).")
  }

  groups_results <- list()
  for (nm in names(fits)) {
    entry <- fits[[nm]]
    sr <- drdid_cell_short_results(entry$fit, entry$D, entry$deltaY, entry$X, entry$w)
    groups_results[[nm]] <- list(sr)
  }

  model <- list(
    info = list(model = "drdid", target = NULL),
    results = list(main = NULL, groups = groups_results),
    coefs = list(main = NULL, groups = lapply(groups_results, combine.cross.fits))
  )
  class(model) <- "dml"
  model
}
