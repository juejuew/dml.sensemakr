# Bridges did::att_gt() output into the short.results structure that
# prep_bounds()/bounds() expect, for a single (group, t) cell.
#
# theta.s/psi.theta.s are read directly from att_gt()'s own output
# (verified to exactly reproduce it -- see the one-cell empirical
# cross-check in project notes). sigma2.s/nu2.s are the MAIN-PAPER,
# control-population scale parameters
#
#   sigma_0s^2 := E[Var(deltaY | X, D=0) | D=0]
#   nu_0s^2    := E[(pi(X)/(1-pi(X))) / (p/(1-p)))^2 | D=0]     (p = P(D=1))
#
# -- NOT the pooled/unconditional sigma_s^2, nu_s^2 of the appendix
# parameterization that ate.plm()/ate.npm() (R/short-parameters.R) compute
# for the package's own DML fits. An earlier version of this file computed
# the pooled quantities here too (via an extra treated-arm outcome model,
# mu1(X) = E[deltaY|X,D=1]); that was a real error for this adapter's use
# case -- ovb4did's Theorem 1 (the bound this package feeds into via
# bounds()/dml_bounds()) is stated in terms of the control-population
# sigma_0s^2/nu_0s^2, not the pooled ones. See did_cell_scale_nuisances()/
# did_cell_sigma2_nu2() below for the corrected, Neyman-orthogonal,
# cross-fitted estimators, and the git history of this file for the
# removed pooled-parameterization code.
#
# Architecture (unchanged by the above correction):
#  - theta.s/psi.theta.s: external, from att_gt() (did_cell_short_results()
#    never re-estimates the ATT itself).
#  - did_cell_sample(): reconstructs the exact (D, deltaY, X, w) 2x2
#    comparison sample att_gt() used internally for one (group, t) cell,
#    from DIDparams$data alone.
#  - did_cell_scale_nuisances(): L-fold cross-fitted (pi_hat, p_hat, g_hat)
#    for the sigma_0s^2/nu_0s^2 estimating equations. Unlike theta.s's own
#    nuisances (which att_gt()/DRDID fit once on the full cell, no
#    cross-fitting -- see did_cell_nuisances() below), these are
#    genuinely cross-fitted DML nuisances, per the main paper's own
#    construction for these components.
#  - did_cell_sigma2_nu2(): solves the paper's Neyman-orthogonal
#    estimating equations for sigma_0s^2/nu_0s^2 from those nuisances, and
#    returns their influence functions -- feeding the SAME short.results
#    fields (sigma2.s, nu2.s, psi.sigma2.s.cell, psi.nu2.s.cell) the rest
#    of the package (prep_bounds()/bounds() in R/bias-bounds.R) already
#    expects, so nothing downstream of this file needed to change: S =
#    sqrt(sigma2.s*nu2.s), psi.S2 = sigma2.s*psi.nu2.s + nu2.s*psi.sigma2.s,
#    and every bound/SE bias-bounds.R derives from those two are already
#    generic in (theta.s, sigma2.s, nu2.s, their psis) and don't assume
#    anything about how sigma2.s/nu2.s were estimated.
#  - did_cell_nuisances()/recompute_drdid_att() (drdid-adapter.R): a
#    SEPARATE, non-cross-fitted, full-cell (propensity + control-arm
#    outcome only, no treated-arm model) refit used ONLY to verify that a
#    caller-supplied sample reproduces fit$ATT -- this is a consistency
#    check on theta.s, not part of the sigma_0s^2/nu_0s^2 estimator.
#
# Scope, and what has actually been verified so far:
#  - panel data, est_method = "dr", no clustering, fix_weights unset --
#    validate_did_scope() enforces these and errors otherwise.
#  - unit weights only (all w == 1) for the cross-fitted sigma_0s^2/
#    nu_0s^2 estimator -- validate_unit_weights() enforces this and errors
#    otherwise. theta.s/psi.theta.s (external, from att_gt()/DRDID) are
#    unaffected and work with whatever weights att_gt()/drdid_panel() were
#    themselves given.
#  - control_group = "nevertreated": base-period resolution, treated/
#    control filtering, and the resulting theta.s/psi.theta.s were
#    verified to exactly reproduce att_gt()'s own numbers on a real
#    (mpdta) cell.
#  - control_group = "notyettreated" and base_period = "universal" (and
#    anticipation > 0) have ALSO been verified the same way -- exact
#    reproduction of att_gt()'s own att/se on real mpdta cells for each.
#  - The (n/n1) rescaling and zero-padding used to embed a cell's own
#    psi.sigma2.s.cell/psi.nu2.s.cell into the full mp$n-length vector
#    (did_cell_short_results(), below) was cross-checked directly against
#    did's own source (run_DRDID(): "if_x <- (n/n1) * attgt$att.inf.func"
#    for panel data) and against how did lays out unit rows: data is
#    sorted by (tname, gname, idname) in did_standardization(), then
#    time_invariant_data <- data[1:id_count] and cohort/control index
#    blocks (get_did_cohort_index()) are built on that SAME sorted table
#    -- so row i of DIDparams$time_invariant_data really is unit i of
#    mp$inffunc's index space, and match(sample$ids,
#    time_invariant_data[[idname]]) (used below) look up that row by ID,
#    not by an assumed row-order coincidence.
#
# sensemakr() on did_to_dml()'s output (results$main = NULL, groups-only):
#  - sensemakr(), print(<dml.sensemakr>), summary(<dml.sensemakr>),
#    robustness_value(), extreme_robustness_value(), and confidence_bounds()
#    all work as expected. summary.dml()/print.summary_dml()
#    (R/print-summary-dml.R) now detect the absence of object$fits (no real
#    cross-fitting happened) and skip the ML-method/tuning/R^2 reporting,
#    showing only the group (g,t) ATE table.
#  - plot(<dml.sensemakr>) with default args errors clearly, pointing at
#    group = TRUE, group.number = k (it needs a "main" target to plot by
#    default, which this adapter doesn't have -- there's no single natural
#    default across a grid of (g,t) cells). plot(group = TRUE,
#    group.number = k) works and plots that specific cell.
#  - dml_benchmark()/benchmark_covariates are NOT available at all for this
#    adapter's output, at any level: bench_fun() (R/benchmarks.R) reads
#    exclusively from model$results$main, which this adapter always leaves
#    NULL (results live in model$results$groups instead). sensemakr(...,
#    benchmark_covariates = ...) detects this and skips with a warning
#    rather than erroring.

# Checks that a did::att_gt() "MP" object is within the scope this adapter
# supports, erroring clearly and immediately rather than silently computing
# something built on unverified assumptions.
validate_did_scope <- function(mp) {
  if (!inherits(mp, "MP")) {
    stop("Expected an object of class 'MP' (the output of did::att_gt()).")
  }
  dp <- mp$DIDparams
  if (!isTRUE(dp$panel)) {
    stop("This adapter currently only supports panel data (DIDparams$panel = TRUE). ",
        "Repeated-cross-section fits are not yet supported.")
  }
  if (!identical(dp$est_method, "dr")) {
    stop("This adapter requires att_gt() to have been fit with est_method = \"dr\" ",
        "(doubly robust). Got: '", dp$est_method, "'. The 'reg' and 'ipw' variants ",
        "don't fit both nuisance models, so nu2.s/sigma2.s can't be recovered.")
  }
  if (!is.null(dp$clustervars) && length(dp$clustervars) > 0) {
    stop("This adapter does not yet support clustered standard errors ",
        "(DIDparams$clustervars is set). The sigma2.s/nu2.s influence functions ",
        "computed here are not cluster-robust.")
  }
  if (!is.null(dp$fix_weights)) {
    stop("This adapter does not yet support fix_weights (DIDparams$fix_weights is set). ",
        "The weight-period resolution implemented here assumes the default ",
        "(w_idx = min(pret, t)) rule.")
  }
  if (is.null(mp$inffunc)) {
    stop("This adapter requires att_gt() to have been called with compute_inffunc = TRUE.")
  }
  invisible(TRUE)
}

# Resolves the base (pre-treatment) period for a single (group, t) cell,
# following run_att_gt_estimation()'s exact logic:
#   - base_period == "universal": always the group's own last
#     pre-treatment period (anticipation-adjusted), regardless of t.
#   - base_period == "varying": the period immediately preceding t, UNLESS
#     t is already post-treatment for this group (group <= t), in which
#     case it also falls back to the group's own last pre-treatment period.
resolve_did_base_period <- function(dp, group, t) {
  time_periods <- dp$time_periods
  anticipation <- dp$anticipation

  pre_idx <- which((time_periods + anticipation) < group)
  if (length(pre_idx) == 0L) {
    stop("No pre-treatment period available for group '", group, "'.")
  }
  pret_g <- time_periods[pre_idx[length(pre_idx)]]

  if (identical(dp$base_period, "universal") || group <= t) {
    return(pret_g)
  }

  t_idx <- match(t, time_periods)
  if (is.na(t_idx) || t_idx < 2L) {
    stop("Cannot resolve a preceding period for t = ", t, " under base_period = 'varying'.")
  }
  time_periods[t_idx - 1L]
}

# Resolves the control-group unit ids for a single (group, t, pret) cell,
# following get_did_cohort_index()'s exact logic. "nevertreated": only
# units with gname == Inf (did's own internal never-treated encoding,
# regardless of the raw input's own sentinel value). "notyettreated":
# additionally includes units whose own first-treatment period is strictly
# later than max(t, pret) + anticipation.
#
# Both branches have been independently validated by exact reproduction of
# att_gt()'s own att/se on real mpdta cells (see file header).
resolve_did_control_ids <- function(dp, group, t, pret) {
  d <- dp$data
  idname <- dp$idname; gname <- dp$gname

  if (!identical(dp$control_group, "notyettreated")) {
    return(unique(d[[idname]][d[[gname]] == Inf]))
  }

  time_periods <- dp$time_periods
  tfac <- if (identical(dp$base_period, "universal")) 0L else 1L
  m_idx <- match(max(t, pret), time_periods)
  cutoff_idx <- m_idx + tfac
  cutoff <- if (cutoff_idx > length(time_periods)) {
    max(time_periods)  # no later observed period -- only never-treated (Inf) can qualify
  } else {
    time_periods[cutoff_idx] + dp$anticipation
  }
  unique(d[[idname]][d[[gname]] > cutoff])
}

# Reconstructs the (D, deltaY, X, w) sample for one (group, t) cell of a
# did::att_gt() "MP" object, from DIDparams$data alone -- verified to
# exactly reproduce att_gt()'s own ATT/SE for the "nevertreated" case (see
# file header). Weights are read from period min(pret, t), matching
# run_att_gt_estimation()'s own w_idx <- min(pret, t) (the default,
# fix_weights = NULL, rule -- see validate_did_scope()).
did_cell_sample <- function(mp, group, t) {
  dp <- mp$DIDparams
  d  <- dp$data
  idname <- dp$idname; tname <- dp$tname; gname <- dp$gname; yname <- dp$yname

  pret <- resolve_did_base_period(dp, group, t)
  treated_ids <- unique(d[[idname]][d[[gname]] == group])
  control_ids <- resolve_did_control_ids(dp, group, t, pret)
  keep_ids <- c(treated_ids, control_ids)

  d_t <- d[d[[tname]] == t & d[[idname]] %in% keep_ids, ]
  d_0 <- d[d[[tname]] == pret & d[[idname]] %in% keep_ids, ]
  d_t <- d_t[match(keep_ids, d_t[[idname]]), ]
  d_0 <- d_0[match(keep_ids, d_0[[idname]]), ]

  if (anyNA(d_t[[idname]]) || anyNA(d_0[[idname]])) {
    stop("Unbalanced panel or missing observations for group '", group,
        "' at periods ", pret, "/", t, " -- not currently supported.")
  }

  w_period <- min(pret, t)
  d_w <- if (w_period == t) d_t else if (w_period == pret) d_0 else {
    dw <- d[d[[tname]] == w_period & d[[idname]] %in% keep_ids, ]
    dw[match(keep_ids, dw[[idname]]), ]
  }

  D <- as.integer(d_t[[idname]] %in% treated_ids)
  deltaY <- d_t[[yname]] - d_0[[yname]]
  X <- model.matrix(dp$xformla, data = d_0)
  w <- d_w[[".w"]]

  list(D = D, deltaY = deltaY, X = X, w = w, ids = keep_ids,
      group = group, t = t, pret = pret, n1 = length(keep_ids))
}

# Fits the propensity score (weighted logit) and the control-arm outcome-
# change regression (mu0) on a single cell's reconstructed sample. Mirrors
# drdid_panel()'s own fastglm_fit() calls (full sample, no cross-fitting),
# via base R's glm.fit()/solve() -- verified numerically equivalent for the
# validated cell (see file header). USED ONLY as the consistency check in
# recompute_drdid_att() (drdid-adapter.R) that a caller-supplied sample
# reproduces fit$ATT -- NOT used to estimate sigma2.s/nu2.s (see
# did_cell_scale_nuisances()/did_cell_sigma2_nu2() for that). There is
# therefore no treated-arm outcome model (mu1) here: recompute_drdid_att()
# never needed one (DR-DiD's own att formula doesn't use it), and the old
# pooled sigma2.s that did need one has been removed (see file header).
#
# `trim.level` defaults to drdid_panel()'s own default (0.995); pass the
# actual trim.level a given fit was produced with when it might differ
# (drdid-adapter.R does, reading it from fit$argu$trim.level) -- otherwise
# recompute_drdid_att()'s consistency check would replicate the WRONG
# trimming rule and could spuriously fail (or, worse, spuriously pass by
# coincidence) for a fit produced with a non-default trim.level.
did_cell_nuisances <- function(sample, trim.level = 0.995) {
  D <- sample$D; deltaY <- sample$deltaY; X <- sample$X; w <- sample$w
  n <- length(D)

  fit_ps <- glm.fit(x = X, y = D, weights = w, family = binomial())
  p_hat  <- pmin(fit_ps$fitted.values, 1 - 1e-6)

  # trimming -- needed here so that recompute_drdid_att() replicates
  # drdid_panel()'s own att formula exactly, including its trimming.
  # Treated units are never trimmed, matching drdid_panel().
  trim_ps <- rep(1, n)
  trim_ps[D == 0] <- as.numeric(p_hat[D == 0] < trim.level)

  beta0 <- as.vector(solve(
    crossprod(X[D == 0, , drop = FALSE], w[D == 0] * X[D == 0, , drop = FALSE]),
    crossprod(X[D == 0, , drop = FALSE], w[D == 0] * deltaY[D == 0])
  ))

  list(p_hat = p_hat, trim_ps = trim_ps, beta0 = beta0)
}

# Errors if propensity-score trimming is ACTUALLY active for this cell's
# analysis population -- i.e. at least one control unit's fitted propensity
# is >= trim.level, so drdid_panel() (called directly, or internally by
# did::att_gt(est_method = "dr")) excludes it from the external ATT via its
# own trim.ps[D==0] <- (ps.fit[D==0] < trim.level) rule (see DRDID::
# drdid_panel()'s source). When that happens, theta.s/psi.theta.s (external,
# defined on the TRIMMED sample) and this adapter's sigma2.s/nu2.s (defined
# on the FULL reconstructed cell -- did_cell_scale_nuisances() does not
# implement any trimming, see file header "Scope") would no longer target
# the same analysis population -- the external-estimator compatibility
# assumption this whole adapter relies on (see file header) would be
# violated. Paper-consistent sensitivity analysis for a trimmed ATT target
# has not been derived, so such cells are rejected rather than silently
# computed on mismatched populations.
#
# trim.level < 1 alone is NOT grounds for rejection -- drdid_panel()'s
# default (0.995) is virtually always set, but in most datasets no
# observation actually crosses it. Only ACTIVE trimming (the indicator
# actually firing for >= 1 control) is rejected.
validate_no_active_trimming <- function(D, p_hat, trim.level) {
  n_trimmed <- sum(D == 0 & p_hat >= trim.level)
  if (n_trimmed > 0) {
    stop(n_trimmed, " control unit(s) have an estimated propensity score ",
        ">= trim.level (", trim.level, ") and are therefore TRIMMED by the ",
        "external DR-DiD ATT estimator (DRDID::drdid_panel(), used directly ",
        "or internally by did::att_gt(est_method = \"dr\")). When trimming ",
        "is active, the external theta.s and this adapter's sigma2.s/nu2.s ",
        "no longer correspond to the same analysis population -- ",
        "paper-consistent sensitivity analysis for a trimmed ATT target has ",
        "not yet been derived/supported. This cell cannot be used with this ",
        "adapter.")
  }
  invisible(TRUE)
}

# Errors if any weight differs from 1 -- the cross-fitted sigma_0s^2/
# nu_0s^2 estimating equations below have only been derived for unit
# weights (see file header, "Scope").
validate_unit_weights <- function(w) {
  if (!isTRUE(all.equal(unname(w), rep(1, length(w)), tolerance = 1e-8, check.attributes = FALSE))) {
    stop("The cross-fitted sigma2.s/nu2.s (sigma_0s^2/nu_0s^2) estimator ",
        "currently only supports unit weights (all w == 1). Non-unit ",
        "sampling weights have not been separately derived for this ",
        "estimator -- see R/did-adapter.R.")
  }
  invisible(TRUE)
}

# Assigns each of n1 cell observations to one of L folds for the
# cross-fitted scale-nuisance estimation below. A single common partition
# is used for all three nuisances (g_0s, pi, p), per the main paper's
# construction.
make_cv_folds <- function(n1, L, seed) {
  if (!is.null(seed)) set.seed(seed)
  sample(rep_len(seq_len(L), n1))
}

# L-fold cross-fitted nuisances for the paper's control-population scale
# estimators (see file header): g_0s(X) = E[deltaY|X,D=0] (weighted-least-
# squares outcome regression, fit ONLY on D=0 units of the training fold),
# pi(X) = P(D=1|X) (weighted logistic regression, fit on the full training
# fold), and p_l = mean(D) in the training fold (the fold's own estimate of
# the marginal treatment probability p = P(D=1)). Every held-out
# observation is scored using nuisances trained on the OTHER folds only.
#
# Why the marginal p_l (not the conditional pi_hat(X)) belongs in
# sigma_0s^2/nu_0s^2's estimating equations (did_cell_sigma2_nu2(), below):
# for any h(X), E[(1-D)/(1-p) * h(X)] = E_X[h(X) * (1-pi(X))] / (1-p) =
# E[h(X) | D=0] by Bayes' rule (f(X|D=0) = f(X)*(1-pi(X))/(1-p)) -- i.e.
# using the MARGINAL p as the reweighting denominator is exactly what
# converts a whole-population moment into the control-population moment
# the paper's sigma_0s^2/nu_0s^2 are defined as. Using the CONDITIONAL
# pi_hat(X) there instead would collapse the reweighting to the trivial
# identity E[(1-D)/(1-pi(X))*h(X)] = E_X[h(X)] -- the WHOLE-population
# mean of h(X), not the control-population one. pi_hat(X) is still needed,
# but only inside r_i = odds-ratio(X_i) itself (did_cell_sigma2_nu2()).
did_cell_scale_nuisances <- function(sample, L = 5, seed = 1) {
  D <- sample$D; deltaY <- sample$deltaY; X <- sample$X; w <- sample$w
  n1 <- length(D)
  validate_unit_weights(w)

  fold <- make_cv_folds(n1, L, seed)

  pi_hat <- rep(NA_real_, n1)
  p_hat  <- rep(NA_real_, n1)
  g_hat  <- rep(NA_real_, n1)

  for (l in seq_len(L)) {
    train <- which(fold != l)
    held  <- which(fold == l)

    D_tr <- D[train]; X_tr <- X[train, , drop = FALSE]; dY_tr <- deltaY[train]

    if (sum(D_tr) < 1L || sum(D_tr) > length(D_tr) - 1L) {
      stop("Fold ", l, " of the cross-fitted scale-estimation split has no ",
          "variation in D in its training set (all treated or all control) ",
          "-- cannot fit a propensity model. Try a smaller `cf.folds`.")
    }
    X0_tr <- X_tr[D_tr == 0, , drop = FALSE]
    if (nrow(X0_tr) <= ncol(X0_tr)) {
      stop("Fold ", l, " of the cross-fitted scale-estimation split has too ",
          "few control units in its training set to fit g_0s(X) (need more ",
          "than ", ncol(X0_tr), " control observations). Try a smaller ",
          "`cf.folds`.")
    }

    fit_pi  <- glm.fit(x = X_tr, y = D_tr, family = binomial())
    eta_hld <- as.vector(X[held, , drop = FALSE] %*% fit_pi$coefficients)
    pi_hat[held] <- 1 / (1 + exp(-eta_hld))

    p_hat[held] <- mean(D_tr)

    dY0_tr <- dY_tr[D_tr == 0]
    beta_g <- as.vector(solve(crossprod(X0_tr), crossprod(X0_tr, dY0_tr)))
    g_hat[held] <- as.vector(X[held, , drop = FALSE] %*% beta_g)
  }

  eps <- 1e-6
  pi_hat <- pmin(pmax(pi_hat, eps), 1 - eps)
  p_hat  <- pmin(pmax(p_hat, eps), 1 - eps)

  list(pi_hat = pi_hat, p_hat = p_hat, g_hat = g_hat, fold = fold)
}

# Solves the paper's Neyman-orthogonal estimating equations for
# sigma_0s^2/nu_0s^2 from cross-fitted nuisances (did_cell_scale_nuisances()
# above), and returns their influence functions -- at the SINGLE cell's own
# sample size (rescaled to full-sample by the caller). Output field names
# match the pooled-parameterization code this replaces (sigma2.s, nu2.s,
# psi.sigma2.s.cell, psi.nu2.s.cell), so downstream embedding/rescaling
# (did_cell_short_results()) and every consumer in R/bias-bounds.R are
# unchanged. See file header for what's verified and the derivation.
did_cell_sigma2_nu2 <- function(sample, nuis) {
  D <- sample$D; deltaY <- sample$deltaY
  pi_hat <- nuis$pi_hat; p_hat <- nuis$p_hat; g_hat <- nuis$g_hat

  # r_i = odds-ratio(X_i) = [pi(X_i)/(1-pi(X_i))] / [p/(1-p)]
  r <- (pi_hat / (1 - pi_hat)) / (p_hat / (1 - p_hat))

  # --- sigma_0s^2 := E[Var(deltaY | X, D=0) | D=0] ---
  a_sigma <- -(1 - D) / (1 - p_hat)
  b_sigma <- (1 - D) / (1 - p_hat) * (deltaY - g_hat)^2
  mean_a_sigma <- mean(a_sigma)
  sigma2_hat <- -mean(b_sigma) / mean_a_sigma
  psi_sigma_hat <- a_sigma * sigma2_hat + b_sigma
  phi_sigma <- -psi_sigma_hat / mean_a_sigma

  # --- nu_0s^2 := E[(pi(X)/(1-pi(X)) / (p/(1-p)))^2 | D=0] ---
  a_nu <- -2 * D / p_hat + (1 - D) / (1 - p_hat)
  b_nu <- 2 * D / p_hat * r - (1 - D) / (1 - p_hat) * r^2
  mean_a_nu <- mean(a_nu)  # -> -1 asymptotically, but NOT hard-coded: see
                            # item 5 of the request this implements -- p_hat
                            # is fold-estimated, so this need not equal -1
                            # exactly in finite samples.
  nu2_hat <- -mean(b_nu) / mean_a_nu
  psi_nu_hat <- a_nu * nu2_hat + b_nu
  phi_nu <- -psi_nu_hat / mean_a_nu

  list(sigma2.s = sigma2_hat, nu2.s = nu2_hat,
      psi.sigma2.s.cell = phi_sigma, psi.nu2.s.cell = phi_nu)
}

# Assembles a short.results-shaped object (matching what ate.plm()/
# ate.npm() return) for a single (group, t) cell of a did::att_gt() "MP"
# object, ready to feed into prep_bounds()/bounds(). theta.s/psi.theta.s
# are read directly from att_gt()'s own output; sigma2.s/nu2.s/their psis
# are the cross-fitted control-population estimators (did_cell_scale_
# nuisances()/did_cell_sigma2_nu2(), see file header) rescaled to the same
# full-sample, zero-padded convention att.inf.func already uses (see
# did:::run_DRDID()'s own (n/n1)-rescaling -- verified directly against
# did's source, see file header).
#
# `cf.folds`/`cf.seed` control the L-fold cross-fitting used for sigma2.s/
# nu2.s only (theta.s/psi.theta.s are untouched, since they are external).
did_cell_short_results <- function(mp, group, t, cf.folds = 5, cf.seed = 1) {
  validate_did_scope(mp)

  k <- which(mp$group == group & mp$t == t)
  if (length(k) != 1L) {
    stop("Could not find a unique (group, t) = (", group, ", ", t,
        ") cell in this att_gt() object.")
  }

  theta.s <- unname(mp$att[k])
  psi.theta.s.full <- as.numeric(as.matrix(mp$inffunc)[, k])
  n <- mp$n

  sample <- did_cell_sample(mp, group, t)

  # did::att_gt() never exposes a trim.level argument -- run_DRDID() always
  # calls drdid_panel() (for est_method = "dr", the only est_method this
  # adapter supports) without one, so it always uses DRDID's own hardcoded
  # default (0.995), confirmed against did's source. Refit the same
  # full-sample (non-cross-fitted) propensity model drdid_panel() itself
  # would fit, purely to check whether that trimming is actually active for
  # this cell -- see validate_no_active_trimming().
  nuis_val <- did_cell_nuisances(sample, trim.level = 0.995)
  validate_no_active_trimming(sample$D, nuis_val$p_hat, trim.level = 0.995)

  nuis   <- did_cell_scale_nuisances(sample, L = cf.folds, seed = cf.seed)
  comp   <- did_cell_sigma2_nu2(sample, nuis)

  scale <- n / sample$n1
  idx <- match(sample$ids, mp$DIDparams$time_invariant_data[[mp$DIDparams$idname]])
  psi.sigma2.s <- rep(0, n)
  psi.nu2.s    <- rep(0, n)
  psi.sigma2.s[idx] <- scale * comp$psi.sigma2.s.cell
  psi.nu2.s[idx]    <- scale * comp$psi.nu2.s.cell

  sigma2.s <- comp$sigma2.s
  nu2.s    <- comp$nu2.s
  psi.S2   <- sigma2.s * psi.nu2.s + nu2.s * psi.sigma2.s

  list(
    estimates = list(theta.s = theta.s,
                     se.theta.s = psi.sd(psi.theta.s.full),
                     sigma2.s = sigma2.s,
                     nu2.s = nu2.s,
                     S2 = sigma2.s * nu2.s,
                     se.S2 = psi.sd(psi.S2),
                     cov.theta.S2 = mean(psi.theta.s.full * psi.S2) / n),
    psis = list(psi.theta.s = psi.theta.s.full,
               psi.sigma2.s = psi.sigma2.s,
               psi.nu2.s = psi.nu2.s,
               psi.S2 = psi.S2)
  )
}

# Slot name for one (group, t) cell, e.g. "g2004_t2006".
did_slot_name <- function(group, t) paste0("g", group, "_t", t)

# Loops over every (group, t) cell of a did::att_gt() "MP" object and
# assembles them into a single pseudo-"dml" object, ready for
# confidence_bounds()/robustness_value()/extreme_robustness_value().
#
# Every cell goes into results$groups/coefs$groups; results$main/
# coefs$main are left NULL, and info$target is left NULL with info$model
# set to "did" (anything other than "plm"). This was verified (not just
# assumed) to work correctly with the existing coef.dml()/se.dml()/
# confint.dml()/row_only_list()/dml_bounds()/robustness_value()/
# extreme_robustness_value() machinery unchanged -- see the project notes
# for the verification, which also uncovered and fixed a real,
# previously-untriggered bug in coef.dml()/se.dml() (sapply(NULL, f)
# returns list(), silently coercing the whole coefficient vector into a
# list when object$coefs$main is NULL).
#
# Cells where att_gt() itself returned NA -- either because of an overlap/
# rank-condition failure (att = NA), or because it's the trivial
# "comparison period equals the group's own base period" cell that
# base_period = "universal" produces for every group (att = 0 exactly, but
# se = NA, since att_gt() special-cases it as a self-comparison and never
# actually calls a DiD estimator for it -- verified against a real
# att_gt() cell, see project notes) -- are skipped with a warning, not
# silently dropped. Cells where this adapter's own nuisance refit fails
# (e.g. a singular design matrix in a very small subgroup) are also
# skipped with a warning, rather than aborting the whole grid.
did_to_dml <- function(mp, cf.folds = 5, cf.seed = 1) {
  validate_did_scope(mp)

  groups_results <- list()
  for (k in seq_along(mp$group)) {
    group <- mp$group[k]; t <- mp$t[k]
    slot <- did_slot_name(group, t)

    if (is.na(mp$att[k]) || is.na(mp$se[k])) {
      warning("Skipping (group, t) = (", group, ", ", t, "): att_gt() itself ",
              "returned att = NA and/or se = NA for this cell (e.g. an overlap/",
              "rank-condition failure, or a trivial self-comparison cell under ",
              "base_period = \"universal\").",
              call. = FALSE)
      next
    }

    sr <- tryCatch(did_cell_short_results(mp, group, t, cf.folds = cf.folds, cf.seed = cf.seed), error = function(e) {
      warning("Skipping (group, t) = (", group, ", ", t, "): ", conditionMessage(e),
              call. = FALSE)
      NULL
    })
    if (is.null(sr)) next

    groups_results[[slot]] <- list(sr)
  }

  if (length(groups_results) == 0L) {
    stop("No (group, t) cells could be computed.")
  }

  model <- list(
    info = list(model = "did", target = NULL),
    results = list(main = NULL, groups = groups_results),
    coefs = list(main = NULL, groups = lapply(groups_results, combine.cross.fits))
  )
  class(model) <- "dml"
  model
}
