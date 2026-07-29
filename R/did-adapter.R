# Bridges did::att_gt() output into the short.results structure that
# prep_bounds()/bounds() expect, for a single (group, t) cell.
#
# theta.s/psi.theta.s are read directly from att_gt()'s own output
# (verified to exactly reproduce it -- see the one-cell empirical
# cross-check in project notes). sigma2.s/nu2.s and their influence
# functions are computed here, following the same non-cross-fitted
# M-estimation construction did/DRDID use internally for theta.s (nuisances
# fit once on the full 2x2 comparison sample), extended with a new mu1(X)
# fit (the treated-arm outcome-change regression, which DR-DiD's own
# estimator never needs but sigma2.s's "each unit's own arm" definition
# does).
#
# Scope, and what has actually been verified so far:
#  - panel data, est_method = "dr", no clustering, fix_weights unset --
#    validate_did_scope() enforces these and errors otherwise.
#  - control_group = "nevertreated": base-period resolution, treated/
#    control filtering, and the resulting theta.s/psi.theta.s were
#    verified to exactly reproduce att_gt()'s own numbers on a real
#    (mpdta) cell. sigma2.s needs no sandwich correction (verified via
#    bootstrap: an envelope-theorem argument, since beta0/beta1 solve
#    their own weighted-least-squares first-order condition). nu2.s's
#    sandwich correction (for gamma/c1/c0's estimation uncertainty) was
#    verified via bootstrap on both an unweighted and a weighted
#    simulated sample.
#  - control_group = "notyettreated" and base_period = "universal" (and
#    anticipation > 0) have ALSO been verified the same way -- exact
#    reproduction of att_gt()'s own att/se on real mpdta cells for each.
#  - Trimming (trim.level = 0.995, matching drdid_panel()'s default): the
#    nu2.s sandwich correction does NOT differentiate through the trim
#    indicator (standard argument in this literature -- Sant'Anna & Zhao,
#    2020; Crump, Hotz, Imbens & Mitnik, 2009 -- that the set of trimmed
#    units doesn't change under small perturbations of gamma, as long as
#    the true propensity has no probability mass exactly at the
#    threshold). This WAS checked by simulation (two independent
#    bootstrap comparisons, one with a single boundary-case trimmed unit,
#    one with a larger n and a stable set of ~4 trimmed units): the
#    corrected formula captures the great majority of the true SE
#    (~95% of the bootstrap SD, vs. ~30-40% for the naive, uncorrected
#    one) but leaves a small, consistent residual understatement of
#    around 5% when trimming is genuinely active -- attributable to
#    exactly that omitted derivative. See row_is_undefined()-style
#    warnings elsewhere in this package for the precedent this follows:
#    did_cell_nuisances() warns when trimming is active on a cell, so
#    users know that cell's SE may be modestly understated.

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
# NOTE: only the "nevertreated" branch has been independently validated
# (see file header); "notyettreated" is traced from source but unverified.
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

# Fits the propensity score (weighted logit) and both outcome-change
# regressions (control-arm mu0, matching drdid_panel()'s own fit; and the
# new treated-arm mu1) on a single cell's reconstructed sample. Mirrors
# drdid_panel()'s own fastglm_fit() calls (full sample, no cross-fitting),
# via base R's glm.fit()/solve() -- verified numerically equivalent for the
# validated cell (see file header).
did_cell_nuisances <- function(sample) {
  D <- sample$D; deltaY <- sample$deltaY; X <- sample$X; w <- sample$w
  n <- length(D)

  # nu2.s's sandwich correction includes a 1/c1^3 term (c1 = mean of the
  # treated-arm weight), which blows up as the treated group shrinks --
  # confirmed by simulation (see file header) to produce a real, roughly
  # 12% analytical-vs-bootstrap SE gap on a real cell with only 20 treated
  # units, resolved on a much larger, otherwise-identical sample. The
  # threshold here (<= 20) is that one confirmed problem case (exactly 20
  # treated units), not a precisely calibrated cutoff -- it's a heuristic,
  # not a guarantee that 21+ treated units is always fine or that 20
  # always fails.
  n_treated <- sum(D)
  if (n_treated <= 20) {
    warning("Only ", n_treated, " treated unit(s) in this cell. nu2.s's ",
            "sandwich correction includes a term that grows as the treated ",
            "group shrinks, and was found, by simulation, to meaningfully ",
            "understate the true SE on a similarly small cell (see ",
            "R/did-adapter.R) -- treat this cell's nu2.s-derived precision ",
            "with extra caution.", call. = FALSE)
  }

  fit_ps <- glm.fit(x = X, y = D, weights = w, family = binomial())
  p_hat  <- pmin(fit_ps$fitted.values, 1 - 1e-6)

  # trimming (trim.level = 0.995, matching drdid_panel()'s default); see
  # file header for what this is verified to cost. Treated units are
  # never trimmed, matching drdid_panel().
  trim_ps <- rep(1, n)
  trim_ps[D == 0] <- as.numeric(p_hat[D == 0] < 0.995)
  n_trimmed <- sum(trim_ps[D == 0] == 0)
  if (n_trimmed > 0) {
    warning(n_trimmed, " control unit(s) trimmed (propensity score >= 0.995) ",
            "for this cell. nu2.s's standard error does not account for the ",
            "trimming indicator's own estimation uncertainty and was found, ",
            "by simulation, to understate the true SE by roughly 5% when ",
            "trimming is active (see R/did-adapter.R) -- treat this cell's ",
            "precision with a little extra caution.", call. = FALSE)
  }

  beta0 <- as.vector(solve(
    crossprod(X[D == 0, , drop = FALSE], w[D == 0] * X[D == 0, , drop = FALSE]),
    crossprod(X[D == 0, , drop = FALSE], w[D == 0] * deltaY[D == 0])
  ))
  beta1 <- as.vector(solve(
    crossprod(X[D == 1, , drop = FALSE], w[D == 1] * X[D == 1, , drop = FALSE]),
    crossprod(X[D == 1, , drop = FALSE], w[D == 1] * deltaY[D == 1])
  ))

  list(p_hat = p_hat, trim_ps = trim_ps, beta0 = beta0, beta1 = beta1)
}

# Computes sigma2.s/psi.sigma2.s and nu2.s/psi.nu2.s (both at the SINGLE
# cell's own sample size -- rescaled to full-sample by the caller) for one
# (group, t) cell. See file header for what's verified.
did_cell_sigma2_nu2 <- function(sample, nuis) {
  D <- sample$D; deltaY <- sample$deltaY; X <- sample$X; w <- sample$w
  n <- length(D)
  p_hat <- nuis$p_hat; trim_ps <- nuis$trim_ps
  beta0 <- nuis$beta0; beta1 <- nuis$beta1

  # --- sigma2.s: no sandwich correction needed (envelope theorem) ---
  gs <- as.vector(ifelse(D == 1, X %*% beta1, X %*% beta0))
  sigma2_hat <- mean(w * (deltaY - gs)^2)
  psi_sigma2_cell <- w * (deltaY - gs)^2 - sigma2_hat

  # --- nu2.s: sandwich correction needed for gamma/c1/c0 ---
  w1 <- w * D * trim_ps
  w0 <- w * (1 - D) * p_hat / (1 - p_hat) * trim_ps
  c1 <- mean(w1); c0 <- mean(w0)
  a  <- w1 / c1 - w0 / c0
  nu2_hat <- mean(a^2)
  psi_nu2_naive <- a^2 - nu2_hat

  score_mat <- (w * (D - p_hat)) * X
  Wd <- w * p_hat * (1 - p_hat)
  Hbar <- crossprod(X, X * Wd) / n
  IF_gamma <- score_mat %*% solve(Hbar)

  Ew1sq <- mean(w1^2)
  Ew0sq <- mean(w0^2)
  Ew0X  <- colMeans(w0 * X)

  IF_c0 <- (w0 - c0) + as.vector(IF_gamma %*% Ew0X)
  corr_gamma <- (2 / c0^2) * as.vector(IF_gamma %*% colMeans(w0^2 * X))
  corr_c1    <- -(2 / c1^3) * Ew1sq * (w1 - c1)
  corr_c0    <- -(2 / c0^3) * Ew0sq * IF_c0

  psi_nu2_cell <- psi_nu2_naive + corr_gamma + corr_c1 + corr_c0

  list(sigma2.s = sigma2_hat, nu2.s = nu2_hat,
      psi.sigma2.s.cell = psi_sigma2_cell, psi.nu2.s.cell = psi_nu2_cell)
}

# Assembles a short.results-shaped object (matching what ate.plm()/
# ate.npm() return) for a single (group, t) cell of a did::att_gt() "MP"
# object, ready to feed into prep_bounds()/bounds(). theta.s/psi.theta.s
# are read directly from att_gt()'s own output; sigma2.s/nu2.s/their psis
# are computed above and rescaled to the same full-sample, zero-padded
# convention att.inf.func already uses (see did:::run_DRDID()'s own
# (n/n1)-rescaling, verified to reproduce the same SE as computing directly
# on the cell's own n1-length psi).
did_cell_short_results <- function(mp, group, t) {
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
  nuis   <- did_cell_nuisances(sample)
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
did_to_dml <- function(mp) {
  validate_did_scope(mp)

  if (identical(mp$DIDparams$control_group, "notyettreated")) {
    warning("control_group = \"notyettreated\" is implemented from did's own ",
            "source logic but has not been independently verified the same way ",
            "as \"nevertreated\" (see R/did-adapter.R). Treat results with extra caution.",
            call. = FALSE)
  }

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

    sr <- tryCatch(did_cell_short_results(mp, group, t), error = function(e) {
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
