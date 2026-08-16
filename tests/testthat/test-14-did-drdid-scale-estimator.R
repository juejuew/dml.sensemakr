# Tests for the corrected, paper-native, cross-fitted control-population
# scale estimator shared by the did/DRDID adapters (did_cell_scale_
# nuisances()/did_cell_sigma2_nu2(), R/did-adapter.R): sigma2.s/nu2.s now
# target
#
#   sigma_0s^2 := E[Var(deltaY | X, D=0) | D=0]
#   nu_0s^2    := E[(pi(X)/(1-pi(X)) / (p/(1-p)))^2 | D=0]
#
# (the ovb4did main-paper, control-population parameterization), NOT the
# pooled/whole-population quantities an earlier version of this adapter
# computed. These tests exercise did_cell_scale_nuisances()/
# did_cell_sigma2_nu2() directly on synthetic (D, deltaY, X, w) samples with
# a known DGP -- independent of did::att_gt()/DRDID::drdid_panel() -- so the
# true sigma_0s^2/nu_0s^2 (and a clearly WRONG pooled/whole-population
# alternative) can be pinned down exactly and compared to the estimator's
# output.
skip_if_not_installed("DRDID")

# --- Shared DGP for the target tests (A, B, C) -----------------------------
# X = (1, x1), x1 ~ Uniform(-1, 1); pi(X) = plogis(g0 + g1*x1); D ~
# Bernoulli(pi(X)); deltaY = b0 + b1*x1 + tau*D + eps, eps ~ N(0, sigma1^2)
# if D=1, N(0, sigma0^2) if D=0 -- i.e. HETEROSKEDASTIC BY ARM
# (sigma1 != sigma0), so sigma_0s^2 = sigma0^2 exactly (homoskedastic
# within the control arm, so it doesn't depend on X), while the pooled
# "each unit's own arm" quantity the old estimator computed would be close
# to p*sigma1^2 + (1-p)*sigma0^2 -- a different, and here very different,
# number.
sim_scale_target_dgp <- function(n, seed, g0 = 0, g1 = 1.2, b0 = 1, b1 = 0.5,
                                 sigma0 = 1, sigma1 = 2, tau = 0.3) {
  set.seed(seed)
  x1 <- runif(n, -1, 1)
  X <- cbind(1, x1)
  pi_x <- plogis(g0 + g1 * x1)
  D <- rbinom(n, 1, pi_x)
  eps <- ifelse(D == 1, rnorm(n, 0, sigma1), rnorm(n, 0, sigma0))
  deltaY <- b0 + b1 * x1 + tau * D + eps
  list(D = D, deltaY = deltaY, X = X, w = rep(1, n))
}

# Exact (numerically-integrated) population targets for the DGP above.
sim_scale_true_targets <- function(g0 = 0, g1 = 1.2, sigma0 = 1, sigma1 = 2) {
  pi_fun <- function(x) plogis(g0 + g1 * x)
  p <- integrate(function(x) pi_fun(x) * 0.5, -1, 1)$value
  r_fun <- function(x) (pi_fun(x) / (1 - pi_fun(x))) / (p / (1 - p))

  nu0s2 <- integrate(function(x) r_fun(x)^2 * (1 - pi_fun(x)) * 0.5, -1, 1)$value / (1 - p)
  nu2_whole_pop <- integrate(function(x) r_fun(x)^2 * 0.5, -1, 1)$value

  sigma0s2 <- sigma0^2
  sigma2_pooled_wrong <- p * sigma1^2 + (1 - p) * sigma0^2

  list(p = p, sigma0s2 = sigma0s2, sigma2_pooled_wrong = sigma2_pooled_wrong,
      nu0s2 = nu0s2, nu2_whole_pop_wrong = nu2_whole_pop)
}

# === A: sigma target ========================================================
test_that("sigma2.s converges to the control-only sigma_0s^2, not the pooled mix", {
  truth <- sim_scale_true_targets()
  sample <- sim_scale_target_dgp(n = 3000, seed = 42)

  nuis <- dml.sensemakr:::did_cell_scale_nuisances(sample, L = 5, seed = 1)
  comp <- dml.sensemakr:::did_cell_sigma2_nu2(sample, nuis)

  expect_equal(comp$sigma2.s, truth$sigma0s2, tolerance = 0.15)
  # must be much closer to the true control-only value than to the pooled one
  expect_true(abs(comp$sigma2.s - truth$sigma0s2) <
             0.3 * abs(truth$sigma2_pooled_wrong - truth$sigma0s2))
})

# === B: nu target ============================================================
test_that("nu2.s converges to E[(OX/O)^2|D=0], not the whole-population moment", {
  truth <- sim_scale_true_targets()
  sample <- sim_scale_target_dgp(n = 3000, seed = 42)

  nuis <- dml.sensemakr:::did_cell_scale_nuisances(sample, L = 5, seed = 1)
  comp <- dml.sensemakr:::did_cell_sigma2_nu2(sample, nuis)

  expect_equal(comp$nu2.s, truth$nu0s2, tolerance = 0.15)
  expect_true(abs(comp$nu2.s - truth$nu0s2) <
             0.5 * abs(truth$nu2_whole_pop_wrong - truth$nu0s2))
})

# === C: score (estimating-equation) tests ===================================
# sigma2.s/nu2.s are DEFINED by solving mean(a*beta_hat + b) = 0, so the
# resulting psi's should average to (numerically) zero regardless of the
# DGP -- a construction check, not a statistical one.
test_that("mean(psi.sigma2.s)/mean(psi.nu2.s) are numerically zero at the solved estimate", {
  sample <- sim_scale_target_dgp(n = 800, seed = 7)
  nuis <- dml.sensemakr:::did_cell_scale_nuisances(sample, L = 5, seed = 1)
  comp <- dml.sensemakr:::did_cell_sigma2_nu2(sample, nuis)

  expect_equal(mean(comp$psi.sigma2.s.cell), 0, tolerance = 1e-8)
  expect_equal(mean(comp$psi.nu2.s.cell), 0, tolerance = 1e-8)
})

# === D: external ATT preservation ===========================================
# (see also test-12/test-13's exact-reproduction tests -- this is a direct,
# self-contained check tied to this file's DGP/adapter path.)
test_that("theta.s/psi.theta.s are untouched: they equal fit$ATT/fit$att.inf.func exactly", {
  skip_if_not_installed("DRDID")
  sample <- sim_scale_target_dgp(n = 500, seed = 11)
  fit <- DRDID::drdid_panel(y1 = sample$deltaY, y0 = rep(0, length(sample$deltaY)),
                            D = sample$D, covariates = sample$X, i.weights = sample$w,
                            inffunc = TRUE)
  sr <- dml.sensemakr:::drdid_cell_short_results(fit, D = sample$D, deltaY = sample$deltaY,
                                                 X = sample$X, w = sample$w)
  expect_identical(sr$estimates$theta.s, fit$ATT)
  expect_identical(sr$psis$psi.theta.s, as.numeric(fit$att.inf.func))
})

# === Scope: unit weights only ==============================================
test_that("did_cell_scale_nuisances() errors clearly on non-unit weights", {
  sample <- sim_scale_target_dgp(n = 200, seed = 1)
  sample$w <- runif(200, 0.5, 1.5)
  expect_error(dml.sensemakr:::did_cell_scale_nuisances(sample, L = 5, seed = 1),
              "unit weights")
})

# === E & F: Monte Carlo SE calibration, including the joint (external ATT +
# scale) covariance in the bounds() theta.m/theta.p SEs ======================
# Repeats the whole did-adapter pipeline (DRDID::drdid_panel() for theta.s/
# psi.theta.s, the cross-fitted scale estimator for sigma2.s/nu2.s, then
# bounds()) across many simulated cells, and compares each analytical
# (plug-in influence-function) SE to the empirical SD across replications --
# including theta.m/theta.p, whose SEs bounds() builds from the JOINT
# psi.theta.s -/+ c*psi.S2 combination (see R/bias-bounds.R), not from
# psi.theta.s and psi.S2's variances added independently.
test_that("analytical SEs for sigma2.s/nu2.s/theta.m/theta.p are well-calibrated against Monte Carlo SD", {
  skip_if_not_installed("DRDID")

  n <- 600; M <- 300
  g0 <- 0; g1 <- 1.0; b0 <- 1; b1 <- 0.5
  sigma0 <- 1; sigma1 <- 1.3; tau <- 0.2

  one_rep <- function(seed) {
    set.seed(seed)
    x1 <- runif(n, -1, 1)
    X <- cbind(1, x1)
    pi_x <- plogis(g0 + g1 * x1)
    D <- rbinom(n, 1, pi_x)
    eps <- ifelse(D == 1, rnorm(n, 0, sigma1), rnorm(n, 0, sigma0))
    y1 <- b0 + b1 * x1 + tau * D + eps
    y0 <- rep(0, n)
    w <- rep(1, n)

    fit <- DRDID::drdid_panel(y1 = y1, y0 = y0, D = D, covariates = X,
                              i.weights = w, inffunc = TRUE)
    sr <- dml.sensemakr:::drdid_cell_short_results(fit, D = D, deltaY = y1 - y0,
                                                   X = X, w = w, cf.folds = 5, cf.seed = 1)
    b <- dml.sensemakr:::bounds(sr, cf.y = 0.03, cf.d = 0.04, rho2 = 1)
    c(sigma2.s = sr$estimates$sigma2.s,
      se.sigma2.s = dml.sensemakr:::psi.sd(sr$psis$psi.sigma2.s),
      nu2.s = sr$estimates$nu2.s,
      se.nu2.s = dml.sensemakr:::psi.sd(sr$psis$psi.nu2.s),
      theta.m = b$estimates$theta.m, se.theta.m = b$estimates$se.theta.m,
      theta.p = b$estimates$theta.p, se.theta.p = b$estimates$se.theta.p)
  }

  res <- t(vapply(1:M, one_rep, numeric(8)))

  calib <- function(estcol, secol) mean(res[, secol]) / sd(res[, estcol])

  ratio.sigma <- calib("sigma2.s", "se.sigma2.s")
  ratio.nu    <- calib("nu2.s", "se.nu2.s")
  ratio.thm   <- calib("theta.m", "se.theta.m")
  ratio.thp   <- calib("theta.p", "se.theta.p")

  # A ratio near 1 means the analytical (delta-method/plug-in) SE matches
  # the true sampling SD; generous bounds since this is a single L=5
  # cross-fit per replication, not a fully repeated-cross-fitting average.
  expect_true(ratio.sigma > 0.7 && ratio.sigma < 1.4)
  expect_true(ratio.nu    > 0.7 && ratio.nu    < 1.4)
  expect_true(ratio.thm   > 0.7 && ratio.thm   < 1.4)
  expect_true(ratio.thp   > 0.7 && ratio.thp   < 1.4)
})

# === Active propensity-score trimming: population-compatibility check =====
# DRDID::drdid_panel()'s default trim.level = 0.995 excludes any control
# unit whose fitted propensity is >= 0.995 from the external ATT (its own
# trim.ps[D==0] <- (ps.fit[D==0] < trim.level) rule -- confirmed against its
# source). did::att_gt(est_method = "dr") never exposes a trim.level
# argument, so it always inherits that same 0.995 default internally
# (confirmed against did's run_DRDID() source: it calls drdid_panel()
# without a trim.level argument). When trimming is actually active, theta.s
# (defined on the trimmed sample) and sigma2.s/nu2.s (defined on the full
# reconstructed cell -- did_cell_scale_nuisances() implements no trimming of
# its own) would target different analysis populations, so both adapters
# must detect and reject it (validate_no_active_trimming(),
# R/did-adapter.R) -- NOT merely reject trim.level < 1 on principle, since
# in most cells trim.level = 0.995 is set but never actually binds.
#
# Builds a small synthetic 2-period panel (id/time/group/x1/y) where one
# control ("never-treated") unit has an outlying x1 -- calibrated (via
# `extreme_x1`) to either leave every control's fitted propensity safely
# below 0.995 (no trim) or push exactly that one control's fitted
# propensity above it (active trim), while att_gt() itself still returns a
# non-NA att/se either way (i.e. this is NOT the same thing as did's own
# overlap-condition failure, which is a separate, coarser check).
build_trim_test_panel <- function(extreme_x1, seed = 11) {
  set.seed(seed)
  n_treat <- 30; n_control <- 300; n <- n_treat + n_control
  x1 <- c(rnorm(n_treat, 2, 1), rnorm(n_control - 1, 0, 1), extreme_x1)
  group <- c(rep(2L, n_treat), rep(0L, n_control))  # 0 = never-treated
  dat <- data.frame(id = rep(seq_len(n), 2), time = rep(c(1L, 2L), each = n),
                    group = rep(group, 2), x1 = rep(x1, 2))
  dat$y <- 1 + 0.1 * dat$x1 + ifelse(dat$group == 2 & dat$time == 2, 0.5, 0) +
    rnorm(nrow(dat))
  dat
}

fit_trim_test_mp <- function(dat) {
  suppressWarnings(did::att_gt(
    yname = "y", tname = "time", idname = "id", gname = "group", xformla = ~x1,
    data = dat, est_method = "dr", compute_inffunc = TRUE, bstrap = FALSE,
    print_details = FALSE, control_group = "nevertreated", base_period = "varying"
  ))
}

test_that("did_cell_short_results() succeeds when trim.level=0.995 but no control is actually trimmed", {
  skip_if_not_installed("did")
  mp <- fit_trim_test_mp(build_trim_test_panel(extreme_x1 = 1.5))
  expect_false(is.na(mp$att[mp$group == 2 & mp$t == 2]))  # sanity: not an overlap failure

  sr <- quietly(dml.sensemakr:::did_cell_short_results(mp, group = 2, t = 2))
  expect_true(is.finite(sr$estimates$sigma2.s))
  expect_true(is.finite(sr$estimates$nu2.s))
})

test_that("did_cell_short_results() errors when propensity-score trimming is actually active", {
  skip_if_not_installed("did")
  mp <- fit_trim_test_mp(build_trim_test_panel(extreme_x1 = 6))
  expect_false(is.na(mp$att[mp$group == 2 & mp$t == 2]))  # sanity: not an overlap failure

  expect_error(dml.sensemakr:::did_cell_short_results(mp, group = 2, t = 2),
              "TRIMMED")
})

test_that("drdid_cell_short_results() succeeds when trim.level=0.995 but no control is actually trimmed", {
  skip_if_not_installed("did"); skip_if_not_installed("DRDID")
  mp <- fit_trim_test_mp(build_trim_test_panel(extreme_x1 = 1.5))
  sample <- dml.sensemakr:::did_cell_sample(mp, group = 2, t = 2)
  fit <- DRDID::drdid_panel(y1 = sample$deltaY, y0 = rep(0, length(sample$deltaY)),
                            D = sample$D, covariates = sample$X, i.weights = sample$w,
                            inffunc = TRUE)

  sr <- quietly(dml.sensemakr:::drdid_cell_short_results(
    fit, D = sample$D, deltaY = sample$deltaY, X = sample$X, w = sample$w
  ))
  expect_true(is.finite(sr$estimates$sigma2.s))
  expect_true(is.finite(sr$estimates$nu2.s))
})

test_that("drdid_cell_short_results() errors when propensity-score trimming is actually active", {
  skip_if_not_installed("did"); skip_if_not_installed("DRDID")
  mp <- fit_trim_test_mp(build_trim_test_panel(extreme_x1 = 6))
  sample <- dml.sensemakr:::did_cell_sample(mp, group = 2, t = 2)
  fit <- DRDID::drdid_panel(y1 = sample$deltaY, y0 = rep(0, length(sample$deltaY)),
                            D = sample$D, covariates = sample$X, i.weights = sample$w,
                            inffunc = TRUE)

  expect_error(
    dml.sensemakr:::drdid_cell_short_results(
      fit, D = sample$D, deltaY = sample$deltaY, X = sample$X, w = sample$w
    ),
    "TRIMMED"
  )
})

test_that("drdid_cell_short_results() reads the fit's own trim.level, not a hard-coded default", {
  # A fit built with trim.level = 1 (i.e. no trimming, whatever the
  # propensity) must never be flagged, even if some control's fitted
  # propensity would have crossed the DEFAULT 0.995 threshold.
  skip_if_not_installed("did"); skip_if_not_installed("DRDID")
  mp <- fit_trim_test_mp(build_trim_test_panel(extreme_x1 = 6))
  sample <- dml.sensemakr:::did_cell_sample(mp, group = 2, t = 2)
  fit <- DRDID::drdid_panel(y1 = sample$deltaY, y0 = rep(0, length(sample$deltaY)),
                            D = sample$D, covariates = sample$X, i.weights = sample$w,
                            inffunc = TRUE, trim.level = 1)

  sr <- quietly(dml.sensemakr:::drdid_cell_short_results(
    fit, D = sample$D, deltaY = sample$deltaY, X = sample$X, w = sample$w
  ))
  expect_true(is.finite(sr$estimates$sigma2.s))
})
