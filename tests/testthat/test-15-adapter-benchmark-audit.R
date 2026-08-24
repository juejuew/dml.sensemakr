# Audit tests for did_dml_benchmark()/drdid_dml_benchmark() (R/adapter-
# benchmarks.R) against the six-component hybrid benchmarking invariants:
# same analysis population for full/reduced external fits, active-trimming
# checked on BOTH fits, correctly aligned observation-level IFs for both
# the external ATT and the paper-native scale estimators, and joint (not
# independent) delta-method combination. See R/adapter-benchmarks.R's file
# header and check_same_did_population()/check_no_term_leakage() for the
# invariants these tests exercise.
skip_if_not_installed("did")
skip_if_not_installed("DRDID")

library(did); library(DRDID)

# === A: DRDID continuous covariate benchmark identities =====================
test_that("drdid_dml_benchmark(): Bias/V.g/V.a identities hold exactly, benchmark IFs average ~0", {
  set.seed(11)
  n <- 600
  x1 <- rnorm(n); x2 <- rnorm(n)
  pi_x <- plogis(0.3 * x1 + 0.8 * x2)
  D <- rbinom(n, 1, pi_x)
  y1 <- 1 + 0.2 * x1 + 1.2 * x2 + 0.3 * D + rnorm(n)
  y0 <- rep(0, n)
  X <- cbind(1, x1, x2); colnames(X) <- c("(Intercept)", "x1", "x2")
  Xred <- X[, c("(Intercept)", "x1"), drop = FALSE]
  w <- rep(1, n)

  fit <- DRDID::drdid_panel(y1 = y1, y0 = y0, D = D, covariates = X, i.weights = w, inffunc = TRUE)
  fit.wo <- DRDID::drdid_panel(y1 = y1, y0 = y0, D = D, covariates = Xred, i.weights = w, inffunc = TRUE)

  # independently-computed full/reduced short.results, NOT via the benchmark
  # function -- these are the ground truth the benchmark output is checked against
  sr    <- drdid_cell_short_results(fit, D = D, deltaY = y1 - y0, X = X, w = w)
  sr.wo <- drdid_cell_short_results(fit.wo, D = D, deltaY = y1 - y0, X = Xred, w = w)

  bench <- suppressWarnings(drdid_dml_benchmark(fit, D = D, deltaY = y1 - y0, X = X, w = w,
                                                benchmark_covariates = "x2"))
  b <- bench$benchmarks[["x2"]]
  p <- bench$benchmarks_psis[["x2"]]

  # full ATT genuinely differs from reduced ATT (x2 is a real confounder)
  expect_false(isTRUE(all.equal(sr$estimates$theta.s, sr.wo$estimates$theta.s)))
  expect_equal(b$theta.s, sr$estimates$theta.s)
  expect_equal(b$theta.sj, sr.wo$estimates$theta.s)

  # exact identities (definitions, not estimates -- should hold to floating-point precision)
  expect_equal(b$delta, sr.wo$estimates$theta.s - sr$estimates$theta.s)
  V.g <- sr.wo$estimates$sigma2.s - sr$estimates$sigma2.s
  V.a <- sr$estimates$nu2.s - sr.wo$estimates$nu2.s
  expect_equal(b$gain.Y, V.g / sr$estimates$sigma2.s)
  expect_equal(b$gain.D, V.a / sr.wo$estimates$nu2.s)

  # benchmark IFs (psi.GY/psi.GD/psi.rho) are linear combinations of
  # mean-zero influence functions with fixed (non-random) coefficients, so
  # their own means should be ~0 -- not floored/truncated, genuinely close to 0
  expect_equal(mean(p$psi.GY[[1]]), 0, tolerance = 1e-6)
  expect_equal(mean(p$psi.GD[[1]]), 0, tolerance = 1e-6)
  if (b$rho != 0) expect_equal(mean(p$psi.rho[[1]]), 0, tolerance = 1e-6)
})

# === B: did ATT(g,t) benchmark -- genuine refit, cell-by-cell alignment =====
test_that("did_dml_benchmark(): genuine per-cell refit with correctly aligned full/reduced IFs", {
  data("mpdta", package = "did")
  assign("mpdta", mpdta, envir = .GlobalEnv)

  mp <- att_gt(yname = "lemp", tname = "year", idname = "countyreal", gname = "first.treat",
              xformla = ~lpop, data = mpdta, est_method = "dr", compute_inffunc = TRUE,
              bstrap = FALSE, print_details = FALSE, control_group = "nevertreated",
              base_period = "varying")

  # independently-built "reduced" fit, NOT via the benchmark function
  mp.wo <- att_gt(yname = "lemp", tname = "year", idname = "countyreal", gname = "first.treat",
                  xformla = ~1, data = mpdta, est_method = "dr", compute_inffunc = TRUE,
                  bstrap = FALSE, print_details = FALSE, control_group = "nevertreated",
                  base_period = "varying")

  bench <- quietly(did_dml_benchmark(mp, benchmark_covariates = "lpop", cf.folds = 5, cf.seed = 1))
  expect_equal(names(bench), names(quietly(did_to_dml(mp, cf.folds = 5, cf.seed = 1))$results$groups))

  slot <- "g2004_t2006"
  k <- which(mp$group == 2004 & mp$t == 2006)
  k.wo <- which(mp.wo$group == 2004 & mp.wo$t == 2006)

  # genuine refit: reduced ATT for this cell equals the independently-built
  # xformla=~1 fit's own att_gt() output exactly, and differs from the full ATT
  expect_equal(bench[[slot]]$benchmarks[["lpop"]]$theta.sj, mp.wo$att[k.wo])
  expect_false(isTRUE(all.equal(bench[[slot]]$benchmarks[["lpop"]]$theta.s,
                                bench[[slot]]$benchmarks[["lpop"]]$theta.sj)))

  # full/reduced fits share the SAME comparison, not just the same unit set:
  # identical unit ids/order, D, deltaY, weights, and base-period/time
  # definition (item 3's strengthened check -- xformla only changes X, so
  # every one of these must match exactly, not just approximately).
  samp.full <- did_cell_sample(mp, group = 2004, t = 2006)
  samp.wo   <- did_cell_sample(mp.wo, group = 2004, t = 2006)
  expect_true(setequal(samp.full$ids, samp.wo$ids))
  ord <- match(samp.full$ids, samp.wo$ids)
  expect_identical(samp.full$pret, samp.wo$pret)   # same base period
  expect_identical(samp.full$t, samp.wo$t)         # same evaluation period
  expect_identical(samp.full$D, samp.wo$D[ord])
  expect_equal(samp.full$deltaY, samp.wo$deltaY[ord])
  expect_equal(samp.full$w, samp.wo$w[ord])

  # full/reduced external ATT IFs are aligned correctly: the benchmark's
  # stored psi.theta.s/psi.theta.s.wo for this cell must equal directly-
  # computed did_cell_short_results() calls on mp/mp.wo respectively
  sr.direct    <- quietly(did_cell_short_results(mp, group = 2004, t = 2006, cf.folds = 5, cf.seed = 1))
  sr.wo.direct <- quietly(did_cell_short_results(mp.wo, group = 2004, t = 2006, cf.folds = 5, cf.seed = 1))
  psi <- bench[[slot]]$benchmarks_psis[["lpop"]]
  expect_equal(psi$psi.theta.s[[1]], sr.direct$psis$psi.theta.s)
  expect_equal(psi$psi.theta.s.wo[[1]], sr.wo.direct$psis$psi.theta.s)

  # benchmark identities hold exactly for this cell
  b <- bench[[slot]]$benchmarks[["lpop"]]
  expect_equal(b$delta, b$theta.sj - b$theta.s)
  V.g <- sr.wo.direct$estimates$sigma2.s - sr.direct$estimates$sigma2.s
  V.a <- sr.direct$estimates$nu2.s - sr.wo.direct$estimates$nu2.s
  expect_equal(b$gain.Y, V.g / sr.direct$estimates$sigma2.s)
  expect_equal(b$gain.D, V.a / sr.wo.direct$estimates$nu2.s)
})

# === G: manual six-component delta-method reconstruction (drdid) ===========
# Independently rederives phi_Bias/phi_gainY/phi_gainD/phi_rho from the six
# observation-level IFs (phi_theta_full/reduced, phi_sigma2_full/reduced,
# phi_nu2_full/reduced -- all ALREADY-NORMALIZED influence functions, not
# raw scores: did_cell_sigma2_nu2() (R/did-adapter.R) divides its raw score
# by mean_a_sigma/mean_a_nu before returning psi.sigma2.s.cell/psi.nu2.s.cell,
# and DRDID::drdid_panel()'s att.inf.func is already the fully-constructed
# M-estimation IF) using the manuscript's own delta-method formulas, then
# compares element-by-element against what drdid_dml_benchmark() actually
# stores/reports -- not just approximately-zero IF means.
test_that("drdid_dml_benchmark(): manually reconstructed six-component delta-method IFs match to floating-point precision", {
  set.seed(11)
  n <- 600
  x1 <- rnorm(n); x2 <- rnorm(n)
  D <- rbinom(n, 1, plogis(0.3 * x1 + 0.8 * x2))
  y1 <- 1 + 0.2 * x1 + 1.2 * x2 + 0.3 * D + rnorm(n)
  y0 <- rep(0, n)
  X <- cbind(1, x1, x2); colnames(X) <- c("(Intercept)", "x1", "x2")
  Xred <- X[, c("(Intercept)", "x1"), drop = FALSE]
  w <- rep(1, n)

  fit    <- DRDID::drdid_panel(y1 = y1, y0 = y0, D = D, covariates = X, i.weights = w, inffunc = TRUE)
  fit.wo <- DRDID::drdid_panel(y1 = y1, y0 = y0, D = D, covariates = Xred, i.weights = w, inffunc = TRUE)
  sr    <- drdid_cell_short_results(fit, D = D, deltaY = y1 - y0, X = X, w = w)
  sr.wo <- drdid_cell_short_results(fit.wo, D = D, deltaY = y1 - y0, X = Xred, w = w)

  phi_theta_full <- sr$psis$psi.theta.s;       phi_theta_reduced <- sr.wo$psis$psi.theta.s
  phi_sigma2_full <- sr$psis$psi.sigma2.s;     phi_sigma2_reduced <- sr.wo$psis$psi.sigma2.s
  phi_nu2_full <- sr$psis$psi.nu2.s;           phi_nu2_reduced <- sr.wo$psis$psi.nu2.s
  sigma2.s <- sr$estimates$sigma2.s;           sigma2.s.wo <- sr.wo$estimates$sigma2.s
  nu2.s <- sr$estimates$nu2.s;                 nu2.s.wo <- sr.wo$estimates$nu2.s
  theta.s <- sr$estimates$theta.s;             theta.s.wo <- sr.wo$estimates$theta.s
  V.g <- sigma2.s.wo - sigma2.s; V.a <- nu2.s - nu2.s.wo
  Bias <- theta.s.wo - theta.s
  valid <- V.g > 0 && V.a > 0

  phi_Bias  <- phi_theta_reduced - phi_theta_full
  phi_gainY <- (1 / sigma2.s) * phi_sigma2_reduced - (sigma2.s.wo / sigma2.s^2) * phi_sigma2_full
  phi_gainD <- (1 / nu2.s.wo) * phi_nu2_full - (nu2.s / nu2.s.wo^2) * phi_nu2_reduced
  phi_rho <- if (!valid) rep(0, n) else {
    (1 / sqrt(V.g * V.a)) * (phi_theta_reduced - phi_theta_full) -
      (Bias / (2 * V.g^(3 / 2) * sqrt(V.a))) * (phi_sigma2_reduced - phi_sigma2_full) -
      (Bias / (2 * sqrt(V.g) * V.a^(3 / 2))) * (phi_nu2_full - phi_nu2_reduced)
  }

  bench <- suppressWarnings(drdid_dml_benchmark(fit, D = D, deltaY = y1 - y0, X = X, w = w,
                                                benchmark_covariates = "x2"))
  p <- bench$benchmarks_psis[["x2"]]; b <- bench$benchmarks[["x2"]]

  tol <- 1e-10
  expect_equal(phi_theta_full, p$psi.theta.s[[1]], tolerance = tol)
  expect_equal(phi_theta_reduced, p$psi.theta.s.wo[[1]], tolerance = tol)
  expect_equal(phi_sigma2_full, p$psi.sigma2.s[[1]], tolerance = tol)
  expect_equal(phi_sigma2_reduced, p$psi.sigma2.s.wo[[1]], tolerance = tol)
  expect_equal(phi_nu2_full, p$psi.nu2.s[[1]], tolerance = tol)
  expect_equal(phi_nu2_reduced, p$psi.nu2.s.wo[[1]], tolerance = tol)
  expect_equal(phi_gainY, p$psi.GY[[1]], tolerance = tol)
  expect_equal(phi_gainD, p$psi.GD[[1]], tolerance = tol)
  expect_equal(phi_rho, p$psi.rho[[1]], tolerance = tol)

  # the SE actually reported by summary() equals sqrt(mean(phi^2)/n) computed
  # from the manually reconstructed IFs -- for each of delta/gain.Y/gain.D/rho
  s <- summary(bench)
  expect_equal(unname(s$benchmarks["x2", "se.delta"]), psi.sd(phi_Bias), tolerance = tol)
  expect_equal(unname(s$benchmarks["x2", "se.gain.Y"]), psi.sd(phi_gainY), tolerance = tol)
  expect_equal(unname(s$benchmarks["x2", "se.gain.D"]), psi.sd(phi_gainD), tolerance = tol)
  expect_equal(unname(s$benchmarks["x2", "se.rho"]), psi.sd(phi_rho), tolerance = tol)

  # "final benchmark bound": plugging gain.Y/gain.D/|rho| into bounds() as a
  # (fixed, point-estimate) hypothetical sensitivity scenario -- the same
  # "what if the true confounder were this strong" convention bounds()/
  # confidence_bounds() already use everywhere else in this package. Its own
  # bf()-derived bias.bound IF, built from ONLY the full model's sigma2.s/
  # nu2.s/theta.s (gain.Y/gain.D/rho enter as fixed numbers, not re-propagated
  # -- consistent with cf.y/cf.d/rho2 always being treated as user-supplied
  # constants elsewhere in bounds()), reproduces bounds()'s own reported SE.
  if (b$gain.Y > 0 && b$gain.Y < 1 && b$gain.D > 0 && b$gain.D < 1) {
    bnd <- bounds(sr, cf.y = b$gain.Y, cf.d = b$gain.D, rho2 = b$rho^2)
    bf_manual <- sqrt(b$rho^2 * b$gain.Y * (b$gain.D / (1 - b$gain.D)))
    S_manual  <- sqrt(sigma2.s * nu2.s)
    psi.S2_manual <- sigma2.s * phi_nu2_full + nu2.s * phi_sigma2_full
    psi.bias.bound_manual <- (bf_manual / 2) * (1 / S_manual) * psi.S2_manual
    phi_theta_m <- phi_theta_full - psi.bias.bound_manual
    phi_theta_p <- phi_theta_full + psi.bias.bound_manual
    expect_equal(bnd$estimates$se.theta.m, psi.sd(phi_theta_m), tolerance = tol)
    expect_equal(bnd$estimates$se.theta.p, psi.sd(phi_theta_p), tolerance = tol)
  }
})

# === H: manual six-component delta-method reconstruction (did) =============
test_that("did_dml_benchmark(): manually reconstructed six-component delta-method IFs match to floating-point precision", {
  data("mpdta", package = "did")
  assign("mpdta", mpdta, envir = .GlobalEnv)
  mp <- att_gt(yname = "lemp", tname = "year", idname = "countyreal", gname = "first.treat",
              xformla = ~lpop, data = mpdta, est_method = "dr", compute_inffunc = TRUE,
              bstrap = FALSE, print_details = FALSE, control_group = "nevertreated",
              base_period = "varying")
  mp.wo <- att_gt(yname = "lemp", tname = "year", idname = "countyreal", gname = "first.treat",
                  xformla = ~1, data = mpdta, est_method = "dr", compute_inffunc = TRUE,
                  bstrap = FALSE, print_details = FALSE, control_group = "nevertreated",
                  base_period = "varying")

  sr    <- quietly(did_cell_short_results(mp, group = 2004, t = 2006, cf.folds = 5, cf.seed = 1))
  sr.wo <- quietly(did_cell_short_results(mp.wo, group = 2004, t = 2006, cf.folds = 5, cf.seed = 1))

  phi_theta_full <- sr$psis$psi.theta.s;       phi_theta_reduced <- sr.wo$psis$psi.theta.s
  phi_sigma2_full <- sr$psis$psi.sigma2.s;     phi_sigma2_reduced <- sr.wo$psis$psi.sigma2.s
  phi_nu2_full <- sr$psis$psi.nu2.s;           phi_nu2_reduced <- sr.wo$psis$psi.nu2.s
  sigma2.s <- sr$estimates$sigma2.s;           sigma2.s.wo <- sr.wo$estimates$sigma2.s
  nu2.s <- sr$estimates$nu2.s;                 nu2.s.wo <- sr.wo$estimates$nu2.s
  theta.s <- sr$estimates$theta.s;             theta.s.wo <- sr.wo$estimates$theta.s
  V.g <- sigma2.s.wo - sigma2.s; V.a <- nu2.s - nu2.s.wo
  Bias <- theta.s.wo - theta.s
  valid <- V.g > 0 && V.a > 0

  phi_Bias  <- phi_theta_reduced - phi_theta_full
  phi_gainY <- (1 / sigma2.s) * phi_sigma2_reduced - (sigma2.s.wo / sigma2.s^2) * phi_sigma2_full
  phi_gainD <- (1 / nu2.s.wo) * phi_nu2_full - (nu2.s / nu2.s.wo^2) * phi_nu2_reduced
  phi_rho <- if (!valid) rep(0, length(phi_theta_full)) else {
    (1 / sqrt(V.g * V.a)) * (phi_theta_reduced - phi_theta_full) -
      (Bias / (2 * V.g^(3 / 2) * sqrt(V.a))) * (phi_sigma2_reduced - phi_sigma2_full) -
      (Bias / (2 * sqrt(V.g) * V.a^(3 / 2))) * (phi_nu2_full - phi_nu2_reduced)
  }

  bench <- quietly(did_dml_benchmark(mp, benchmark_covariates = "lpop", cf.folds = 5, cf.seed = 1))
  p <- bench$g2004_t2006$benchmarks_psis[["lpop"]]

  tol <- 1e-10
  expect_equal(phi_theta_full, p$psi.theta.s[[1]], tolerance = tol)
  expect_equal(phi_theta_reduced, p$psi.theta.s.wo[[1]], tolerance = tol)
  expect_equal(phi_sigma2_full, p$psi.sigma2.s[[1]], tolerance = tol)
  expect_equal(phi_sigma2_reduced, p$psi.sigma2.s.wo[[1]], tolerance = tol)
  expect_equal(phi_nu2_full, p$psi.nu2.s[[1]], tolerance = tol)
  expect_equal(phi_nu2_reduced, p$psi.nu2.s.wo[[1]], tolerance = tol)
  expect_equal(phi_gainY, p$psi.GY[[1]], tolerance = tol)
  expect_equal(phi_gainD, p$psi.GD[[1]], tolerance = tol)
  expect_equal(phi_rho, p$psi.rho[[1]], tolerance = tol)

  s <- suppressWarnings(summary(bench$g2004_t2006))
  expect_equal(unname(s$benchmarks["lpop", "se.delta"]), psi.sd(phi_Bias), tolerance = tol)
  expect_equal(unname(s$benchmarks["lpop", "se.gain.Y"]), psi.sd(phi_gainY), tolerance = tol)
  expect_equal(unname(s$benchmarks["lpop", "se.gain.D"]), psi.sd(phi_gainD), tolerance = tol)
  expect_equal(unname(s$benchmarks["lpop", "se.rho"]), psi.sd(phi_rho), tolerance = tol)
})

# === C: Missing-value test -- must not silently let a unit re-enter ========
test_that("did_dml_benchmark() errors clearly when dropping a covariate changes the analysis population", {
  set.seed(1)
  n_treat <- 30; n_control <- 200; n <- n_treat + n_control
  x1 <- rnorm(n); x2 <- rnorm(n)
  na_unit_idx <- n_control + 5  # a treated unit
  x2[na_unit_idx] <- NA
  group <- c(rep(0L, n_control), rep(2L, n_treat))
  dat <- data.frame(id = rep(seq_len(n), 2), time = rep(c(1L, 2L), each = n),
                    group = rep(group, 2), x1 = rep(x1, 2), x2 = rep(x2, 2))
  dat$y <- 1 + 0.1 * dat$x1 + ifelse(is.na(dat$x2), 0, 0.1 * dat$x2) +
    ifelse(dat$group == 2 & dat$time == 2, 0.5, 0) + rnorm(nrow(dat))
  assign("dat", dat, envir = .GlobalEnv)

  mp <- suppressWarnings(att_gt(yname = "y", tname = "time", idname = "id", gname = "group",
                                xformla = ~x1 + x2, data = dat, est_method = "dr",
                                compute_inffunc = TRUE, bstrap = FALSE, print_details = FALSE,
                                control_group = "nevertreated", base_period = "varying"))

  # sanity: att_gt() itself really did drop the NA unit from the full fit
  expect_lt(mp$n, n)

  expect_error(
    suppressWarnings(did_dml_benchmark(mp, benchmark_covariates = "x2", cf.folds = 5, cf.seed = 1)),
    "analysis population"
  )
})

test_that("drdid_dml_benchmark() errors clearly on missing covariate values", {
  set.seed(2)
  n <- 100
  D <- rbinom(n, 1, 0.4)
  X <- cbind(1, rnorm(n)); colnames(X) <- c("(Intercept)", "x1")
  deltaY <- rnorm(n)
  fit <- DRDID::drdid_panel(y1 = deltaY, y0 = rep(0, n), D = D, covariates = X,
                            i.weights = rep(1, n), inffunc = TRUE)
  X_na <- X
  X_na[3, "x1"] <- NA
  expect_error(
    drdid_dml_benchmark(fit, D = D, deltaY = deltaY, X = X_na, benchmark_covariates = "x1"),
    "missing values"
  )
})

# === D: Trimming asymmetry -- reduced-X trims even though full-X does not ==
test_that("drdid_dml_benchmark() errors when dropping a covariate newly triggers active trimming", {
  set.seed(1)
  n <- 500
  x1 <- rnorm(n); x2 <- 0.9 * x1 + rnorm(n, sd = 0.3)
  x1[n] <- 5; x2[n] <- 0.9 * 5 * 0.3   # anomalous control: high x1, only "typical" x2
  score <- 1.5 * x1 + 1.5 * x2
  D <- c(rbinom(n - 1, 1, plogis(score[1:(n - 1)])), 0)
  y1 <- 1 + 0.3 * x1 + 0.3 * x2 + 0.2 * D + rnorm(n)
  y0 <- rep(0, n)
  X <- cbind(1, x1, x2); colnames(X) <- c("(Intercept)", "x1", "x2")
  w <- rep(1, n)

  fit <- DRDID::drdid_panel(y1 = y1, y0 = y0, D = D, covariates = X, i.weights = w, inffunc = TRUE)

  # sanity: the FULL fit itself does not trim
  expect_no_error(quietly(drdid_cell_short_results(fit, D = D, deltaY = y1 - y0, X = X, w = w)))

  # dropping x2 must trigger the (now-active) trimming check and error --
  # not silently compute sigma2.s/nu2.s on a population theta.s no longer matches
  expect_error(
    drdid_dml_benchmark(fit, D = D, deltaY = y1 - y0, X = X, w = w, benchmark_covariates = "x2"),
    "TRIMMED"
  )
})

# === E: Factor benchmark (did and drdid) ====================================
test_that("did_dml_benchmark() drops a factor term's dummy columns together", {
  data("mpdta", package = "did")
  set.seed(3)
  county_region <- setNames(sample(c("north", "south", "east"),
                                   length(unique(mpdta$countyreal)), replace = TRUE),
                            unique(mpdta$countyreal))
  mpdta$region <- factor(county_region[as.character(mpdta$countyreal)])
  assign("mpdta", mpdta, envir = .GlobalEnv)

  mp <- att_gt(yname = "lemp", tname = "year", idname = "countyreal", gname = "first.treat",
              xformla = ~lpop + region, data = mpdta, est_method = "dr", compute_inffunc = TRUE,
              bstrap = FALSE, print_details = FALSE, control_group = "nevertreated",
              base_period = "varying")

  bench <- quietly(did_dml_benchmark(mp, benchmark_covariates = "region", cf.folds = 5, cf.seed = 1))
  expect_true(length(bench) > 0)
  b <- bench[[1]]$benchmarks[["region"]]
  expect_true(is.finite(b$theta.sj))
  expect_false(isTRUE(all.equal(b$theta.s, b$theta.sj)))
})

test_that("drdid_dml_benchmark() can drop a factor's dummy columns together via a named group", {
  set.seed(4)
  n <- 400
  x1 <- rnorm(n)
  region <- factor(sample(c("north", "south", "east"), n, replace = TRUE))
  D <- rbinom(n, 1, plogis(0.3 * x1 + ifelse(region == "south", 0.5, 0)))
  y1 <- 1 + 0.2 * x1 + ifelse(region == "south", 0.8, 0) + 0.3 * D + rnorm(n)
  y0 <- rep(0, n)
  X <- model.matrix(~x1 + region)
  w <- rep(1, n)

  fit <- DRDID::drdid_panel(y1 = y1, y0 = y0, D = D, covariates = X, i.weights = w, inffunc = TRUE)
  region.cols <- setdiff(colnames(X), c("(Intercept)", "x1"))
  expect_true(length(region.cols) >= 2)  # confirms the factor really expanded to >1 dummy column

  # grouped: all dummy columns dropped together (the "consistent with did" semantics)
  bench.grouped <- suppressWarnings(drdid_dml_benchmark(
    fit, D = D, deltaY = y1 - y0, X = X, w = w,
    benchmark_covariates = list(region = region.cols)
  ))
  expect_equal(names(bench.grouped$benchmarks), "region")

  # default (ungrouped): each dummy benchmarked on its own -- documented,
  # deliberately different default from did's automatic whole-term dropping
  bench.single <- suppressWarnings(drdid_dml_benchmark(
    fit, D = D, deltaY = y1 - y0, X = X, w = w,
    benchmark_covariates = region.cols[1]
  ))
  expect_equal(names(bench.single$benchmarks), region.cols[1])
})

# === F: Interaction/transformation -- must not silently under-remove ======
test_that("did_dml_benchmark() rejects benchmarking a main effect that still leaks through an interaction", {
  data("mpdta", package = "did")
  set.seed(3)
  mpdta$x2 <- rnorm(nrow(mpdta))
  assign("mpdta", mpdta, envir = .GlobalEnv)

  mp <- att_gt(yname = "lemp", tname = "year", idname = "countyreal", gname = "first.treat",
              xformla = ~lpop + x2 + lpop:x2, data = mpdta, est_method = "dr",
              compute_inffunc = TRUE, bstrap = FALSE, print_details = FALSE,
              control_group = "nevertreated", base_period = "varying")

  expect_error(
    did_dml_benchmark(mp, benchmark_covariates = "lpop", cf.folds = 5, cf.seed = 1),
    "still leaves"
  )
  expect_error(
    did_dml_benchmark(mp, benchmark_covariates = "x2", cf.folds = 5, cf.seed = 1),
    "still leaves"
  )

  # benchmarking the interaction term itself is well-defined and unambiguous
  bench <- quietly(did_dml_benchmark(mp, benchmark_covariates = "lpop:x2", cf.folds = 5, cf.seed = 1))
  expect_true(length(bench) > 0)
})
