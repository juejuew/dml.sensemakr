# Test DML PLM fitting on a small dataset for speed
# Uses a small subset of pension data with minimal folds/reps.

library(testthat)
library(dml.sensemakr)

# Fit once and reuse across tests (small dataset, minimal settings)
setup_plm <- local({
  data("pension", package = "dml.sensemakr")
  set.seed(42)
  idx <- sample(nrow(pension), 500)  # small subset for speed
  y <- pension$net_tfa[idx]
  d <- pension$e401[idx]
  x <- model.matrix(~ -1 + age + inc + educ + fsize + marr + twoearn + pira + hown,
                     data = pension[idx, ])
  fit <- dml(y, d, x, model = "plm", cf.folds = 2, cf.reps = 2, verbose = FALSE)
  list(fit = fit, y = y, d = d, x = x)
})

# === dml() object structure ===
test_that("dml() returns correct class and structure for PLM", {
  fit <- setup_plm$fit
  expect_s3_class(fit, "dml")
  expect_true(all(c("data", "call", "info", "fits", "results", "coefs") %in% names(fit)))
  expect_equal(fit$info$model, "plm")
  expect_equal(fit$info$target, "ate")
  expect_equal(fit$info$cf.folds, 2)
  expect_equal(fit$info$cf.reps, 2)
})

test_that("dml() always stores a resolved cf.seed in its call, even when not supplied", {
  # cf.seed was not passed to dml() in setup_plm, so this checks that it
  # gets auto-resolved to a concrete value and stamped into fit$call --
  # this is what lets dml_benchmark()'s leave-one-out refit (built by
  # re-evaluating that same call with only `x` swapped) share the same
  # cross-fitting fold partition as the original fit, instead of drawing
  # an independent one.
  fit <- setup_plm$fit
  expect_true("cf.seed" %in% names(fit$call))
  expect_true(is.numeric(fit$call$cf.seed))
  expect_length(fit$call$cf.seed, 1)
})

test_that("dml() stores data correctly", {
  fit <- setup_plm$fit
  expect_equal(fit$data$y, setup_plm$y)
  expect_equal(fit$data$d, setup_plm$d)
  expect_equal(fit$data$x, setup_plm$x)
})

test_that("dml() PLM fits contain predictions", {
  fit <- setup_plm$fit
  expect_length(fit$fits, 2)  # cf.reps = 2
  for (i in 1:2) {
    preds <- fit$fits[[i]]$preds
    expect_true(all(c("dhat", "yhat") %in% names(preds)))
    expect_length(preds$dhat, length(setup_plm$y))
    expect_length(preds$yhat, length(setup_plm$y))
  }
})

test_that("dml() PLM produces main coefficients", {
  fit <- setup_plm$fit
  expect_true(!is.null(fit$coefs$main))
  # Main coefs are keyed by target (e.g., "all" for ATE)
  expect_true(length(fit$coefs$main) > 0)
  # Each element should be a 2x2 matrix (mean/median x estimate/se)
  first_coef <- fit$coefs$main[[1]]
  expect_true(is.matrix(first_coef))
  expect_equal(nrow(first_coef), 2)
  expect_equal(rownames(first_coef), c("mean", "median"))
  expect_equal(colnames(first_coef), c("estimate", "se"))
})

# === coef, se, confint ===
test_that("coef.dml returns named vector", {
  cf <- coef(setup_plm$fit)
  expect_type(cf, "double")
  expect_true(length(cf) > 0)
  # PLM with target="ate" produces names like "ate.all"
  expect_true(any(grepl("^ate", names(cf))))
})

test_that("coef.dml works with mean and median methods", {
  cf_median <- coef(setup_plm$fit, combine.method = "median")
  cf_mean <- coef(setup_plm$fit, combine.method = "mean")
  expect_type(cf_median, "double")
  expect_type(cf_mean, "double")
  expect_equal(length(cf_median), length(cf_mean))
})

test_that("se.dml returns positive standard errors", {
  s <- se(setup_plm$fit)
  expect_type(s, "double")
  expect_true(all(s > 0))
  expect_equal(names(s), names(coef(setup_plm$fit)))
})

test_that("confint.dml returns matrix with correct dimensions", {
  ci <- confint(setup_plm$fit, level = 0.95)
  expect_true(is.matrix(ci))
  expect_equal(ncol(ci), 2)
  # Lower bound should be less than upper bound
  expect_true(all(ci[, 1] < ci[, 2]))
})

test_that("confint.dml respects level parameter", {
  ci_95 <- confint(setup_plm$fit, level = 0.95)
  ci_99 <- confint(setup_plm$fit, level = 0.99)
  # 99% CI should be wider than 95%
  width_95 <- ci_95[1, 2] - ci_95[1, 1]
  width_99 <- ci_99[1, 2] - ci_99[1, 1]
  expect_true(width_99 > width_95)
})

# === summary ===
test_that("summary.dml returns correct class", {
  s <- summary(setup_plm$fit)
  expect_s3_class(s, "summary_dml")
})

test_that("print.dml does not error", {
  expect_output(print(setup_plm$fit), "Debiased Machine Learning")
})

test_that("print.summary_dml does not error", {
  expect_output(print(summary(setup_plm$fit)), "Debiased Machine Learning")
})

# === dml_bounds ===
test_that("dml_bounds returns correct class and structure", {
  bounds <- dml_bounds(setup_plm$fit, cf.y = 0.04, cf.d = 0.03, rho2 = 1)
  expect_s3_class(bounds, "dml.bounds")
  expect_equal(bounds$info$cf.y, 0.04)
  expect_equal(bounds$info$cf.d, 0.03)
  expect_equal(bounds$info$rho2, 1)
  expect_true(!is.null(bounds$coefs$main))
  expect_true(!is.null(bounds$dml.fit))
})

test_that("coef.dml.bounds returns matrix", {
  bounds <- dml_bounds(setup_plm$fit, cf.y = 0.04, cf.d = 0.03, rho2 = 1)
  cf <- coef(bounds)
  expect_true(is.matrix(cf))
  expect_true(nrow(cf) > 0)
  expect_true("theta.s" %in% rownames(cf))
  expect_true("bias.bound" %in% rownames(cf))
  expect_true("theta.m" %in% rownames(cf))
  expect_true("theta.p" %in% rownames(cf))
})

test_that("se.dml.bounds returns matrix with positive values", {
  bounds <- dml_bounds(setup_plm$fit, cf.y = 0.04, cf.d = 0.03, rho2 = 1)
  s <- se(bounds)
  expect_true(is.matrix(s))
  expect_true(all(s > 0))
})

test_that("confint.dml.bounds returns a list", {
  bounds <- dml_bounds(setup_plm$fit, cf.y = 0.04, cf.d = 0.03, rho2 = 1)
  ci <- confint(bounds)
  expect_type(ci, "list")
  expect_true(length(ci) > 0)
})

test_that("summary.dml.bounds does not error", {
  bounds <- dml_bounds(setup_plm$fit, cf.y = 0.04, cf.d = 0.03, rho2 = 1)
  expect_output(print(summary(bounds)), "Debiased Machine Learning")
})

# === confidence_bounds on dml object ===
test_that("confidence_bounds.dml returns correct structure", {
  cb <- confidence_bounds(setup_plm$fit, cf.y = 0.04, cf.d = 0.03, rho2 = 1)
  expect_s3_class(cb, "confidence.bounds")
  expect_true(is.matrix(cb))
  expect_equal(ncol(cb), 2)
  expect_equal(colnames(cb), c("lwr", "upr"))
  expect_true(all(cb[, "lwr"] < cb[, "upr"]))
})

# === confidence_bounds on dml.bounds object ===
test_that("confidence_bounds.dml.bounds returns correct structure", {
  bounds <- dml_bounds(setup_plm$fit, cf.y = 0.04, cf.d = 0.03, rho2 = 1)
  cb <- confidence_bounds(bounds)
  expect_s3_class(cb, "confidence.bounds")
  expect_true(is.matrix(cb))
  expect_equal(ncol(cb), 2)
  expect_true(all(cb[, "lwr"] < cb[, "upr"]))
})

# === dml_bounds()/prep_dml_bounds() caching equivalence, on a real fit ===
test_that("dml_bounds() gives identical output with and without a precomputed `fixed`", {
  fixed <- dml.sensemakr:::prep_dml_bounds(setup_plm$fit)
  without <- dml_bounds(setup_plm$fit, cf.y = 0.04, cf.d = 0.03, rho2 = 1)
  with    <- dml_bounds(setup_plm$fit, cf.y = 0.04, cf.d = 0.03, rho2 = 1, fixed = fixed)
  expect_equal(with, without)
})

test_that("confidence_bounds.dml gives identical output with and without a precomputed `fixed`", {
  fixed <- dml.sensemakr:::prep_dml_bounds(setup_plm$fit)
  without <- confidence_bounds(setup_plm$fit, cf.y = 0.04, cf.d = 0.03, rho2 = 1)
  with    <- confidence_bounds(setup_plm$fit, cf.y = 0.04, cf.d = 0.03, rho2 = 1, fixed = fixed)
  expect_equal(with, without)
})

# === robustness_value ===
test_that("robustness_value.dml returns named numeric", {
  rv <- robustness_value(setup_plm$fit)
  expect_type(rv, "double")
  expect_true(length(rv) > 0)
  expect_true(all(rv >= 0 & rv <= 1))
})

test_that("robustness_value.dml.bounds returns named numeric", {
  bounds <- dml_bounds(setup_plm$fit, cf.y = 0.04, cf.d = 0.03, rho2 = 1)
  rv <- robustness_value(bounds)
  expect_type(rv, "double")
  expect_true(length(rv) > 0)
  expect_true(all(rv >= 0 & rv <= 1))
})

# === extreme_robustness_value ===
test_that("extreme_robustness_value.dml returns named numeric", {
  xrv <- extreme_robustness_value(setup_plm$fit)
  expect_type(xrv, "double")
  expect_true(length(xrv) > 0)
  expect_true(all(xrv >= 0 & xrv <= 1))
})

test_that("extreme_robustness_value is never larger than robustness_value", {
  # XRV fixes cf.y = 1 (its maximum) and searches only cf.d, while RV
  # searches cf.y = cf.d jointly. Solving sqrt(cf.y*cf.d/(1-cf.d)) = k for a
  # fixed target k: RV solves r^2/(1-r) = k^2 (r = cf.y = cf.d), XRV solves
  # r/(1-r) = k^2 (r = cf.d, cf.y = 1 fixed). Since r^2 < r on (0,1), the
  # solution to the first equation is always >= the solution to the second
  # for the same k -- i.e. XRV <= RV always, for the same theta/alpha/rho2.
  rv  <- robustness_value(setup_plm$fit, theta = 0, alpha = 0.05)
  xrv <- extreme_robustness_value(setup_plm$fit, theta = 0, alpha = 0.05)
  expect_true(all(xrv <= rv + 1e-6))
})

test_that("extreme_robustness_value's alpha = 1 closed-form path runs without error", {
  xrv <- extreme_robustness_value(setup_plm$fit, theta = 0, alpha = 1)
  expect_type(xrv, "double")
  expect_true(all(is.finite(xrv)))
  expect_true(all(xrv >= 0 & xrv <= 1))
})

test_that("extreme_robustness_value's alpha = 1 closed form agrees with the general search", {
  # confidence_bounds.numeric() clamps level = 1 - alpha to a floor of 0.5
  # (level[level < 0.5] <- 0.5), so qnorm(level) = 0 for ANY alpha >= 0.5 --
  # the SE-adjustment term vanishes and the general optim()-based search
  # degenerates to exactly the same "point bound only" condition the closed
  # form solves. So alpha = 1 and, e.g., alpha = 0.7 should give the same
  # answer (up to optim()'s finite convergence tolerance), even though they
  # take completely different code paths internally.
  xrv_closed  <- extreme_robustness_value(setup_plm$fit, theta = 0, alpha = 1)
  xrv_general <- extreme_robustness_value(setup_plm$fit, theta = 0, alpha = 0.7)
  expect_equal(unname(xrv_closed), unname(xrv_general), tolerance = 1e-3)
})

# === NA short-circuit when nu2.s < 0 (row_is_undefined()) ===
# nu2.s does not depend on cf.y/cf.d, so if it's negative for even one
# cf.rep, S = sqrt(sigma2.s*nu2.s) is NaN for that rep regardless of
# cf.y/cf.d, and combine.mean()/combine.median() (no na.rm) propagate that
# into an NA combined bound for EVERY cf.y/cf.d -- the row's
# confidence_bounds() is a constant-NA function over the whole search
# domain. optim(method = "Brent") does not error or return NA on a
# constant-NA objective: it deterministically converges to its `upper`
# bound (verified empirically outside this test suite), which looks like an
# ordinary, meaningful RV/XRV but carries no real information. These tests
# confirm robustness_value()/extreme_robustness_value() detect this case up
# front (via row_is_undefined()) and return an honest NA instead.
test_that("row_is_undefined() detects a negative nu2.s in any cf.rep", {
  good_rep <- list(nu2.s = 0.5)
  bad_rep  <- list(nu2.s = -0.5)
  fixed <- list(main = list(all = list(good_rep, good_rep)),
               groups = list(low = list(good_rep, bad_rep)))
  expect_false(dml.sensemakr:::row_is_undefined(fixed, list(main = "all", groups = NULL)))
  expect_true(dml.sensemakr:::row_is_undefined(fixed, list(main = NULL, groups = "low")))
})

test_that("robustness_value()/extreme_robustness_value() return NA (not a misleading number) when nu2.s < 0", {
  # Hand-built PLM-shaped dml object (mirroring ate.plm()'s actual output
  # structure in R/short-parameters.R) with two cf.reps for a single "all"
  # slot: one well-behaved, one with a deliberately negative nu2.s. This
  # forces the pathological case directly and deterministically, rather
  # than relying on incidental sampling noise from a real ranger fit (which
  # is what actually triggers this in practice, but isn't reproducible
  # enough to hang a test on).
  n <- 200
  set.seed(321)
  x <- rnorm(n); d <- x + rnorm(n); y <- 2 * d + x + rnorm(n)
  make_rep <- function(nu2.s) {
    yhat <- x; dhat <- x
    resY <- y - yhat; resD <- d - dhat
    RRs <- resD / mean(resD^2)
    theta.s <- mean(resY * RRs)
    eresY <- resY - theta.s * resD
    psi.theta.s  <- eresY * RRs
    sigma2.s     <- mean(eresY^2)
    psi.sigma2.s <- eresY^2 - sigma2.s
    psi.nu2.s    <- rnorm(n, sd = 0.01)  # arbitrary; only nu2.s's sign matters here
    list(
      psis = list(psi.theta.s = psi.theta.s, psi.sigma2.s = psi.sigma2.s, psi.nu2.s = psi.nu2.s),
      estimates = list(theta.s = theta.s, se.theta.s = dml.sensemakr:::psi.sd(psi.theta.s),
                       sigma2.s = sigma2.s, nu2.s = nu2.s)
    )
  }
  good_rep <- make_rep(nu2.s = 1.2)
  bad_rep  <- make_rep(nu2.s = -0.8)

  model <- list(info = list(model = "plm", target = "ate"),
               results = list(main = list(all = list(good_rep, bad_rep))))
  model$coefs$main <- list(all = dml.sensemakr:::combine.cross.fits(model$results$main$all))
  class(model) <- "dml"

  expect_warning(rv <- robustness_value(model), "nu\\^2 is negative")
  expect_warning(xrv <- extreme_robustness_value(model), "nu\\^2 is negative")

  expect_true(is.na(rv[["ate.all"]]))
  expect_true(is.na(xrv[["ate.all"]]))
})

# === sensemakr ===
test_that("sensemakr.dml returns correct class", {
  sens <- sensemakr(setup_plm$fit, cf.y = 0.04, cf.d = 0.03)
  expect_s3_class(sens, "dml.sensemakr")
  expect_true(!is.null(sens$info))
  expect_true(!is.null(sens$model))
  expect_true(!is.null(sens$sensitivity_stats))
  expect_true(!is.null(sens$conf.bounds))
})

test_that("print.dml.sensemakr does not error", {
  sens <- sensemakr(setup_plm$fit, cf.y = 0.04, cf.d = 0.03)
  expect_output(print(sens), "Sensitivity Analysis")
})

test_that("summary.dml.sensemakr does not error", {
  sens <- sensemakr(setup_plm$fit, cf.y = 0.04, cf.d = 0.03)
  expect_output(print(summary(sens)), "Sensitivity Analysis")
})
