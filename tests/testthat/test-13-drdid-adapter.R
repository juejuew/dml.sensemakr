# Tests for the standalone DRDID::drdid_panel() -> short.results adapter
# (R/drdid-adapter.R). Unlike the did::att_gt() adapter (R/did-adapter.R),
# DRDID's own return object does not carry its data along with it, so the
# caller must supply the raw (D, deltaY, X, w) sample themselves -- these
# tests exercise both the happy path and the consistency check that guards
# against a mismatched sample.
skip_if_not_installed("did")
skip_if_not_installed("DRDID")

library(did)

# Fixture: reuse the exact mpdta cell (group = 2004, t = 2006,
# control_group = "nevertreated", base_period = "varying") already
# validated against att_gt() elsewhere (test-12), then build the standalone
# DRDID::drdid_panel() fit a user would pass to this adapter.
setup_drdid <- local({
  data("mpdta", package = "did")
  mp <- att_gt(yname = "lemp", tname = "year", idname = "countyreal",
              gname = "first.treat", xformla = ~lpop, data = mpdta,
              est_method = "dr", compute_inffunc = TRUE, bstrap = FALSE,
              print_details = FALSE, control_group = "nevertreated",
              base_period = "varying")
  sample <- dml.sensemakr:::did_cell_sample(mp, group = 2004, t = 2006)
  fit <- DRDID::drdid_panel(
    y1 = sample$deltaY, y0 = rep(0, length(sample$deltaY)),
    D = sample$D, covariates = sample$X, i.weights = sample$w,
    inffunc = TRUE
  )
  list(mp = mp, sample = sample, fit = fit)
})

test_that("validate_drdid_fit() accepts a valid panel fit with inffunc = TRUE", {
  expect_true(dml.sensemakr:::validate_drdid_fit(setup_drdid$fit))
})

test_that("validate_drdid_fit() rejects non-'drdid' objects", {
  expect_error(dml.sensemakr:::validate_drdid_fit(list()), "class 'drdid'")
})

test_that("validate_drdid_fit() rejects a fit built without inffunc = TRUE", {
  s <- setup_drdid$sample
  fit_noif <- DRDID::drdid_panel(y1 = s$deltaY, y0 = rep(0, length(s$deltaY)),
                                 D = s$D, covariates = s$X, i.weights = s$w,
                                 inffunc = FALSE)
  expect_error(dml.sensemakr:::validate_drdid_fit(fit_noif), "inffunc = TRUE")
})

test_that("validate_drdid_fit() rejects repeated-cross-section (argu$panel = FALSE) fits", {
  fit_rc <- setup_drdid$fit
  fit_rc$argu$panel <- FALSE
  expect_error(dml.sensemakr:::validate_drdid_fit(fit_rc), "panel data")
})

test_that("drdid_cell_short_results() reproduces att_gt()'s theta.s/se.theta.s exactly", {
  sr <- quietly(dml.sensemakr:::drdid_cell_short_results(
    setup_drdid$fit, D = setup_drdid$sample$D, deltaY = setup_drdid$sample$deltaY,
    X = setup_drdid$sample$X, w = setup_drdid$sample$w
  ))
  k <- which(setup_drdid$mp$group == 2004 & setup_drdid$mp$t == 2006)
  expect_equal(sr$estimates$theta.s, setup_drdid$mp$att[k], tolerance = 1e-8)
  expect_equal(sr$estimates$se.theta.s, setup_drdid$mp$se[k], tolerance = 1e-6)
})

test_that("drdid_cell_short_results()'s sigma2.s/nu2.s exactly match the did-adapter's own values for the same cell", {
  sr_drdid <- suppressWarnings(dml.sensemakr:::drdid_cell_short_results(
    setup_drdid$fit, D = setup_drdid$sample$D, deltaY = setup_drdid$sample$deltaY,
    X = setup_drdid$sample$X, w = setup_drdid$sample$w
  ))
  sr_did <- quietly(dml.sensemakr:::did_cell_short_results(setup_drdid$mp, group = 2004, t = 2006))
  expect_equal(sr_drdid$estimates$sigma2.s, sr_did$estimates$sigma2.s)
  expect_equal(sr_drdid$estimates$nu2.s, sr_did$estimates$nu2.s)
})

test_that("drdid_cell_short_results()'s psi.theta.s equals fit$att.inf.func at the cell's own (unrescaled) scale", {
  sr <- quietly(dml.sensemakr:::drdid_cell_short_results(
    setup_drdid$fit, D = setup_drdid$sample$D, deltaY = setup_drdid$sample$deltaY,
    X = setup_drdid$sample$X, w = setup_drdid$sample$w
  ))
  expect_equal(sr$psis$psi.theta.s, as.numeric(setup_drdid$fit$att.inf.func))
  expect_length(sr$psis$psi.theta.s, setup_drdid$sample$n1)

  # did-adapter's own psi.theta.s for the same cell is this same influence
  # function rescaled by (n/n1) and zero-padded to the full panel -- the two
  # are consistent representations, not independent quantities (see
  # R/drdid-adapter.R's header).
  sr_did <- quietly(dml.sensemakr:::did_cell_short_results(setup_drdid$mp, group = 2004, t = 2006))
  idx <- match(setup_drdid$sample$ids,
              setup_drdid$mp$DIDparams$time_invariant_data[[setup_drdid$mp$DIDparams$idname]])
  n <- setup_drdid$mp$n
  n1 <- setup_drdid$sample$n1
  expect_equal(sr_did$psis$psi.theta.s[idx], sr$psis$psi.theta.s * (n / n1))
})

test_that("drdid_cell_short_results() rejects a sample that doesn't reproduce fit$ATT when refit", {
  s <- setup_drdid$sample
  bad_D <- s$D
  bad_D[1] <- 1L - bad_D[1]
  expect_error(
    suppressWarnings(dml.sensemakr:::drdid_cell_short_results(
      setup_drdid$fit, D = bad_D, deltaY = s$deltaY, X = s$X, w = s$w
    )),
    "does not reproduce fit\\$ATT"
  )
})

test_that("drdid_cell_short_results() rejects a sample whose length doesn't match fit$att.inf.func", {
  s <- setup_drdid$sample
  expect_error(
    dml.sensemakr:::drdid_cell_short_results(
      setup_drdid$fit, D = s$D[-1], deltaY = s$deltaY[-1], X = s$X[-1, , drop = FALSE], w = s$w[-1]
    ),
    "doesn't match"
  )
})

test_that("drdid_to_dml() builds a working dml object from a single fit", {
  s <- setup_drdid$sample
  model <- suppressWarnings(dml.sensemakr:::drdid_to_dml(list(
    g2004_t2006 = list(fit = setup_drdid$fit, D = s$D, deltaY = s$deltaY, X = s$X, w = s$w)
  )))
  expect_s3_class(model, "dml")
  expect_equal(model$info$model, "drdid")
  expect_null(model$results$main)
  expect_equal(names(model$results$groups), "g2004_t2006")

  rv <- robustness_value(model)
  xrv <- extreme_robustness_value(model)
  expect_true(is.finite(rv) && is.finite(xrv))
  expect_true(xrv <= rv + 1e-6)

  sens <- sensemakr(model, cf.y = 0.03, cf.d = 0.04)
  expect_s3_class(sens, "dml.sensemakr")
  expect_output(print(sens), "Robustness Values")
})

test_that("drdid_to_dml() assembles multiple independent fits into separate groups", {
  data("mpdta", package = "did")
  mp2 <- setup_drdid$mp
  s1 <- setup_drdid$sample
  fit1 <- setup_drdid$fit

  s2 <- dml.sensemakr:::did_cell_sample(mp2, group = 2006, t = 2006)
  fit2 <- DRDID::drdid_panel(y1 = s2$deltaY, y0 = rep(0, length(s2$deltaY)),
                             D = s2$D, covariates = s2$X, i.weights = s2$w,
                             inffunc = TRUE)

  model <- suppressWarnings(dml.sensemakr:::drdid_to_dml(list(
    g2004_t2006 = list(fit = fit1, D = s1$D, deltaY = s1$deltaY, X = s1$X, w = s1$w),
    g2006_t2006 = list(fit = fit2, D = s2$D, deltaY = s2$deltaY, X = s2$X, w = s2$w)
  )))
  expect_equal(names(model$results$groups), c("g2004_t2006", "g2006_t2006"))

  k2 <- which(mp2$group == 2006 & mp2$t == 2006)
  expect_equal(model$results$groups$g2006_t2006[[1]]$estimates$theta.s, mp2$att[k2], tolerance = 1e-8)

  cf <- coef(model)
  expect_length(cf, 2)
  rv <- robustness_value(model)
  expect_length(rv, 2)
})

test_that("drdid_to_dml() validates its input shape", {
  expect_error(dml.sensemakr:::drdid_to_dml(list()), "non-empty named list")
  expect_error(dml.sensemakr:::drdid_to_dml(list(1, 2)), "non-empty named list")
  expect_error(dml.sensemakr:::drdid_to_dml(list(a = list(), "")), "non-empty named list")
})
