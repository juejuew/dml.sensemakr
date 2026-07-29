# Tests for the did::att_gt() -> short.results adapter (R/did-adapter.R).
#
# These require both `did` and `DRDID` to be installed; skipped otherwise
# since they're Suggests, not hard dependencies.
skip_if_not_installed("did")
skip_if_not_installed("DRDID")

library(did)

# Fixture: did's own mpdta example, matching the exact cell (group = 2004,
# t = 2006, control_group = "nevertreated", base_period = "varying") that
# was independently validated by hand against att_gt()'s own reconstructed
# sample earlier in this project (theta.s = -0.1404483, se = 0.03537815).
setup_did <- local({
  data("mpdta", package = "did")
  mp <- att_gt(yname = "lemp", tname = "year", idname = "countyreal",
              gname = "first.treat", xformla = ~lpop, data = mpdta,
              est_method = "dr", compute_inffunc = TRUE, bstrap = FALSE,
              print_details = FALSE, control_group = "nevertreated",
              base_period = "varying")
  list(mp = mp)
})

test_that("validate_did_scope() accepts a valid panel/dr/unclustered MP object", {
  expect_true(dml.sensemakr:::validate_did_scope(setup_did$mp))
})

test_that("validate_did_scope() rejects non-MP objects", {
  expect_error(dml.sensemakr:::validate_did_scope(list()), "class 'MP'")
})

test_that("validate_did_scope() rejects est_method other than 'dr'", {
  data("mpdta", package = "did")
  mp_reg <- att_gt(yname = "lemp", tname = "year", idname = "countyreal",
                   gname = "first.treat", xformla = ~lpop, data = mpdta,
                   est_method = "reg", compute_inffunc = TRUE, bstrap = FALSE,
                   print_details = FALSE)
  expect_error(dml.sensemakr:::validate_did_scope(mp_reg), "est_method")
})

test_that("did_cell_short_results() reproduces att_gt()'s own theta.s/se.theta.s exactly", {
  sr <- dml.sensemakr:::did_cell_short_results(setup_did$mp, group = 2004, t = 2006)
  expect_equal(sr$estimates$theta.s, -0.1404483, tolerance = 1e-6)
  expect_equal(sr$estimates$se.theta.s, 0.03537815, tolerance = 1e-6)
})

test_that("did_cell_short_results() produces positive sigma2.s/nu2.s and finite SEs", {
  sr <- dml.sensemakr:::did_cell_short_results(setup_did$mp, group = 2004, t = 2006)
  expect_true(sr$estimates$sigma2.s > 0)
  expect_true(sr$estimates$nu2.s > 0)
  expect_true(is.finite(sr$estimates$se.S2))
  expect_true(is.finite(sr$estimates$cov.theta.S2))
})

test_that("did_cell_short_results() psi vectors are all full-sample length, zero outside the cell", {
  sr <- dml.sensemakr:::did_cell_short_results(setup_did$mp, group = 2004, t = 2006)
  n <- setup_did$mp$n
  expect_length(sr$psis$psi.theta.s, n)
  expect_length(sr$psis$psi.sigma2.s, n)
  expect_length(sr$psis$psi.nu2.s, n)
  # units outside this cell's sample must be exactly zero in the new psis,
  # matching att.inf.func's own zero-padding convention
  zero_in_theta <- which(sr$psis$psi.theta.s == 0)
  expect_true(all(sr$psis$psi.sigma2.s[zero_in_theta] == 0))
  expect_true(all(sr$psis$psi.nu2.s[zero_in_theta] == 0))
})

test_that("did_cell_short_results() output feeds cleanly into prep_bounds()/bounds()", {
  sr <- dml.sensemakr:::did_cell_short_results(setup_did$mp, group = 2004, t = 2006)
  fixed <- dml.sensemakr:::prep_bounds(sr)
  expect_equal(fixed$theta.s, sr$estimates$theta.s)

  b <- dml.sensemakr:::bounds(sr, cf.y = 0.03, cf.d = 0.04, rho2 = 1)
  expect_true(b$estimates$theta.m < b$estimates$theta.s)
  expect_true(b$estimates$theta.s < b$estimates$theta.p)
  expect_equal(b$estimates$theta.s - b$estimates$theta.m,
              b$estimates$theta.p - b$estimates$theta.s)
  expect_true(all(c(b$estimates$se.theta.s, b$estimates$se.bias.bound,
                    b$estimates$se.theta.m, b$estimates$se.theta.p) > 0))
})
