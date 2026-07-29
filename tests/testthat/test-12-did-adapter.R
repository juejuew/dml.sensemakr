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

# === Full (g,t) grid: did_to_dml() ===
# Every cell goes into results$groups (results$main/info$target left
# empty) -- verified against coef.dml()/se.dml()/confint.dml()/
# row_only_list()/dml_bounds()/robustness_value()/
# extreme_robustness_value() before this was written, and also uncovered a
# real, previously-untriggered bug in coef.dml()/se.dml() (sapply(NULL, f)
# silently coercing the whole coefficient vector into a list) fixed
# alongside this feature -- see R/print-summary-dml.R.
test_that("did_to_dml() builds a group per (group,t) cell, matching att_gt()'s own att exactly", {
  model <- dml.sensemakr:::did_to_dml(setup_did$mp)
  expect_equal(length(model$results$groups), length(setup_did$mp$group))
  expect_null(model$results$main)

  cf <- coef(model)
  expected <- setNames(setup_did$mp$att,
                       paste0("gate.g", setup_did$mp$group, "_t", setup_did$mp$t))
  expect_equal(unname(cf[names(expected)]), unname(expected))
})

test_that("did_to_dml() skips no cells on a grid with no failed att_gt() cells", {
  # The mpdta fixture's group=2004 cohort has exactly 20 treated units,
  # which legitimately triggers the small-treated-cohort warning added to
  # did_cell_nuisances() -- that's expected, informational, and not a
  # "skip"/failure. What this test actually checks is that every cell
  # still computes successfully (no "Skipping (group, t) = ..." warnings),
  # i.e. did_to_dml() ends up with one group per att_gt() cell, none lost.
  expect_false(any(is.na(setup_did$mp$att)))
  msgs <- c()
  model <- withCallingHandlers(
    dml.sensemakr:::did_to_dml(setup_did$mp),
    warning = function(w) { msgs <<- c(msgs, conditionMessage(w)); invokeRestart("muffleWarning") }
  )
  expect_false(any(grepl("^Skipping", msgs)))
  expect_equal(length(model$results$groups), length(setup_did$mp$group))
})

test_that("did_to_dml()'s confint has one row per cell and matches coef.dml() length", {
  model <- dml.sensemakr:::did_to_dml(setup_did$mp)
  ci <- confint(model)
  expect_equal(nrow(ci), length(coef(model)))
})

test_that("robustness_value()/extreme_robustness_value() run across the full grid, XRV <= RV everywhere", {
  model <- dml.sensemakr:::did_to_dml(setup_did$mp)
  rv  <- robustness_value(model, theta = 0, alpha = 0.05)
  xrv <- extreme_robustness_value(model, theta = 0, alpha = 0.05)
  expect_length(rv, length(model$results$groups))
  expect_length(xrv, length(model$results$groups))
  expect_true(all(rv >= 0 & rv <= 1))
  expect_true(all(xrv >= 0 & xrv <= 1))
  expect_true(all(xrv <= rv + 1e-6))
})

test_that("did_to_dml()'s dml_bounds(only=) for one group matches the single-cell adapter directly", {
  model <- dml.sensemakr:::did_to_dml(setup_did$mp)
  target_slot <- "g2004_t2006"

  full <- dml_bounds(model, cf.y = 0.03, cf.d = 0.04, rho2 = 1)
  restricted <- dml_bounds(model, cf.y = 0.03, cf.d = 0.04, rho2 = 1,
                           only = list(main = NULL, groups = target_slot))
  expect_equal(restricted$results$groups[[target_slot]], full$results$groups[[target_slot]])

  sr_direct <- dml.sensemakr:::did_cell_short_results(setup_did$mp, group = 2004, t = 2006)
  expect_equal(sr_direct$estimates$theta.s,
              model$results$groups[[target_slot]][[1]]$estimates$theta.s)
})

test_that("did_to_dml() warns clearly for control_group = 'notyettreated' (unverified path)", {
  data("mpdta", package = "did")
  mp_nyt <- att_gt(yname = "lemp", tname = "year", idname = "countyreal",
                   gname = "first.treat", xformla = ~lpop, data = mpdta,
                   est_method = "dr", compute_inffunc = TRUE, bstrap = FALSE,
                   print_details = FALSE, control_group = "notyettreated",
                   base_period = "varying")
  expect_warning(dml.sensemakr:::did_to_dml(mp_nyt), "not been independently verified")
})

# === Cross-checks for the three previously-unverified settings, following
# the same "reconstruct from DIDparams$data, call DRDID::drdid_panel()
# directly, compare to att_gt()'s own att/se" method that validated
# "nevertreated"/"varying" originally. ===

test_that("resolve_did_control_ids() for 'notyettreated' reproduces att_gt() exactly", {
  data("mpdta", package = "did")
  mp_nyt <- att_gt(yname = "lemp", tname = "year", idname = "countyreal",
                   gname = "first.treat", xformla = ~lpop, data = mpdta,
                   est_method = "dr", compute_inffunc = TRUE, bstrap = FALSE,
                   print_details = FALSE, control_group = "notyettreated",
                   base_period = "varying")
  sr <- dml.sensemakr:::did_cell_short_results(mp_nyt, group = 2004, t = 2004)
  k <- which(mp_nyt$group == 2004 & mp_nyt$t == 2004)
  expect_equal(sr$estimates$theta.s, mp_nyt$att[k], tolerance = 1e-8)
  expect_equal(sr$estimates$se.theta.s, mp_nyt$se[k], tolerance = 1e-8)
})

test_that("resolve_did_base_period() for 'universal', post-treatment cell, reproduces att_gt() exactly", {
  data("mpdta", package = "did")
  mp_uni <- att_gt(yname = "lemp", tname = "year", idname = "countyreal",
                   gname = "first.treat", xformla = ~lpop, data = mpdta,
                   est_method = "dr", compute_inffunc = TRUE, bstrap = FALSE,
                   print_details = FALSE, control_group = "nevertreated",
                   base_period = "universal")
  sr <- dml.sensemakr:::did_cell_short_results(mp_uni, group = 2004, t = 2004)
  k <- which(mp_uni$group == 2004 & mp_uni$t == 2004)
  expect_equal(sr$estimates$theta.s, mp_uni$att[k], tolerance = 1e-8)
  expect_equal(sr$estimates$se.theta.s, mp_uni$se[k], tolerance = 1e-8)
})

test_that("did_to_dml() skips base_period = 'universal's trivial self-comparison cells (se = NA)", {
  data("mpdta", package = "did")
  mp_uni <- att_gt(yname = "lemp", tname = "year", idname = "countyreal",
                   gname = "first.treat", xformla = ~lpop, data = mpdta,
                   est_method = "dr", compute_inffunc = TRUE, bstrap = FALSE,
                   print_details = FALSE, control_group = "nevertreated",
                   base_period = "universal")
  n_na <- sum(is.na(mp_uni$se))
  expect_true(n_na > 0)  # this fixture is expected to have such cells

  expect_warning(model <- dml.sensemakr:::did_to_dml(mp_uni), "trivial self-comparison")
  expect_equal(length(model$results$groups), length(mp_uni$group) - n_na)

  rv <- robustness_value(model)
  xrv <- extreme_robustness_value(model)
  expect_true(all(is.finite(rv)) && all(is.finite(xrv)))
  expect_true(all(xrv <= rv + 1e-6))
})

test_that("resolve_did_base_period() with anticipation > 0 reproduces att_gt() exactly", {
  data("mpdta", package = "did")
  mp_ant <- att_gt(yname = "lemp", tname = "year", idname = "countyreal",
                   gname = "first.treat", xformla = ~lpop, data = mpdta,
                   est_method = "dr", compute_inffunc = TRUE, bstrap = FALSE,
                   print_details = FALSE, control_group = "nevertreated",
                   base_period = "varying", anticipation = 1)
  sr <- dml.sensemakr:::did_cell_short_results(mp_ant, group = 2006, t = 2007)
  k <- which(mp_ant$group == 2006 & mp_ant$t == 2007)
  expect_equal(sr$estimates$theta.s, mp_ant$att[k], tolerance = 1e-8)
  expect_equal(sr$estimates$se.theta.s, mp_ant$se[k], tolerance = 1e-8)
})

# === Trimming warning ===
# nu2.s's sandwich correction doesn't differentiate through the trim
# indicator; verified by simulation (see R/did-adapter.R) to leave a
# small (~5%) residual SE understatement when trimming is genuinely
# active. did_cell_nuisances() warns whenever that happens.
test_that("did_cell_nuisances() warns when trimming is active, and not otherwise", {
  set.seed(99)
  n <- 500
  X <- cbind(1, c(runif(n - 1, -2, 2), 8))  # one extreme covariate value
  D <- c(rbinom(n - 1, 1, plogis(X[1:(n-1), 2])), 0)  # force the extreme unit to be a control
  deltaY <- rnorm(n)
  w <- rep(1, n)
  sample_trim <- list(D = D, deltaY = deltaY, X = X, w = w)
  expect_warning(dml.sensemakr:::did_cell_nuisances(sample_trim), "trimmed")

  sample_notrim <- list(D = D[1:(n-1)], deltaY = deltaY[1:(n-1)],
                        X = X[1:(n-1), , drop = FALSE], w = w[1:(n-1)])
  expect_no_warning(dml.sensemakr:::did_cell_nuisances(sample_notrim))
})

# === Small treated-cohort warning ===
# nu2.s's sandwich correction includes a 1/c1^3 term that grows as the
# treated group shrinks; confirmed by simulation (see R/did-adapter.R) to
# produce a real ~12% analytical-vs-bootstrap SE gap on a real cell with
# exactly 20 treated units. did_cell_nuisances() warns at n_treated <= 20.
test_that("did_cell_nuisances() warns for a small treated cohort, and not for a larger one", {
  set.seed(5)
  n <- 400
  X <- cbind(1, rnorm(n))
  D_small <- c(rep(1, 15), rep(0, n - 15))  # 15 treated: below the threshold
  D_large <- c(rep(1, 100), rep(0, n - 100))  # 100 treated: well above it
  deltaY <- rnorm(n)
  w <- rep(1, n)

  expect_warning(dml.sensemakr:::did_cell_nuisances(list(D = D_small, deltaY = deltaY, X = X, w = w)),
                "treated unit")
  expect_no_warning(dml.sensemakr:::did_cell_nuisances(list(D = D_large, deltaY = deltaY, X = X, w = w)))
})

# === Categorical / multiple covariates ===
# mpdta's own xformla (~lpop) only exercises a single continuous
# covariate; model.matrix()'s dummy-encoding of a factor covariate had
# never been tested end-to-end.
test_that("did_cell_short_results() reproduces att_gt() exactly with a factor covariate", {
  data("mpdta", package = "did")
  set.seed(3)
  county_region <- setNames(sample(c("north","south","east"),
                                   length(unique(mpdta$countyreal)), replace = TRUE),
                            unique(mpdta$countyreal))
  mpdta$region <- factor(county_region[as.character(mpdta$countyreal)])

  mp_factor <- att_gt(yname = "lemp", tname = "year", idname = "countyreal",
                      gname = "first.treat", xformla = ~region, data = mpdta,
                      est_method = "dr", compute_inffunc = TRUE, bstrap = FALSE,
                      print_details = FALSE, control_group = "nevertreated",
                      base_period = "varying")
  sr <- dml.sensemakr:::did_cell_short_results(mp_factor, group = 2004, t = 2006)
  k <- which(mp_factor$group == 2004 & mp_factor$t == 2006)
  expect_equal(sr$estimates$theta.s, mp_factor$att[k], tolerance = 1e-8)
  expect_equal(sr$estimates$se.theta.s, mp_factor$se[k], tolerance = 1e-8)
})

test_that("did_to_dml() runs end-to-end with mixed continuous + factor covariates", {
  data("mpdta", package = "did")
  set.seed(3)
  county_region <- setNames(sample(c("north","south","east"),
                                   length(unique(mpdta$countyreal)), replace = TRUE),
                            unique(mpdta$countyreal))
  mpdta$region <- factor(county_region[as.character(mpdta$countyreal)])

  mp_multi <- att_gt(yname = "lemp", tname = "year", idname = "countyreal",
                     gname = "first.treat", xformla = ~lpop + region, data = mpdta,
                     est_method = "dr", compute_inffunc = TRUE, bstrap = FALSE,
                     print_details = FALSE, control_group = "nevertreated",
                     base_period = "varying")
  sr <- dml.sensemakr:::did_cell_short_results(mp_multi, group = 2004, t = 2006)
  k <- which(mp_multi$group == 2004 & mp_multi$t == 2006)
  expect_equal(sr$estimates$theta.s, mp_multi$att[k], tolerance = 1e-8)
  expect_equal(sr$estimates$se.theta.s, mp_multi$se[k], tolerance = 1e-8)

  model <- dml.sensemakr:::did_to_dml(mp_multi)
  expect_equal(length(model$results$groups), length(mp_multi$group))
  rv <- robustness_value(model)
  xrv <- extreme_robustness_value(model)
  expect_true(all(xrv <= rv + 1e-6))
})

# sensemakr()/print()/summary()/plot() on did_to_dml()'s groups-only output
# (results$main = NULL). See R/did-adapter.R's header for what is and isn't
# expected to work here.
test_that("sensemakr() and print() work on did_to_dml() output", {
  model <- dml.sensemakr:::did_to_dml(setup_did$mp)
  sens <- suppressWarnings(sensemakr(model, cf.y = 0.03, cf.d = 0.04))
  expect_s3_class(sens, "dml.sensemakr")
  expect_equal(nrow(sens$sensitivity_stats), length(setup_did$mp$group))
  expect_equal(nrow(sens$conf.bounds), length(setup_did$mp$group))
  expect_null(sens$bench.bounds)
  expect_output(print(sens), "Robustness Values")
})

test_that("summary() works on did_to_dml() output, skipping cross-fitting-specific reporting", {
  model <- dml.sensemakr:::did_to_dml(setup_did$mp)
  sens <- suppressWarnings(sensemakr(model, cf.y = 0.03, cf.d = 0.04))
  expect_output(summary(sens), "Group Average Treatment Effect")
  expect_no_error(suppressWarnings(summary(sens)))
  # no cross-fitting info exists for a did-adapter object -- these sections
  # must be silently skipped rather than erroring on missing object$fits.
  expect_false(any(grepl("Cross-Fitting|ML Method", capture.output(suppressWarnings(summary(sens))))))
})

test_that("plot() on did_to_dml() output errors clearly without group=TRUE, works with it", {
  model <- dml.sensemakr:::did_to_dml(setup_did$mp)
  sens <- suppressWarnings(sensemakr(model, cf.y = 0.03, cf.d = 0.04))
  expect_error(plot(sens), "group = TRUE")
  expect_no_error(plot(sens, group = TRUE, group.number = 1))
})

test_that("sensemakr(benchmark_covariates=) warns rather than silently skipping for did output", {
  model <- dml.sensemakr:::did_to_dml(setup_did$mp)
  expect_warning(
    sens <- sensemakr(model, cf.y = 0.03, cf.d = 0.04, benchmark_covariates = "lpop"),
    "benchmarking is not available"
  )
  expect_null(sens$bench.bounds)
})
