# Test benchmark functions
#
# Note: dml_benchmark() re-evaluates the dml() call via eval(model.call).
# The eval happens in the package namespace, so variables from the stored call
# must be in the global environment.

library(testthat)
library(dml.sensemakr)

data("pension", package = "dml.sensemakr")
set.seed(77)
idx <- sample(nrow(pension), 500)
# Assign to global environment so eval(model.call) inside dml_benchmark can find them
assign("y", pension$net_tfa[idx], envir = .GlobalEnv)
assign("d", pension$e401[idx], envir = .GlobalEnv)
assign("x", model.matrix(~ -1 + age + inc + educ + fsize + marr + twoearn + pira + hown,
                           data = pension[idx, ]), envir = .GlobalEnv)
bench_fit <- dml(y, d, x, model = "plm", cf.folds = 2, cf.reps = 1, verbose = FALSE)

test_that("dml_benchmark returns correct class", {
  bench <- dml_benchmark(bench_fit, benchmark_covariates = c("inc"))
  expect_s3_class(bench, "dml_benchmark")
})

test_that("dml_benchmark with multiple covariates", {
  bench <- dml_benchmark(bench_fit, benchmark_covariates = c("inc", "pira"))
  expect_s3_class(bench, "dml_benchmark")
})

test_that("summary.dml_benchmark does not error", {
  bench <- dml_benchmark(bench_fit, benchmark_covariates = c("inc"))
  expect_output(print(summary(bench)), regexp = NULL)
})

# === benchmark_gain_y / benchmark_gain_d / benchmark_rho ===
# These are pure helpers (no model fitting needed), extracted from bench_fun()
# specifically so the rho=0 fallback branch (previously untestable without
# forcing a real, noisy model refit to happen to produce V.g<=0 or V.a<=0)
# can be tested directly and deterministically.
test_that("benchmark_gain_y computes the relative gain and floors at 0", {
  expect_equal(dml.sensemakr:::benchmark_gain_y(sigma.sq = 2, V.g = 1), 0.5)
  expect_equal(dml.sensemakr:::benchmark_gain_y(sigma.sq = 2, V.g = -1), 0)
})

test_that("benchmark_gain_d computes the relative gain and floors at 0", {
  expect_equal(dml.sensemakr:::benchmark_gain_d(nu.sq.wo = 4, V.a = 2), 0.5)
  expect_equal(dml.sensemakr:::benchmark_gain_d(nu.sq.wo = 4, V.a = -2), 0)
})

test_that("benchmark_rho computes the signed correlation when valid", {
  # |Bias|/sqrt(V.g*V.a) = |2|/sqrt(4*1) = 1
  expect_equal(dml.sensemakr:::benchmark_rho(Bias = 2, V.g = 4, V.a = 1, valid = TRUE), 1)
  expect_equal(dml.sensemakr:::benchmark_rho(Bias = -2, V.g = 4, V.a = 1, valid = TRUE), -1)
})

test_that("benchmark_rho caps the magnitude at 1", {
  # |Bias|/sqrt(V.g*V.a) = 10/sqrt(1*1) = 10, well above 1
  expect_equal(dml.sensemakr:::benchmark_rho(Bias = 10, V.g = 1, V.a = 1, valid = TRUE), 1)
  expect_equal(dml.sensemakr:::benchmark_rho(Bias = -10, V.g = 1, V.a = 1, valid = TRUE), -1)
})

test_that("benchmark_rho forces 0 when valid is FALSE, regardless of Bias", {
  # this is the branch that's hit whenever dropping a covariate doesn't
  # actually reduce outcome fit / treatment precision (V.g<=0 or V.a<=0) --
  # previously impossible to test without engineering a real model refit to
  # coincidentally land here
  expect_equal(dml.sensemakr:::benchmark_rho(Bias = 5, V.g = -1, V.a = 1, valid = FALSE), 0)
  expect_equal(dml.sensemakr:::benchmark_rho(Bias = 5, V.g = 1, V.a = -1, valid = FALSE), 0)
})

test_that("benchmark_rho handles a mixed-validity vector element-by-element", {
  Bias  <- c(2, 3, -4)
  V.g   <- c(4, -1, 1)
  V.a   <- c(1, 2, 4)
  valid <- V.g > 0 & V.a > 0
  result <- dml.sensemakr:::benchmark_rho(Bias, V.g, V.a, valid)
  expect_equal(result[1], 1)    # |2|/sqrt(4*1) = 1
  expect_equal(result[2], 0)    # forced to 0: V.g <= 0
  expect_equal(result[3], -1)   # |-4|/sqrt(1*4) = 2 -> capped at 1, signed negative
})

# === replication of Table 3 (Chernozhukov et al., 2026) ===
# Table 3 ("Decomposition of the observed confounding strength", Appendix C)
# reports, per covariate: Bias = theta_{s,j} - theta_{s,empty}, Alignment
# (rho), Trend (= sqrt(Gain.Y)) and Imbalance (= sqrt(Gain.D)), where
# "empty" denotes the estimate using NO covariates at all. That compares the
# full covariate set against the *empty* set, which is a different
# comparison than dml_benchmark() ever makes (it always compares the full
# set against "all covariates minus one"), so this cannot be replicated via
# an actual dml_benchmark() call. But benchmark_gain_y()/benchmark_gain_d()/
# benchmark_rho() are generic formulas that don't know or care which two
# nested covariate sets produced their inputs -- so they can still be
# checked directly against the manuscript's "Xall" row (dropping every
# covariate at once), for which the footnote below Table 3 gives the
# empty-covariate-set baseline exactly: sigma0s,empty^2 = 0.024 and
# nu0s,empty^2 = 1 ("because there are no covariates").
#
# Sign convention: benchmark_rho()'s formula, rho = Bias/sqrt(V.g*V.a) with
# NO leading minus sign, only reproduces the correct (true-correlation) sign
# of rho when Bias = theta.short.wo - theta.short (without minus with). The
# manuscript's own formula for this same quantity is theta_{s,j} -
# theta_{s,empty} = -rho*Trend*Imbalance*Scale -- opposite subtraction order
# AND an explicit leading minus sign, which nets out to the same rho. So the
# manuscript's raw reported Bias must be negated before it matches what
# benchmark_rho() expects as input (this is exactly what `bench_fun()` does
# internally: it computes Bias in this same without-minus-with direction for
# all of Cor/psi.rho, and only negates it when building the user-facing
# `delta` column).
test_that("benchmark_gain_y/benchmark_gain_d/benchmark_rho replicate Table 3's Xall row", {
  # empty-covariate-set baseline, given exactly in the note below Table 3
  sigma.sq.empty <- 0.024
  nu.sq.empty    <- 1

  # manuscript's reported values for the "Xall" row (dropping every
  # covariate at once); bias.manuscript is theta_{s,all} - theta_{s,empty}
  bias.manuscript       <- -0.0089
  alignment.manuscript  <- 0.194
  trend.manuscript      <- 0.117   # = sqrt(Gain.Y)
  imbalance.manuscript  <- 2.548   # = sqrt(Gain.D)

  # V.a follows directly, since the empty-set nu^2 is exactly 1
  V.a <- imbalance.manuscript^2 * nu.sq.empty

  # V.g and sigma.sq (using all covariates) solve simultaneously from
  # Gain.Y = V.g / sigma.sq  and  V.g = sigma.sq.empty - sigma.sq
  sigma.sq <- sigma.sq.empty / (1 + trend.manuscript^2)
  V.g <- sigma.sq.empty - sigma.sq

  expect_equal(dml.sensemakr:::benchmark_gain_y(sigma.sq, V.g),
              trend.manuscript^2, tolerance = 0.01)
  expect_equal(dml.sensemakr:::benchmark_gain_d(nu.sq.empty, V.a),
              imbalance.manuscript^2, tolerance = 0.01)

  bias.internal <- -bias.manuscript  # see sign note above
  expect_equal(dml.sensemakr:::benchmark_rho(bias.internal, V.g, V.a, valid = TRUE),
              alignment.manuscript, tolerance = 0.01)
})

# === additional coverage: input validation, structure, and correctness ===
# Shared fixtures (computed once, reused across the tests below) so we don't
# re-run dml_benchmark()'s expensive refits for every single assertion.
bench_single <- dml_benchmark(bench_fit, benchmark_covariates = "inc")
bench_multi  <- dml_benchmark(bench_fit, benchmark_covariates = c("inc", "pira"))

# a second base fit (npm, cf.reps = 3) covers two gaps at once: benchmarking
# on a nonparametric model, and summary.dml_benchmark's aggregation across
# cross-fitting repetitions -- neither is exercised by the tests above,
# which only use a plm model with cf.reps = 1.
npm_fit   <- dml(y, d, x, model = "npm", cf.folds = 2, cf.reps = 3, verbose = FALSE)
bench_npm <- dml_benchmark(npm_fit, benchmark_covariates = "inc")

test_that("dml_benchmark errors when a benchmark covariate is not found", {
  expect_error(dml_benchmark(bench_fit, benchmark_covariates = "not_a_covariate"),
              "not_a_covariate")
})

test_that("dml_benchmark$benchmarks is keyed by the requested covariates, in order", {
  expect_named(bench_multi$benchmarks, c("inc", "pira"))
  expect_named(bench_multi$benchmarks_psis, c("inc", "pira"))
})

test_that("dml_benchmark$benchmarks has the expected columns", {
  expect_equal(colnames(bench_single$benchmarks$inc),
              c("gain.Y", "gain.D", "rho", "theta.s", "theta.sj", "delta"))
})

test_that("dml_benchmark$benchmarks_psis has the expected structure", {
  psis <- bench_single$benchmarks_psis$inc
  expect_named(psis, c("psi.theta.s", "psi.sigma2.s", "psi.nu2.s",
                      "psi.theta.s.wo", "psi.sigma2.s.wo", "psi.nu2.s.wo",
                      "psi.GY", "psi.GD", "psi.rho"))
  # one influence-function vector per cross-fitting repetition (cf.reps = 1
  # for bench_fit), each with one value per observation
  expect_length(psis$psi.rho, 1)
  expect_length(psis$psi.rho[[1]], length(y))
})

test_that("delta equals the exact difference between the with/without short estimates", {
  # delta is reported in the manuscript's direction: theta.s (with the
  # covariate) minus theta.sj (without it) -- see the sign-convention note
  # above the Table 3 replication test.
  b <- bench_single$benchmarks$inc
  expect_equal(b$delta, b$theta.s - b$theta.sj)
})

test_that("gain.Y and gain.D are non-negative for every benchmark covariate", {
  for (covar in names(bench_multi$benchmarks)) {
    b <- bench_multi$benchmarks[[covar]]
    expect_true(all(b$gain.Y >= 0))
    expect_true(all(b$gain.D >= 0))
  }
})

test_that("rho is bounded in [-1, 1] for every benchmark covariate", {
  for (covar in names(bench_multi$benchmarks)) {
    expect_true(all(abs(bench_multi$benchmarks[[covar]]$rho) <= 1))
  }
})

test_that("benchmarks$theta.s matches the original model's own short estimate", {
  # guards against bench_fun() accidentally reading from the wrong model
  # (e.g. model.wo) when extracting the with-covariate estimate
  expect_equal(bench_single$benchmarks$inc$theta.s,
              dml.sensemakr:::extract_estimate(bench_fit$results$main[[1]], "theta.s"))
})

test_that("different benchmark covariates produce different results", {
  # guards against an indexing bug where the wrong column gets dropped
  # regardless of which covariate was requested
  inc  <- bench_multi$benchmarks$inc
  pira <- bench_multi$benchmarks$pira
  expect_false(isTRUE(all.equal(inc$delta, pira$delta)))
})

test_that("psi.GY and psi.GD have the same shape as psi.rho", {
  psis <- bench_single$benchmarks_psis$inc
  expect_length(psis$psi.GY, 1)
  expect_length(psis$psi.GY[[1]], length(y))
  expect_length(psis$psi.GD, 1)
  expect_length(psis$psi.GD[[1]], length(y))
})

test_that("all benchmark influence functions are finite", {
  for (covar in names(bench_multi$benchmarks_psis)) {
    psis <- bench_multi$benchmarks_psis[[covar]]
    expect_true(all(is.finite(unlist(psis$psi.GY))))
    expect_true(all(is.finite(unlist(psis$psi.GD))))
    expect_true(all(is.finite(unlist(psis$psi.rho))))
  }
})

test_that("print.dml_benchmark does not error", {
  expect_output(print(bench_single), regexp = NULL)
})

test_that("dml_benchmark works with a nonparametric (npm) model", {
  expect_s3_class(bench_npm, "dml_benchmark")
  expect_equal(colnames(bench_npm$benchmarks$inc),
              c("gain.Y", "gain.D", "rho", "theta.s", "theta.sj", "delta"))
})

test_that("summary.dml_benchmark returns a classed object with aggregated values", {
  s <- summary(bench_npm, combine.method = "mean")
  expect_s3_class(s, "summary_dml_benchmark")
  expect_equal(s$combine.method, "mean")
  # summary() must still carry the influence functions along, unlike before
  # the fix (where it returned a bare matrix and dropped everything else)
  expect_true(!is.null(s$benchmarks_psis))
})

test_that("summary.dml_benchmark aggregates across cross-fitting repetitions correctly", {
  raw <- bench_npm$benchmarks$inc  # one row per cf.rep (cf.reps = 3)

  s_mean   <- summary(bench_npm, combine.method = "mean")
  s_median <- summary(bench_npm, combine.method = "median")

  expect_equal(unname(s_mean$benchmarks["inc", "gain.Y"]), mean(raw$gain.Y))
  expect_equal(unname(s_mean$benchmarks["inc", "delta"]), mean(raw$delta))
  expect_equal(unname(s_median$benchmarks["inc", "gain.Y"]), median(raw$gain.Y))
})

test_that("summary.dml_benchmark defaults to median, matching summary.dml/summary.dml.bounds", {
  s_default <- summary(bench_npm)
  s_median  <- summary(bench_npm, combine.method = "median")
  expect_equal(s_default$combine.method, "median")
  expect_equal(s_default$benchmarks, s_median$benchmarks)
})

test_that("print.dml_benchmark defaults to median, matching print.dml", {
  expect_output(print(bench_npm), "median")
})

test_that("print.summary_dml_benchmark does not error and reports the combine method", {
  expect_output(print(summary(bench_single, combine.method = "mean")),
               "mean")
})

# Clean up global environment
rm(y, d, x, envir = .GlobalEnv)
