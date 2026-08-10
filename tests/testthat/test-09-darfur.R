# Integration test: DML NPM on the real darfur dataset (Cinelli & Hazlett
# 2020), exercising dml_benchmark()/robustness_value() on a real, messy
# applied covariate matrix. test-05-dml-npm.R covers the NPM machinery
# itself on synthetic pension data; this file's job is to check the same
# pipeline holds up on real applied data.
#
# village (a factor with hundreds of levels) is deliberately left out: its
# sparse per-level dummy columns make the ranger treatment-model fit
# unstable (nu^2 going negative and robustness_value() coming back NA on
# some fold splits even after median-combining several cf.reps) -- a real
# propensity-model limitation, not a bug, but not what this file is testing.
library(testthat)
library(dml.sensemakr)

# dml_benchmark() re-evaluates the dml() call via eval(model.call), which
# happens in the package namespace -- so darfur/x must live in .GlobalEnv,
# not inside a local()'s private environment (see test-06-benchmarks.R).
assign("darfur", { data("darfur", package = "sensemakr"); darfur }, envir = .GlobalEnv)
assign("x", model.matrix(~ age + farmer_dar + herder_dar + pastvoted + hhsize_darfur + female,
                         data = darfur), envir = .GlobalEnv)

set.seed(12)
setup_darfur <- list(
  fit = quietly(dml(as.numeric(darfur$peacefactor), darfur$directlyharmed, x,
                    model = "npm", cf.folds = 2, cf.reps = 2, verbose = FALSE))
)

test_that("dml() NPM fits the darfur data with a finite ATE", {
  fit <- setup_darfur$fit
  expect_s3_class(fit, "dml")
  expect_equal(fit$info$model, "npm")
  expect_true(is.finite(coef(fit)))
  expect_true(se(fit) > 0)
})

test_that("summary() reports the darfur ATE without error", {
  expect_output(print(summary(setup_darfur$fit)), "Debiased Machine Learning")
})

test_that("dml_benchmark() produces a benchmark for a real covariate", {
  bench <- quietly(dml_benchmark(setup_darfur$fit, benchmark_covariates = c("female")))
  expect_s3_class(bench, "dml_benchmark")
  expect_equal(names(bench$benchmarks), "female")
})

test_that("robustness_value() returns a value in [0, 1]", {
  rv <- robustness_value(setup_darfur$fit)
  expect_true(is.finite(rv))
  expect_true(rv >= 0 && rv <= 1)
})
