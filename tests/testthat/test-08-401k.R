# Integration test: DML PLM with income-quartile groups (GATE), on real
# pension data, exercising the full reporting pipeline end to end
# (coef/se/confint, summary, confidence_bounds, dml_benchmark, dml_bounds,
# ovb_contour_plot, plot.dml.bounds). test-04-dml-plm.R covers the PLM
# machinery itself in isolation without groups; this file's job is to check
# that groups=, dml_benchmark(), dml_bounds(), and the various plot()
# methods all work together on one real fit.
library(testthat)
library(dml.sensemakr)

# dml_benchmark() re-evaluates the dml() call via eval(model.call), which
# happens in the package namespace -- so y/d/x/g must live in .GlobalEnv,
# not inside a local()'s private environment (see test-06-benchmarks.R).
data("pension", package = "dml.sensemakr")
set.seed(10)
idx <- sample(nrow(pension), 500)  # small subset for speed
pension_401k <- pension[idx, ]
assign("y", pension_401k$net_tfa, envir = .GlobalEnv)
assign("d", pension_401k$e401, envir = .GlobalEnv)
assign("x", model.matrix(~ -1 + age + inc + educ + fsize + marr + twoearn + pira + hown,
                         data = pension_401k), envir = .GlobalEnv)
assign("g", cut(x[, "inc"], quantile(x[, "inc"], c(0, 0.25, .5, .75, 1), na.rm = TRUE),
              labels = c("q1", "q2", "q3", "q4"), include.lowest = TRUE), envir = .GlobalEnv)
setup_401k <- list(
  fit = dml(y, d, x, groups = g, model = "plm", target = "ate", cf.folds = 2, cf.reps = 1, verbose = FALSE),
  y = y, d = d, x = x, g = g
)

test_that("dml() with groups= populates both main ATE and per-group GATEs", {
  fit <- setup_401k$fit
  expect_s3_class(fit, "dml")
  expect_true(!is.null(fit$results$main))
  expect_equal(names(fit$results$groups), levels(setup_401k$g))
})

test_that("coef/se/confint report a finite ATE with a valid confidence interval", {
  fit <- setup_401k$fit
  cf <- coef(fit)
  s <- se(fit)
  ci <- confint(fit)
  expect_true(all(is.finite(cf)))
  expect_true(all(s > 0))
  expect_true(all(ci[, 1] < ci[, 2]))

  cf_mean <- coef(fit, combine.method = "mean")
  expect_equal(length(cf_mean), length(cf))
})

test_that("summary()/print() report the model and don't error", {
  fit <- setup_401k$fit
  expect_output(print(fit), "Debiased Machine Learning")
  expect_output(print(summary(fit)), "Debiased Machine Learning")
})

test_that("confidence_bounds() returns valid sensitivity bounds for the main ATE", {
  cb <- confidence_bounds(setup_401k$fit, cf.y = 0.03, cf.d = 0.04)
  expect_s3_class(cb, "confidence.bounds")
  expect_true(all(cb[, "lwr"] < cb[, "upr"]))

  cb_mean <- confidence_bounds(setup_401k$fit, cf.y = 0.03, cf.d = 0.04, rho2 = 0,
                               combine.method = "mean", level = 0.975)
  expect_true(all(cb_mean[, "lwr"] < cb_mean[, "upr"]))
})

test_that("ovb_contour_plot() runs without error for the default and a custom bound", {
  expect_no_error(quietly(ovb_contour_plot(setup_401k$fit)))
  expect_no_error(quietly(ovb_contour_plot(setup_401k$fit, cf.y = 0.03, cf.d = 0.04, rho2 = 1,
                                           which.bound = "lwr",
                                           bound.label = "Max Match (3x years)",
                                           col.contour = "blue")))
})

test_that("dml_benchmark() produces a benchmark for a real covariate", {
  bench <- quietly(dml_benchmark(setup_401k$fit, benchmark_covariates = c("inc")))
  expect_s3_class(bench, "dml_benchmark")
  expect_equal(names(bench$benchmarks), "inc")
})

test_that("dml_bounds() produces ordered bias bounds around the point estimate", {
  bounds <- dml_bounds(model = setup_401k$fit, cf.y = 0.03, cf.d = 0.04, rho2 = 1)
  expect_s3_class(bounds, "dml.bounds")
  cf <- coef(bounds)
  expect_true(all(cf["theta.m", ] <= cf["theta.s", ]))
  expect_true(all(cf["theta.s", ] <= cf["theta.p", ]))
})

test_that("plot.dml.bounds(type = 'all') runs without error", {
  bounds <- dml_bounds(model = setup_401k$fit, cf.y = 0.03, cf.d = 0.04, rho2 = 1)
  expect_no_error(plot(bounds, type = "all"))
})
