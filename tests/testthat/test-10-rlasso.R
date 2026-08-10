# Test dml() with yreg = "rlasso" (hdm::rlasso), on real pension data.
# rlasso is registered as a caret-style Regression-only model (see
# R/model-list.R), so it can serve as the outcome model (yreg) but not the
# treatment model (dreg) for a binary treatment -- dreg is left as "ranger".
library(testthat)
library(dml.sensemakr)

skip_if_not_installed("hdm")

setup_rlasso <- local({
  data("pension", package = "dml.sensemakr")
  set.seed(11)
  idx <- sample(nrow(pension), 500)  # small subset for speed
  y <- pension$net_tfa[idx]
  d <- pension$e401[idx]
  x <- model.matrix(~ -1 + age + inc + educ + fsize + marr + twoearn + pira + hown,
                    data = pension[idx, ])
  fit <- quietly(dml(y, d, x, model = "plm", yreg = "rlasso", dreg = "ranger",
                     cf.folds = 2, cf.reps = 1, verbose = FALSE))
  list(fit = fit, y = y, d = d, x = x)
})

test_that("dml() with yreg = 'rlasso' fits successfully and records the method", {
  fit <- setup_rlasso$fit
  expect_s3_class(fit, "dml")
  expect_equal(attr(fit$info$yreg$yreg0$method, "name"), "rlasso")
})

test_that("dml() with yreg = 'rlasso' produces a finite ATE with a positive SE", {
  cf <- coef(setup_rlasso$fit)
  s <- se(setup_rlasso$fit)
  expect_true(all(is.finite(cf)))
  expect_true(all(s > 0))
})

test_that("dml() with yreg = 'rlasso' fits contain outcome predictions", {
  fit <- setup_rlasso$fit
  preds <- fit$fits[[1]]$preds
  expect_length(preds$yhat, length(setup_rlasso$y))
  expect_true(all(is.finite(preds$yhat)))
})
