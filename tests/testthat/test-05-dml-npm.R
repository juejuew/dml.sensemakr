# Test DML NPM fitting on a small dataset for speed.

library(testthat)
library(dml.sensemakr)

# Fit once and reuse across tests (small dataset, minimal settings)
setup_npm <- local({
  data("pension", package = "dml.sensemakr")
  set.seed(99)
  idx <- sample(nrow(pension), 500)  # small subset for speed
  y <- pension$net_tfa[idx]
  d <- pension$e401[idx]
  x <- model.matrix(~ -1 + age + inc + educ + fsize + marr + twoearn + pira + hown,
                     data = pension[idx, ])
  g <- cut(x[, "inc"], quantile(x[, "inc"], c(0, 0.5, 1), na.rm = TRUE),
           labels = c("low", "high"), include.lowest = TRUE)
  fit <- quietly(dml(y, d, x, model = "npm", groups = g, cf.folds = 2, cf.reps = 2, verbose = FALSE))
  list(fit = fit, y = y, d = d, x = x, groups = g)
})

# === dml() NPM structure ===
test_that("dml() NPM returns correct structure", {
  fit <- setup_npm$fit
  expect_s3_class(fit, "dml")
  expect_equal(fit$info$model, "npm")
})

test_that("dml() NPM fits have yhat0 and yhat1 predictions", {
  fit <- setup_npm$fit
  for (i in 1:2) {
    preds <- fit$fits[[i]]$preds
    expect_true(all(c("dhat", "yhat0", "yhat1", "phat") %in% names(preds)))
    expect_length(preds$dhat, length(setup_npm$y))
    expect_length(preds$yhat0, length(setup_npm$y))
    expect_length(preds$yhat1, length(setup_npm$y))
  }
})

# === Group results ===
test_that("dml() with groups produces group coefficients", {
  fit <- setup_npm$fit
  expect_true(!is.null(fit$coefs$groups))
  expect_true(length(fit$coefs$groups) > 0)
})

test_that("coef.dml returns both ATE and GATE", {
  cf <- coef(setup_npm$fit)
  expect_true(any(grepl("^ate", names(cf))))
  expect_true(any(grepl("^gate", names(cf))))
})

test_that("se.dml returns SEs for ATE and GATE", {
  s <- se(setup_npm$fit)
  expect_true(any(grepl("^ate", names(s))))
  expect_true(any(grepl("^gate", names(s))))
  expect_true(all(s > 0))
})

test_that("confint.dml with groups has correct number of rows", {
  ci <- confint(setup_npm$fit)
  cf <- coef(setup_npm$fit)
  expect_equal(nrow(ci), length(cf))
})

# === DML bounds with groups ===
test_that("dml_bounds with groups produces group bounds", {
  bounds <- quietly(dml_bounds(setup_npm$fit, cf.y = 0.04, cf.d = 0.03, rho2 = 1))
  expect_s3_class(bounds, "dml.bounds")
  expect_true(!is.null(bounds$coefs$groups))
})

test_that("coef.dml.bounds with groups returns matrix with GATE rows", {
  bounds <- quietly(dml_bounds(setup_npm$fit, cf.y = 0.04, cf.d = 0.03, rho2 = 1))
  cf <- coef(bounds)
  expect_true(is.matrix(cf))
  expect_true(any(grepl("^gate", colnames(cf))))
  expect_true(any(grepl("^ate", colnames(cf))))
})

# === confidence_bounds with groups ===
test_that("confidence_bounds with groups returns multi-row matrix", {
  bounds <- quietly(dml_bounds(setup_npm$fit, cf.y = 0.04, cf.d = 0.03, rho2 = 1))
  cb <- quietly(confidence_bounds(bounds))
  expect_true(nrow(cb) > 1)
  # Check that bounds exist (may have NaN with small samples, so just check structure)
  expect_equal(ncol(cb), 2)
  expect_equal(colnames(cb), c("lwr", "upr"))
})

# === sensitivity analysis with groups ===
test_that("sensemakr with groups works correctly", {
  sens <- quietly(sensemakr(setup_npm$fit, cf.y = 0.04, cf.d = 0.03))
  expect_s3_class(sens, "dml.sensemakr")
  expect_output(print(summary(sens)), "Sensitivity Analysis")
})

# === plotting ===
test_that("plot.dml does not error for NPM", {
  expect_silent(plot(setup_npm$fit))
})

test_that("plot.dml.bounds does not error", {
  bounds <- quietly(dml_bounds(setup_npm$fit, cf.y = 0.04, cf.d = 0.03, rho2 = 1))
  expect_silent(plot(bounds))
})

test_that("coef_plot returns a ggplot object", {
  # coef_plot takes numeric vectors, not a dml object
  cf <- coef(setup_npm$fit)
  ci <- confint(setup_npm$fit)
  p <- coef_plot(estimate = cf, lwr1 = ci[, 1], upr1 = ci[, 2],
                 labels = names(cf))
  expect_true(inherits(p, "gg") || inherits(p, "ggplot"))
})

test_that("ovb_contour_plot does not error", {
  expect_no_error(quietly(ovb_contour_plot(setup_npm$fit)))
})

# === closed-form XRV / caching equivalence: NPM, multi-target, and groups ===
# All robustness_value()/extreme_robustness_value() tests so far used a PLM
# fixture, which only ever exercises the trivial `main.slots <- rep("all", ...)`
# branch. A non-PLM model exercises the other branch,
# `unname(.target_to_slot[model$info$target])` -- the one structurally
# closest to the bug that was actually fixed (indexing something by a
# target-derived key) and, until now, never run by any test.
test_that("extreme_robustness_value's alpha = 1 closed form agrees with the general search (NPM)", {
  # Deliberately NOT setup_npm$fit: that fixture has groups, and group rows
  # never get the alpha = 1 closed-form shortcut (they fall back to the
  # general optim() search regardless of alpha -- see
  # extreme_robustness_value.dml()'s fallback condition). Combined with the
  # separate, pre-existing inefficiency that dml_bounds() recomputes every
  # slot AND every group on every single optim() evaluation (not something
  # this round's changes address), running this specific alpha = 1 vs.
  # alpha = 0.7 comparison against a grouped fit multiplies out to an
  # impractically large number of redundant bounds() calls -- and, with
  # ~125 observations per group per fold, a lot of accompanying
  # `nu2.s < 0` warnings. A small, single-target, group-free fit exercises
  # the same branch this test targets without that blowup.
  quietly({
    npm_fit <- dml(setup_npm$y, setup_npm$d, setup_npm$x, model = "npm",
                  cf.folds = 2, cf.reps = 1, verbose = FALSE)
    xrv_closed  <- extreme_robustness_value(npm_fit, theta = 0, alpha = 1)
    xrv_general <- extreme_robustness_value(npm_fit, theta = 0, alpha = 0.7)
  })
  expect_equal(unname(xrv_closed), unname(xrv_general), tolerance = 1e-3)
})

test_that("dml_bounds() gives identical output with and without a precomputed `fixed` (with groups)", {
  # setup_npm$fit was fit with groups = g, so model$results$groups is
  # populated -- this exercises the fixed$groups branch of dml_bounds(),
  # which no earlier test (all on ungrouped PLM fixtures) ever ran.
  expect_true(!is.null(setup_npm$fit$results$groups))
  quietly({
    fixed <- dml.sensemakr:::prep_dml_bounds(setup_npm$fit)
    without <- dml_bounds(setup_npm$fit, cf.y = 0.04, cf.d = 0.03, rho2 = 1)
    with    <- dml_bounds(setup_npm$fit, cf.y = 0.04, cf.d = 0.03, rho2 = 1, fixed = fixed)
  })
  expect_equal(with, without)
})

test_that("confidence_bounds.dml gives identical output with and without a precomputed `fixed` (with groups)", {
  quietly({
    fixed <- dml.sensemakr:::prep_dml_bounds(setup_npm$fit)
    without <- confidence_bounds(setup_npm$fit, cf.y = 0.04, cf.d = 0.03, rho2 = 1)
    with    <- confidence_bounds(setup_npm$fit, cf.y = 0.04, cf.d = 0.03, rho2 = 1, fixed = fixed)
  })
  expect_equal(with, without)
})

test_that("dml_bounds() gives identical output with and without a precomputed `fixed` (multi-target)", {
  # No existing fixture has more than one target; build one here reusing
  # setup_npm's data, since model$results$main having >1 slot is exactly
  # what exercises dml_bounds()'s Map()-based multi-slot threading.
  quietly({
    multi_fit <- dml(setup_npm$y, setup_npm$d, setup_npm$x, model = "npm",
                     target = c("ate", "att"), cf.folds = 2, cf.reps = 1, verbose = FALSE)
  })
  expect_true(length(multi_fit$results$main) >= 2)

  quietly({
    fixed <- dml.sensemakr:::prep_dml_bounds(multi_fit)
    without <- dml_bounds(multi_fit, cf.y = 0.04, cf.d = 0.03, rho2 = 1)
    with    <- dml_bounds(multi_fit, cf.y = 0.04, cf.d = 0.03, rho2 = 1, fixed = fixed)
  })
  expect_equal(with, without)
})

# === `only` restriction: dml_bounds() should compute a single slot/group ===
# without changing its result for that row, relative to the unrestricted
# (compute-everything) call. This is what lets robustness_value()/
# extreme_robustness_value()'s optim() search skip every OTHER slot/group on
# every evaluation -- the pre-existing inefficiency that caused a warning/
# runtime blowup on grouped fits (see row_only_list()).
test_that("row_only_list() row order matches confint()'s row order (main rows before group rows)", {
  fit <- setup_npm$fit
  conf <- confint(fit)
  row.only <- dml.sensemakr:::row_only_list(fit)
  expect_equal(length(row.only), nrow(conf))

  n.main <- length(fit$results$main)
  for (i in seq_len(n.main)) {
    expect_false(is.null(row.only[[i]]$main))
    expect_true(is.null(row.only[[i]]$groups))
  }
  group.names <- names(fit$results$groups)
  for (j in seq_along(group.names)) {
    i <- n.main + j
    expect_true(is.null(row.only[[i]]$main))
    expect_equal(row.only[[i]]$groups, group.names[j])
  }
})

test_that("dml_bounds() with `only` restricted to a single group matches the unrestricted result for that group", {
  group_name <- names(setup_npm$fit$results$groups)[1]
  quietly({
    full <- dml_bounds(setup_npm$fit, cf.y = 0.04, cf.d = 0.03, rho2 = 1)
    restricted <- dml_bounds(setup_npm$fit, cf.y = 0.04, cf.d = 0.03, rho2 = 1,
                             only = list(main = NULL, groups = group_name))
  })
  expect_null(restricted$results$main)
  expect_equal(restricted$results$groups[[group_name]], full$results$groups[[group_name]])
  expect_equal(restricted$coefs$groups[[group_name]], full$coefs$groups[[group_name]])
})

test_that("dml_bounds() with `only` restricted to the main slot matches the unrestricted result for that slot", {
  quietly({
    full <- dml_bounds(setup_npm$fit, cf.y = 0.04, cf.d = 0.03, rho2 = 1)
    restricted <- dml_bounds(setup_npm$fit, cf.y = 0.04, cf.d = 0.03, rho2 = 1,
                             only = list(main = "all", groups = NULL))
  })
  expect_null(restricted$results$groups)
  expect_equal(restricted$results$main[["all"]], full$results$main[["all"]])
})

test_that("dml_bounds() with `only` AND a precomputed `fixed` stays aligned to the right slot (groups)", {
  # Regression guard: dml_bounds() must subset fixed$main/fixed$groups by
  # NAME (not rely on Map()'s positional matching) so that restricting to
  # one group can't silently pair that group's data with a DIFFERENT
  # group's cached `fixed` quantities.
  group_name <- names(setup_npm$fit$results$groups)[1]
  quietly({
    fixed <- dml.sensemakr:::prep_dml_bounds(setup_npm$fit)
    without_fixed <- dml_bounds(setup_npm$fit, cf.y = 0.04, cf.d = 0.03, rho2 = 1,
                                only = list(main = NULL, groups = group_name))
    with_fixed <- dml_bounds(setup_npm$fit, cf.y = 0.04, cf.d = 0.03, rho2 = 1, fixed = fixed,
                             only = list(main = NULL, groups = group_name))
  })
  expect_equal(with_fixed$results$groups, without_fixed$results$groups)
  expect_equal(with_fixed$coefs$groups, without_fixed$coefs$groups)
})

test_that("confidence_bounds.dml with `only` matches the unrestricted computation, row by row (groups)", {
  fit <- setup_npm$fit
  quietly({
    full <- confidence_bounds(fit, cf.y = 0.04, cf.d = 0.03, rho2 = 1)
    row.only <- dml.sensemakr:::row_only_list(fit)
    expect_equal(length(row.only), nrow(full))
    for (i in seq_len(nrow(full))) {
      par <- rownames(full)[i]
      restricted <- confidence_bounds(fit, cf.y = 0.04, cf.d = 0.03, rho2 = 1, only = row.only[[i]])
      expect_equal(unname(restricted[par, ]), unname(full[par, ]), info = par)
    }
  })
})

test_that("confidence_bounds.dml with `only` matches the unrestricted computation, row by row (multi-target + groups)", {
  quietly({
    combo_fit <- dml(setup_npm$y, setup_npm$d, setup_npm$x, model = "npm",
                     target = c("ate", "att"), groups = setup_npm$groups,
                     cf.folds = 2, cf.reps = 1, verbose = FALSE)
  })
  expect_true(length(combo_fit$results$main) >= 2)
  expect_true(!is.null(combo_fit$results$groups))

  quietly({
    full <- confidence_bounds(combo_fit, cf.y = 0.04, cf.d = 0.03, rho2 = 1)
    row.only <- dml.sensemakr:::row_only_list(combo_fit)
    expect_equal(length(row.only), nrow(full))
    for (i in seq_len(nrow(full))) {
      par <- rownames(full)[i]
      restricted <- confidence_bounds(combo_fit, cf.y = 0.04, cf.d = 0.03, rho2 = 1, only = row.only[[i]])
      expect_equal(unname(restricted[par, ]), unname(full[par, ]), info = par)
    }
  })
})

test_that("robustness_value() and extreme_robustness_value() run correctly on a combined multi-target + groups fixture", {
  quietly({
    combo_fit <- dml(setup_npm$y, setup_npm$d, setup_npm$x, model = "npm",
                     target = c("ate", "att"), groups = setup_npm$groups,
                     cf.folds = 2, cf.reps = 1, verbose = FALSE)
    rv  <- robustness_value(combo_fit, theta = 0, alpha = 0.05)
    xrv <- extreme_robustness_value(combo_fit, theta = 0, alpha = 0.05)
  })
  expect_equal(names(rv), names(xrv))
  expect_true(all(is.na(rv) | (rv >= 0 & rv <= 1)))
  expect_true(all(is.na(xrv) | (xrv >= 0 & xrv <= 1)))
  # XRV <= RV always (same theta/alpha/rho2) -- see the algebraic argument
  # in the PLM-fixture test of the same fact.
  expect_true(all(is.na(rv) | is.na(xrv) | xrv <= rv + 1e-6))
})
