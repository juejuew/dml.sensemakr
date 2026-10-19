parametric_fixture <- function(n = 160L) {
  skip_if_not_installed("DRDID")
  set.seed(10615)
  X <- cbind(`(Intercept)` = 1, x = runif(n, -1, 1), z = runif(n, -1, 1))
  D <- rbinom(n, 1, plogis(0.15 + 0.6 * X[, "x"] - 0.4 * X[, "z"]))
  y0 <- rnorm(n)
  # The quadratic term also exercises nuisance corrections away from the
  # correctly specified outcome model; this test verifies estimator calculus.
  y1 <- y0 + 0.8 * D + 0.5 * X[, "x"] + 0.3 * X[, "z"] +
    0.5 * X[, "x"]^2 + rnorm(n)
  source <- DRDID::drdid_panel(y1, y0, D, X, inffunc = TRUE)
  input <- list(y1 = y1, y0 = y0, D = D, covariates = X)
  list(source = source, input = input,
       model = do.call(dml_from_drdid,
                        c(list(fit = source, method = "parametric"), input)))
}

test_that("parametric mode preserves the source and reconstructs full-sample nuisances", {
  f <- parametric_fixture()
  m <- f$model; r <- m$results$main$treat[[1]]
  X <- f$input$covariates; d <- f$input$D; y <- f$input$y1 - f$input$y0
  or <- lm.fit(X[d == 0, , drop = FALSE], y[d == 0])
  ps <- glm.fit(X, d, family = binomial(),
                 control = glm.control(epsilon = 1e-12, maxit = 100))
  expect_identical(r$estimates$theta.s, as.numeric(f$source$ATT))
  expect_identical(r$psis$psi.theta.s, as.numeric(f$source$att.inf.func))
  expect_equal(unname(se(m)), f$source$se, tolerance = 1e-10)
  expect_equal(m$fits[[1]]$preds$yhat0, as.numeric(X %*% or$coefficients), tolerance = 1e-8)
  expect_equal(m$fits[[1]]$preds$dhat, ps$fitted.values, tolerance = 1e-7)
  expect_equal(r$estimates$sigma2.s, mean(or$residuals^2), tolerance = 1e-10)
  p <- mean(d); e <- ps$fitted.values; h <- e / (1 - e) * (1 - p) / p
  expect_equal(r$estimates$nu2.s, 2 * mean(h[d == 1]) - mean(h[d == 0]^2), tolerance = 1e-7)
  expect_identical(m$info$cf.folds, 0L)
  expect_false(m$info$cross.fitting)
  expect_identical(m$adapter_info$mode, "parametric")
  expect_length(m$adapter_info$fold_id, 0)
  expect_length(m$adapter_info$rep_seeds, 0)
  expect_null(m$adapter_info$dml.att)
  expect_null(m$auxiliary_dml_results)
  expect_length(m$results$main$treat, 1)
  expect_output(print(m), "Full-sample|full-sample")
  set.seed(531)
  old.seed <- .Random.seed
  again <- do.call(dml_from_drdid,
                    c(list(fit = f$source, method = "parametric"), f$input))
  expect_identical(.Random.seed, old.seed)
  expect_identical(again$results, m$results)
})

# Independent oracle: perturb the empirical distribution, refit with base R
# glm.fit/lm.wfit, and differentiate the resulting scalar estimators. It does
# not call adapter moment/IF helpers, and is not a weighted adapter feature.
parametric_weight_oracle <- function(input, w) {
  X <- input$covariates; d <- input$D; y <- input$y1 - input$y0
  ps <- glm.fit(X, d, weights = w, family = quasibinomial(),
                 control = glm.control(epsilon = 1e-13, maxit = 100))
  # Quasibinomial has the identical estimating equation but permits the
  # fractional perturbation weights without binomial noninteger warnings.
  beta <- lm.wfit(X[d == 0, , drop = FALSE], y[d == 0], w[d == 0])$coefficients
  residual <- y - as.numeric(X %*% beta)
  p <- weighted.mean(d, w)
  odds <- ps$fitted.values / (1 - ps$fitted.values)
  theta <- weighted.mean(residual[d == 1], w[d == 1]) -
    weighted.mean(residual[d == 0], (w * odds)[d == 0])
  a <- weighted.mean(residual[d == 0]^2, w[d == 0])
  h <- odds * (1 - p) / p
  b <- 2 * weighted.mean(h[d == 1], w[d == 1]) -
    weighted.mean(h[d == 0]^2, w[d == 0])
  k <- sqrt(0.08 * 0.06 / (1 - 0.06))
  c(theta = theta, sigma2 = a, nu2 = b, S2 = a * b,
    lower = theta - k * sqrt(a * b))
}

test_that("full-sample joint IFs match empirical-weight numerical derivatives", {
  f <- parametric_fixture(96L)
  m <- f$model; r <- m$results$main$treat[[1]]; n <- length(f$input$D)
  w <- rep(1, n)
  step <- 1e-4
  numerical <- t(vapply(seq_len(n), function(i) {
    plus <- minus <- w
    plus[i] <- plus[i] + step; minus[i] <- minus[i] - step
    n * (parametric_weight_oracle(f$input, plus) -
           parametric_weight_oracle(f$input, minus)) / (2 * step)
  }, numeric(5)))
  b <- dml_bounds(m, cf.y = 0.08, cf.d = 0.06)
  analytic <- cbind(theta = r$psis$psi.theta.s,
    sigma2 = r$psis$psi.sigma2.s, nu2 = r$psis$psi.nu2.s,
    S2 = r$psis$psi.S2, lower = b$results$main$treat[[1]]$psis$psi.theta.m)
  cat("\nIF audit maximum absolute errors:\n")
  print(apply(abs(analytic - numerical), 2, max), digits = 5)
  expect_equal(unname(analytic), unname(numerical), tolerance = 2e-5)
  expect_lt(max(abs(colMeans(analytic))), 1e-7)
  expect_gt(max(abs(r$nuisance.correction$psi.nu2.s)), 1e-4)
})

test_that("parametric inference and downstream sensitivity use the same joint IF", {
  f <- parametric_fixture(); m <- f$model; r <- m$results$main$treat[[1]]
  e <- r$estimates
  phi <- cbind(r$psis$psi.theta.s, r$psis$psi.sigma2.s, r$psis$psi.nu2.s)
  n <- nrow(phi)
  V <- cov(phi) * (n - 1) / n^2
  expect_equal(unname(r$joint.covariance), unname(V), tolerance = 1e-12)
  k <- sqrt(0.08 * 0.06 / (1 - 0.06))
  grad <- c(1, -k / 2 * sqrt(e$nu2.s / e$sigma2.s),
            -k / 2 * sqrt(e$sigma2.s / e$nu2.s))
  expected <- sqrt(drop(t(grad) %*% V %*% grad))
  b <- dml_bounds(m, cf.y = 0.08, cf.d = 0.06)
  expect_equal(b$results$main$treat[[1]]$estimates$se.theta.m, expected, tolerance = 1e-9)
  zero <- confidence_bounds(m, cf.y = 0, cf.d = 0)
  expect_equal(unname(zero["att", ]),
    e$theta.s + c(-1, 1) * qnorm(0.95) * f$source$se, tolerance = 1e-9)
  expect_s3_class(sensemakr(m, cf.y = 0.08, cf.d = 0.06), "dml.sensemakr")
  expect_true(all(is.finite(robustness_value(m))))
  expect_s3_class(plot(m), "ggplot")
})

test_that("parametric benchmarks retain their method after serialization", {
  f <- parametric_fixture(); m <- unserialize(serialize(f$model, NULL))
  m$external_fit <- NULL; m$call <- quote(stop("not available"))
  X <- m$data$x[, c("(Intercept)", "x"), drop = FALSE]
  refit <- dml.sensemakr:::.drdid_benchmark_refit(m, X)
  source <- DRDID::drdid_panel(f$input$y1, f$input$y0, f$input$D, X, inffunc = TRUE)
  direct <- dml_from_drdid(source, y1 = f$input$y1, y0 = f$input$y0,
    D = f$input$D, covariates = X, method = "parametric")
  expect_identical(refit$adapter_info$mode, "parametric")
  expect_equal(refit$results, direct$results, tolerance = 1e-10)
  bench <- dml_benchmark(m, "z")
  expect_equal(bench$benchmarks$z$theta.sj, as.numeric(source$ATT))
  expect_equal(bench$benchmarks_psis$z$psi.nu2.s.wo[[1]],
               refit$results$main$treat[[1]]$psis$psi.nu2.s)
  expect_s3_class(benchmark_bounds(m, bench, kY = 0.5, kD = 0.5), "dml_benchmark_bounds")
  expect_error(dml_benchmark(m, "z", dreg = "glm"), "overrides")
})

test_that("intercept-only parametric fits have unit imbalance and zero imbalance IF", {
  f <- parametric_fixture()
  X <- f$input$covariates[, "(Intercept)", drop = FALSE]
  source <- DRDID::drdid_panel(f$input$y1, f$input$y0, f$input$D, X, inffunc = TRUE)
  m <- dml_from_drdid(source, y1 = f$input$y1, y0 = f$input$y0,
    D = f$input$D, covariates = X, method = "parametric")
  r <- m$results$main$treat[[1]]
  expect_equal(r$estimates$nu2.s, 1, tolerance = 1e-10)
  expect_equal(r$psis$psi.nu2.s, rep(0, length(f$input$D)), tolerance = 1e-8)
  expect_equal(unname(coef(m)), source$ATT)
  expect_equal(unname(se(m)), source$se, tolerance = 1e-10)
})

parametric_did_data <- function() {
  skip_if_not_installed("did")
  set.seed(10616)
  n <- 240L
  units <- data.frame(id = seq_len(n), g = rep(c(0, 3, 4), each = n / 3),
                      x = runif(n, -1, 1), z = runif(n, -1, 1), u = rnorm(n))
  do.call(rbind, lapply(1:4, function(t) {
    transform(units, time = t, w = 2,
      y = u + t * (0.2 + 0.4 * x + 0.3 * z) +
        0.8 * (g > 0 & t >= g) + rnorm(n))
  }))
}

test_that("did single-cell and batch paths preserve the parametric method", {
  dat <- parametric_did_data()
  for (fast in c(TRUE, FALSE)) {
    source <- did::att_gt(yname = "y", tname = "time", idname = "id", gname = "g",
      xformla = ~ x + z, data = dat, bstrap = FALSE, cband = FALSE,
      faster_mode = fast, weightsname = "w")
    m <- dml_from_did(source, 3, 4, method = "parametric")
    r <- m$results$main$treat[[1]]; cell <- m$adapter_info$cell
    expect_identical(r$estimates$theta.s, source$att[cell$column])
    expect_equal(r$psis$psi.theta.s,
      as.numeric(source$inffunc[cell$if_rows, cell$column]) * cell$if_scale,
      tolerance = 1e-12)
    expect_equal(unname(se(m)), source$se[cell$column], tolerance = 1e-10)
    expect_output(print(m), "Parametric did ATT")
    batch <- dml_from_did_cells(source,
      cells = data.frame(group = c(3, 3, 3), time = c(3, 2, 4)),
      method = "parametric")
    expect_identical(batch$status$status, c("ok", "error", "ok"))
    expect_equal(batch$models[[3]]$results, m$results)
    expect_identical(did_cell_apply(batch, confidence_bounds, cf.y = 0.05,
      cf.d = 0.05)$status$status, batch$status$status)
    restored <- unserialize(serialize(m, NULL))
    restored$external_fit <- NULL; restored$call <- quote(stop("unavailable"))
    X <- m$data$x[, c("(Intercept)", "x"), drop = FALSE]
    refit <- dml.sensemakr:::.did_benchmark_refit(restored, X)
    expect_identical(refit$adapter_info$mode, "parametric")
    expect_identical(refit$adapter_info$unit_id, m$adapter_info$unit_id)
    bench <- dml_benchmark(restored, "z")
    expect_equal(bench$benchmarks$z$theta.sj, unname(coef(refit)))
  }
})

test_that("parametric mode rejects unsupported settings instead of changing methods", {
  f <- parametric_fixture()
  run <- function(source = f$source, ...) do.call(dml_from_drdid,
    c(list(fit = source, method = "parametric"), f$input, list(...)))
  expect_error(run(dml_args = list(cf.folds = 1)), "empty dml_args")
  expect_error(run(dml_args = list(yreg = "lm")), "empty dml_args")
  expect_error(run(i.weights = seq_along(f$input$D)), "constant positive")
  imp <- do.call(DRDID::drdid_imp_panel, c(f$input, list(inffunc = TRUE)))
  expect_error(run(imp), "estMethod = 'trad'")
  trimmed <- do.call(DRDID::drdid_panel,
                      c(f$input, list(inffunc = TRUE, trim.level = 0.5)))
  expect_error(run(trimmed), "trimming")
  dat <- parametric_did_data(); dat$w <- 1 + dat$id / max(dat$id)
  source <- did::att_gt(yname = "y", tname = "time", idname = "id", gname = "g",
    xformla = ~ x + z, data = dat, bstrap = FALSE, cband = FALSE, weightsname = "w")
  expect_error(dml_from_did(source, 3, 4, method = "parametric"), "equal sampling weights")
  batch <- dml_from_did_cells(source, method = "parametric")
  expect_true(all(batch$status$status == "error"))
  expect_true(all(grepl("equal sampling weights", batch$status$message)))
  expect_error(dml_from_did(source, 3, 4, method = "parametric",
    dml_args = list(cf.seed = 1)), "empty dml_args")
})
