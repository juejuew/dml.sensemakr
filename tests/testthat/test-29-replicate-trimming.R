trimming_fixture <- function(n = 160L, cutoff = 0.55) {
  skip_if_not_installed("DRDID")
  set.seed(31031)
  X <- cbind(`(Intercept)` = 1, x = runif(n, -1, 1), z = runif(n, -1, 1))
  d <- rbinom(n, 1, plogis(0.15 + 0.9 * X[, "x"]))
  y0 <- rnorm(n)
  y1 <- y0 + 0.8 * d + 0.5 * X[, "x"] + 0.3 * X[, "z"] + rnorm(n)
  input <- list(y1 = y1, y0 = y0, D = d, covariates = X)
  source <- do.call(DRDID::drdid_panel, c(input, list(inffunc = TRUE, trim.level = cutoff)))
  model <- do.call(dml_from_drdid, c(list(fit = source, method = "replicate"), input,
                                    list(unit_id = seq_len(n))))
  list(source = source, input = input, model = model, cutoff = cutoff)
}

# Independent base-R oracle; all nuisance models are refitted on every
# empirical-weight perturbation, and the control mask is recomputed each time.
trimming_weight_oracle <- function(input, w, cutoff) {
  X <- input$covariates; d <- input$D; y <- input$y1 - input$y0
  ps <- glm.fit(X, d, weights = w, family = quasibinomial(),
                 control = glm.control(epsilon = 1e-13, maxit = 100))
  beta <- lm.wfit(X[d == 0, , drop = FALSE], y[d == 0], w[d == 0])$coefficients
  r <- y - as.numeric(X %*% beta)
  e <- ps$fitted.values; odds <- e / (1 - e)
  keep <- d == 0 & e < cutoff
  theta <- weighted.mean(r[d == 1], w[d == 1]) -
    weighted.mean(r[keep], (w * odds)[keep])
  p <- weighted.mean(d, w); h <- odds * (1 - p) / p
  a <- weighted.mean(r[d == 0]^2, w[d == 0])
  b <- 2 * weighted.mean(h[d == 1], w[d == 1]) -
    weighted.mean(h[d == 0]^2, w[d == 0])
  k <- sqrt(0.08 * 0.06 / (1 - 0.06))
  c(theta = theta, sigma2 = a, nu2 = b, S2 = a * b, lower = theta - k * sqrt(a * b))
}

test_that("replicate preserves trimmed source estimates and full-sample sensitivity moments", {
  f <- trimming_fixture(); m <- f$model; r <- m$results$main$treat[[1]]
  e <- m$fits[[1]]$preds$dhat; d <- f$input$D; n <- length(d)
  tr <- m$adapter_info$diagnostics$source.trimming
  expect_true(tr$active)
  expect_gt(tr$n.trimmed.controls, 0)
  expect_identical(tr$trimmed.control.indices, which(d == 0 & e >= f$cutoff))
  expect_identical(tr$trimmed.control.ids, m$adapter_info$unit_id[tr$trimmed.control.indices])
  expect_equal(tr$n.retained.controls, sum(d == 0 & e < f$cutoff))
  expect_false(tr$boundary.conditions.verified)
  expect_length(m$fits[[1]]$preds$dhat, n)
  expect_equal(nrow(m$data$x), n)
  expect_identical(r$estimates$theta.s, as.numeric(f$source$ATT))
  expect_identical(r$psis$psi.theta.s, as.numeric(f$source$att.inf.func))
  expect_equal(m$parametric$reconstructed.att.if, as.numeric(f$source$att.inf.func), tolerance = 1e-10)
  expect_equal(unname(se(m)), f$source$se, tolerance = 1e-12)
  untrimmed <- do.call(DRDID::drdid_panel, c(f$input, list(inffunc = TRUE, trim.level = 1)))
  clean <- do.call(dml_from_drdid, c(list(fit = untrimmed, method = "replicate"), f$input))
  cr <- clean$results$main$treat[[1]]
  for (key in c("sigma2.s", "nu2.s", "S2")) {
    expect_equal(r$estimates[[key]], cr$estimates[[key]], tolerance = 1e-12)
    expect_equal(r$psis[[paste0("psi.", key)]], cr$psis[[paste0("psi.", key)]], tolerance = 1e-12)
  }
  expect_gt(abs(f$source$ATT - untrimmed$ATT), 1e-5)
  expect_gt(max(abs(m$parametric$moment.gradients$theta.beta)), 0.05)
  expect_output(print(m), "Source control trimming")
  expect_equal(unname(confidence_bounds(m, cf.y = 0, cf.d = 0)["att", ]),
    f$source$ATT + c(-1, 1) * qnorm(0.95) * f$source$se, tolerance = 1e-10)
  legacy <- do.call(dml_from_drdid, c(list(fit = f$source, method = "parametric"), f$input))
  expect_equal(legacy$results, m$results)
})

test_that("trimmed joint IFs and bounds match independent numerical derivatives", {
  f <- trimming_fixture(96L); m <- f$model; r <- m$results$main$treat[[1]]
  n <- length(f$input$D); w <- rep(1, n); step <- 1e-5
  expect_gt(m$adapter_info$diagnostics$source.trimming$min.control.distance, 1e-5)
  numerical <- t(vapply(seq_len(n), function(i) {
    plus <- minus <- w; plus[i] <- plus[i] + step; minus[i] <- minus[i] - step
    n * (trimming_weight_oracle(f$input, plus, f$cutoff) -
           trimming_weight_oracle(f$input, minus, f$cutoff)) / (2 * step)
  }, numeric(5)))
  b <- dml_bounds(m, cf.y = 0.08, cf.d = 0.06)
  analytic <- cbind(r$psis$psi.theta.s, r$psis$psi.sigma2.s, r$psis$psi.nu2.s,
                    r$psis$psi.S2, b$results$main$treat[[1]]$psis$psi.theta.m)
  cat("\nTrimmed IF numerical-derivative maximum error:", max(abs(analytic - numerical)), "\n")
  expect_equal(unname(analytic), unname(numerical), tolerance = 3e-5)
  V <- cov(analytic[, 1:3]) * (n - 1) / n^2
  expect_equal(unname(r$joint.covariance), unname(V), tolerance = 1e-12)
  k <- sqrt(0.08 * 0.06 / (1 - 0.06))
  grad <- c(1, -k/2 * sqrt(r$estimates$nu2.s/r$estimates$sigma2.s),
                 -k/2 * sqrt(r$estimates$sigma2.s/r$estimates$nu2.s))
  expect_equal(b$results$main$treat[[1]]$estimates$se.theta.m,
                sqrt(drop(t(grad) %*% V %*% grad)), tolerance = 1e-10)
})

test_that("trimming support does not bypass source validation or allow clipping and DML", {
  f <- trimming_fixture(); m <- f$model
  expect_error(do.call(dml_from_drdid, c(list(fit = f$source, method = "dml"), f$input)),
                 "Source propensity trimming")
  bad <- f$source; bad$att.inf.func <- rev(bad$att.inf.func)
  expect_error(do.call(dml_from_drdid, c(list(fit = bad, method = "replicate"), f$input)), "ordered IF")
  expect_error(dml.sensemakr:::.panel_parametric_aux(m$adapter_data, f$source$ATT,
    rev(as.numeric(f$source$att.inf.func)), m$refit_spec$source,
    m$adapter_info$diagnostics), "source ATT influence function")
  e <- m$fits[[1]]$preds$dhat; e[1] <- 1 - 1e-6
  expect_error(dml.sensemakr:::.panel_source_trimming(e, m$adapter_data,
    m$refit_spec$source, TRUE), "clipping")
  spec <- m$refit_spec$source; spec$trim.level <- min(e)/2
  e[1] <- m$fits[[1]]$preds$dhat[1]
  expect_error(dml.sensemakr:::.panel_source_trimming(e, m$adapter_data, spec, TRUE),
                 "no positive control residual weight")
})

test_that("trimmed DRDID benchmarks recompute source fits on the frozen full sample", {
  f <- trimming_fixture(); m <- unserialize(serialize(f$model, NULL))
  m$refit_spec$method <- "parametric" # legacy saved object
  X <- m$data$x[, c("(Intercept)", "z"), drop = FALSE]
  refit <- dml.sensemakr:::.drdid_benchmark_refit(m, X)
  source <- DRDID::drdid_panel(f$input$y1, f$input$y0, f$input$D, X,
                              inffunc = TRUE, trim.level = f$cutoff)
  direct <- dml_from_drdid(source, y1 = f$input$y1, y0 = f$input$y0,
    D = f$input$D, covariates = X, unit_id = m$adapter_info$unit_id, method = "replicate")
  expect_equal(refit$results, direct$results)
  expect_equal(refit$adapter_info$diagnostics, direct$adapter_info$diagnostics)
  expect_identical(refit$adapter_info$unit_id, m$adapter_info$unit_id)
  expect_equal(refit$refit_spec$source$trim.level, f$cutoff)
  expect_false(identical(refit$adapter_info$diagnostics$source.trimming$trimmed.control.indices,
                         m$adapter_info$diagnostics$source.trimming$trimmed.control.indices))
  b <- dml_benchmark(m, "x")
  expect_equal(b$benchmarks$x$theta.sj, as.numeric(source$ATT))
  expect_s3_class(benchmark_bounds(m, b, kY = 0.5, kD = 0.5), "dml_benchmark_bounds")
})

trimming_did_data <- function() {
  skip_if_not_installed("did")
  skip_if_not_installed("DRDID")
  set.seed(31032)
  n <- 6010L
  units <- data.frame(id = seq_len(n), x = c(rep(1, 4010), rep(0, 2000)),
    D = c(rep(1, 4000), rep(0, 10), rep(1, 1000), rep(0, 1000)),
    z = runif(n, -1, 1), u = rnorm(n))
  do.call(rbind, lapply(0:3, function(tt) transform(units, time = tt, g = 3 * D,
    y = u + tt * (0.4 * x + 0.3 * z + 0.2 * D) + 0.8 * D * (tt >= 3) + rnorm(n))))
}

test_that("did accepts active default trimming in both layouts, batches and benchmarks", {
  dat <- trimming_did_data()
  for (fast in c(TRUE, FALSE)) {
    source <- did::att_gt(yname = "y", tname = "time", idname = "id", gname = "g",
      xformla = ~ x + z, data = dat, bstrap = FALSE, cband = FALSE, faster_mode = fast)
    m <- dml_from_did(source, 3, 3, method = "replicate")
    r <- m$results$main$treat[[1]]; cell <- m$adapter_info$cell
    expect_true(m$adapter_info$diagnostics$trim.active)
    expect_equal(m$adapter_info$diagnostics$source.trimming$n.trimmed.controls, 10L)
    expect_equal(m$adapter_info$diagnostics$source.trimming$threshold, 0.995)
    expect_identical(r$estimates$theta.s, source$att[cell$column])
    expect_equal(r$psis$psi.theta.s, as.numeric(source$inffunc[cell$if_rows, cell$column]) * cell$if_scale,
                  tolerance = 1e-12)
    expect_equal(unname(se(m)), source$se[cell$column], tolerance = 1e-10)
    expect_output(print(m), "Source control trimming")
    expect_error(dml_from_did(source, 3, 3, method = "dml"), "Source propensity trimming")
    batch <- dml_from_did_cells(source, method = "replicate")
    expect_identical(batch$status$status, "ok")
    expect_equal(batch$models[[1]]$results, m$results)
    reduced <- dml.sensemakr:::.did_benchmark_refit(m, m$data$x[, c("(Intercept)", "z"), drop = FALSE])
    expect_false(reduced$adapter_info$diagnostics$trim.active)
    expect_equal(nrow(reduced$data$x), nrow(m$data$x))
    expect_identical(reduced$adapter_info$unit_id, m$adapter_info$unit_id)
    b <- dml_benchmark(m, "x")
    expect_equal(b$benchmarks$x$theta.sj, unname(coef(reduced)))
    for (type in c("bias", "scaled")) {
      pre <- did_pretrend_bounds(unserialize(serialize(m, NULL)), c(0, 1),
                                 k = c(0, 0.5, 1), type = type)
      expect_true(pre$pre_model$adapter_info$diagnostics$trim.active)
      expect_identical(pre$pre_model$adapter_info$unit_id, m$adapter_info$unit_id)
      expect_true(all(is.finite(pre$bounds$lwr)))
      pr <- pre$pre_model$results$main$treat[[1]]
      keys <- c("psi.theta.s", "psi.sigma2.s", "psi.nu2.s")
      joint <- cbind(do.call(cbind, r$psis[keys]), do.call(cbind, pr$psis[keys]))
      n <- nrow(joint)
      expect_equal(unname(pre$joint.covariance[[1]]), unname(cov(joint) * (n - 1) / n^2),
                    tolerance = 1e-10)
      expect_equal(pre$bounds$se.minus[1], unname(se(m)), tolerance = 1e-10)
    }
  }
})

test_that("pretrend covariance allows different trimming masks with identical units", {
  f <- trimming_fixture(); m <- f$model; X <- m$data$x
  X[, "x"] <- rev(X[, "x"])
  source <- DRDID::drdid_panel(f$input$y1 - 0.8*f$input$D, f$input$y0,
    f$input$D, X, inffunc = TRUE, trim.level = f$cutoff)
  pre <- dml_from_drdid(source, y1 = f$input$y1 - 0.8*f$input$D,
    y0 = f$input$y0, D = f$input$D, covariates = X,
    unit_id = m$adapter_info$unit_id, method = "replicate")
  expect_false(identical(m$adapter_info$diagnostics$source.trimming$trimmed.control.indices,
                         pre$adapter_info$diagnostics$source.trimming$trimmed.control.indices))
  b <- pretrend_bounds(m, pre, k = 1, type = "scaled")
  expect_length(b$unit_id, length(f$input$D))
  expect_true(all(is.finite(b$joint.covariance[[1]])))
  expect_true(is.finite(b$bounds$se.minus))
})
