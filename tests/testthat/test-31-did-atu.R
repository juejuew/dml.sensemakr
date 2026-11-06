# ATU uses a different outcome arm and target population. Check independent
# cell extraction and hand-written scores, not merely argument forwarding.
atu_panel <- function() {
  skip_if_not_installed("did")
  set.seed(42012)
  n <- 360L
  u <- data.frame(id = sample(seq_len(n) * 11),
    g = rep(c(0, 3, 5), each = n / 3), x = runif(n, -1, 1),
    z = runif(n, -1, 1), u = rnorm(n))
  u$x <- u$x + 0.3 * (u$g == 3)
  do.call(rbind, lapply(0:5, function(tt) transform(u, time = tt,
    y = u + tt * (0.4 * x + 0.3 * z) +
      (1 + 0.5 * x) * (g > 0 & tt >= g) + rnorm(n))))
}

atu_source <- function(dat = atu_panel(), fast = TRUE, control = "nevertreated",
                       anticipation = 0, est = "dr", weightsname = NULL) {
  did::att_gt(yname = "y", tname = "time", idname = "id", gname = "g",
    xformla = ~ x + z, data = dat, panel = TRUE, bstrap = FALSE, cband = FALSE,
    compute_inffunc = FALSE, faster_mode = fast, est_method = est,
    control_group = control, anticipation = anticipation, weightsname = weightsname)
}

atu_args <- function(reps = 1, dirty = FALSE) {
  list(yreg = "lm", dreg = "glm", d.class = TRUE, cf.folds = 3,
    cf.reps = reps, cf.seed = 17, ps.trim = 0, dirty.tuning = dirty)
}

test_that("ATU refits the independently extracted control population and cell", {
  dat <- atu_panel()
  for (fast in c(TRUE, FALSE)) for (control in c("nevertreated", "notyettreated")) {
    fit <- atu_source(dat, fast, control, anticipation = 1)
    fit$att[] <- NA_real_  # source ATT/IF cannot supply an ATU
    m <- dml_from_did(fit, 3, 3, atu_args(2), method = "dml_refit", target = "atu")
    cohorts <- if (control == "nevertreated") c(0, 3) else c(0, 3, 5)
    ids <- sort(unique(dat$id[dat$g %in% cohorts]))
    before <- dat[dat$time == 1, ]; after <- dat[dat$time == 3, ]
    before <- before[match(ids, before$id), ]; after <- after[match(ids, after$id), ]
    X <- model.matrix(~ x + z, before); rownames(X) <- NULL; attr(X, "assign") <- NULL
    expect_identical(m$adapter_data$unit_id, ids)
    expect_equal(m$adapter_data$X, X)
    expect_equal(m$data$y, after$y - before$y)
    expect_equal(m$adapter_data$D, as.numeric(before$g == 3))
    expect_equal(m$adapter_info$target_population$n, sum(before$g != 3))
    expect_equal(m$adapter_info$cell$base_period, 1)
    expect_identical(m$info$target, "atu")
    expect_identical(names(m$results$main), "untr")
    expect_false(m$adapter_info$identification$verified)
    expect_identical(m$adapter_info$original.target, "att")
    native <- do.call(dml, c(list(y = after$y - before$y,
      d = as.numeric(before$g == 3), x = X), m$refit_spec$dml_args))
    expect_identical(m$coefs, native$coefs)
    for (r in seq_len(2)) {
      got <- m$results$main$untr[[r]]
      joint <- got$joint.covariance; got$joint.covariance <- NULL
      expect_identical(got, native$results$main$untr[[r]])
      expect_identical(m$fits[[r]]$preds, native$fits[[r]]$preds)
      expect_true(all(is.na(m$fits[[r]]$preds$yhat0)))
      psi <- cbind(theta.s = got$psis$psi.theta.s,
        sigma2.s = got$psis$psi.sigma2.s, nu2.s = got$psis$psi.nu2.s)
      expect_equal(colMeans(psi), setNames(rep(0, 3), colnames(psi)), tolerance = 1e-10)
      expect_equal(joint, crossprod(psi) / nrow(psi)^2, tolerance = 1e-12)
      expect_equal(unname(diag(joint)),
        unname(unlist(got$estimates[c("se.theta.s", "se.sigma2.s", "se.nu2.s")]))^2)
    }
    expect_null(m$auxiliary_dml_results)
    expect_null(m$refit_spec$source)
    expect_output(print(m), "DML refit did ATU")
    expect_output(print(m), "ordinary ATT parallel trends alone is insufficient")
  }
})

test_that("ATU effect, moments and IFs satisfy hand scores and role symmetry", {
  m <- dml_from_did(atu_source(), 3, 4, atu_args(2), method = "dml_refit", target = "atu")
  y <- m$data$y; d <- m$adapter_data$D
  for (k in seq_len(2)) {
    pred <- m$fits[[k]]$preds
    e <- pred$dhat; p <- pred$phat; g <- pred$yhat1; residual <- y - g
    l <- (1 - d) / (1 - p); lc <- d / p
    q <- d * (1 - e) / (e * (1 - p))
    theta <- mean((q - l) * residual) / mean(l)
    sigma2 <- mean(lc * residual^2) / mean(lc)
    h <- (1 - e) * p / (e * (1 - p)); v <- mean(2 * l - lc)
    nu2 <- mean(2 * l * h - lc * h^2) / v
    ps <- list(psi.theta.s = ((q - l) * residual - theta * l) / mean(l),
      psi.sigma2.s = lc * (residual^2 - sigma2) / mean(lc),
      psi.nu2.s = (2 * l * h - lc * h^2 - (2 * l - lc) * nu2) / v)
    got <- m$results$main$untr[[k]]
    expect_equal(unname(unlist(got$estimates[c("theta.s", "sigma2.s", "nu2.s")])),
      c(theta, sigma2, nu2), tolerance = 1e-12)
    expect_equal(got$psis[names(ps)], ps, tolerance = 1e-12)
    expect_equal(got$psis$psi.S2, nu2 * ps$psi.sigma2.s + sigma2 * ps$psi.nu2.s)
    mirror <- att.npm.cond(y, 1 - d, yhat0 = g, dhat = 1 - e,
      phat = 1 - p, trim = 0, centered = FALSE)
    expect_equal(got$estimates$theta.s, -mirror$estimates$theta.s)
    expect_equal(got$psis$psi.theta.s, -mirror$psis$psi.theta.s)
    for (key in c("sigma2.s", "nu2.s", "S2")) {
      expect_equal(got$estimates[[key]], mirror$estimates[[key]])
      expect_equal(got$psis[[paste0("psi.", key)]], mirror$psis[[paste0("psi.", key)]])
    }
    expect_equal(got$estimates$cov.theta.S2, -mirror$estimates$cov.theta.S2)
  }
})

test_that("ATU benchmarks preserve target, folds, repetitions and gain formulas", {
  m <- dml_from_did(atu_source(), 3, 4, atu_args(2, TRUE), method = "dml_refit", target = "atu")
  m <- unserialize(serialize(m, NULL))
  X <- m$adapter_data$X[, colnames(m$adapter_data$X) != "z", drop = FALSE]
  reduced <- .did_benchmark_refit(m, X)
  expect_identical(reduced$info$target, "atu")
  expect_identical(reduced$adapter_info$fold_id, m$adapter_info$fold_id)
  expect_identical(reduced$adapter_info$rep_seeds, m$adapter_info$rep_seeds)
  expect_true(reduced$info$dirty.tuning)
  b <- dml_benchmark(m, "z")
  get <- function(obj, key) vapply(obj$results$main$untr, function(r) r$estimates[[key]], numeric(1))
  nu <- get(m, "nu2.s"); nr <- get(reduced, "nu2.s")
  sg <- get(m, "sigma2.s"); sr <- get(reduced, "sigma2.s")
  expect_equal(b$benchmarks$z$gain.D, pmax(0, (nu - nr) / nr))
  expect_equal(b$benchmarks$z$gain.Y, pmax(0, (sr - sg) / sg))
  expect_equal(b$benchmarks$z$delta, get(m, "theta.s") - get(reduced, "theta.s"))
  delta <- get(m, "theta.s") - get(reduced, "theta.s")
  valid <- sr > sg & nu > nr
  rho <- rep(0, length(delta))
  rho[valid] <- sign(delta[valid]) * pmin(1,
    abs(delta[valid]) / sqrt((sr[valid] - sg[valid]) * (nu[valid] - nr[valid])))
  expect_equal(b$benchmarks$z$rho, rho)
  expect_s3_class(summary(b), "summary_dml_benchmark")
  expect_s3_class(dml_bounds(m, cf.y = 0.04, cf.d = 0.03), "dml.bounds")
  expect_true(all(is.finite(robustness_value(m))))
  expect_true(all(is.finite(extreme_robustness_value(m))))
  expect_s3_class(benchmark_bounds(m, b), "dml_benchmark_bounds")
  bad <- reduced; bad$info$target <- "att"
  expect_error(.panel_check_benchmark_alignment(m, bad), "Benchmark units")
  bad <- reduced; bad$results$main$untr <- bad$results$main$untr[1]
  expect_error(.panel_check_benchmark_alignment(m, bad), "Benchmark units")
})

test_that("ATU batch conversion retains target and applies sensitivity methods", {
  fit <- atu_source(control = "notyettreated")
  cells <- data.frame(group = c(3, 3, 5), time = c(3, 4, 5))
  batch <- dml_from_did_cells(fit, cells, atu_args(), method = "dml_refit", target = "atu")
  expect_identical(batch$target, "atu")
  expect_output(print(batch), "ATU cells")
  expect_true(all(batch$status$status == "ok"))
  expect_equal(vapply(batch$models, function(m) m$adapter_info$target_population$n, numeric(1)),
    setNames(c(240, 240, 120), names(batch$models)))
  for (i in seq_len(nrow(cells))) {
    single <- dml_from_did(fit, cells$group[i], cells$time[i], atu_args(), method = "dml_refit", target = "atu")
    expect_identical(batch$models[[i]]$coefs, single$coefs)
    expect_identical(batch$models[[i]]$info$target, "atu")
  }
  bounds <- did_cell_apply(batch, dml_bounds, cf.y = 0.04, cf.d = 0.03)
  expect_true(all(bounds$status$status == "ok"))
})

test_that("ATU guards incompatible modes, targets and pre-trend extrapolation", {
  for (method in c("dml", "replicate", "crossfit", "parametric")) {
    expect_error(dml_from_did(NULL, 3, 4, method = method, target = "atu"), "requires method = 'dml_refit'")
    expect_error(dml_from_did_cells(NULL, method = method, target = "atu"), "requires method = 'dml_refit'")
  }
  fit <- atu_source()
  expect_error(dml_from_did(fit, 3, 4, list(target = "att"), method = "dml_refit", target = "atu"), "target is fixed")
  expect_error(dml_from_did_cells(fit, dml_args = list(target = "atu"), method = "dml_refit"), "target is fixed")
  m <- dml_from_did(fit, 3, 4, atu_args(), method = "dml_refit", target = "atu")
  expect_error(did_pretrend_bounds(m, c(0, 1)), "not implemented for ATU")
  expect_error(pretrend_bounds(m, m), "not implemented for ATU")
  expect_error(dml_from_did(fit, 3, 4, modifyList(atu_args(), list(ps.trim = 0.49)),
    method = "dml_refit", target = "atu"), "clipping is active")
})

test_that("ATU preflight checks the treated outcome arm", {
  m <- dml_from_did(atu_source(), 3, 4, atu_args(), method = "dml_refit", target = "atu")
  s <- m$adapter_data
  treated <- s$D == 1
  s$y0[treated] <- 0; s$y1[treated] <- 2
  expect_error(.panel_dml_fit(s, atu_args(), "atu"), "treated-outcome variation")
  s <- m$adapter_data
  s$y0[!treated] <- 0; s$y1[!treated] <- 2
  ok <- .panel_dml_fit(s, atu_args(), "atu")
  expect_identical(ok$fit$info$target, "atu")
  expect_error(.panel_dml_fit(s, atu_args(), "att"), "control-outcome variation")
})

test_that("ATU accepts supported source estimators and rejects unequal weights", {
  for (est in c("reg", "ipw")) {
    fit <- atu_source(est = est)
    m <- dml_from_did(fit, 3, 4, atu_args(), method = "dml_refit", target = "atu")
    expect_identical(m$info$target, "atu")
    expect_null(m$adapter_info$cell$if_rows)
  }
  dat <- atu_panel(); dat$w <- 1 + (dat$id %% 3) / 10
  fit <- atu_source(dat, weightsname = "w")
  expect_error(dml_from_did(fit, 3, 4, atu_args(), method = "dml_refit", target = "atu"), "equal sampling weights")
})

test_that("default ATT remains identical to explicit ATT", {
  fit <- atu_source()
  old <- dml_from_did(fit, 3, 4, atu_args(), method = "dml_refit")
  explicit <- dml_from_did(fit, 3, 4, atu_args(), method = "dml_refit", target = "att")
  expect_identical(old$results, explicit$results)
  expect_identical(old$coefs, explicit$coefs)
  expect_identical(old$fits, explicit$fits)
  expect_output(print(old), "DML refit did ATT")
})
