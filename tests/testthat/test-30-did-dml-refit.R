refit_panel <- function() {
  skip_if_not_installed("did")
  set.seed(41102)
  n <- 360L
  u <- data.frame(id = sample(seq_len(n) * 7),
    g = rep(c(0, 3, 5), each = n / 3), x = runif(n, -1, 1),
    z = runif(n, -1, 1), u = rnorm(n))
  do.call(rbind, lapply(0:5, function(tt) transform(u, time = tt,
    x = x + 0.02 * tt * z,
    y = u + tt * (0.4 * x + 0.3 * z + 0.3 * (g > 0)) +
      0.8 * (g > 0 & tt >= g) + rnorm(n))))
}

refit_source <- function(dat = refit_panel(), fast = TRUE, est = "dr",
                         control = "nevertreated", anticipation = 0,
                         base = "varying", weightsname = NULL, save.if = FALSE) {
  did::att_gt(yname = "y", tname = "time", idname = "id", gname = "g",
    xformla = ~ x + z, data = dat, panel = TRUE, bstrap = FALSE, cband = FALSE,
    compute_inffunc = save.if, faster_mode = fast, est_method = est,
    control_group = control, anticipation = anticipation, base_period = base,
    weightsname = weightsname)
}

refit_args <- function(reps = 2, dirty = FALSE) {
  list(yreg = "lm", dreg = "glm", d.class = TRUE, cf.folds = 3,
    cf.reps = reps, cf.seed = 17, ps.trim = 0, dirty.tuning = dirty)
}

refit_native <- function(m, x = m$adapter_data$X) {
  s <- m$adapter_data
  do.call(dml, c(list(y = s$y1 - s$y0, d = s$D, x = x), m$refit_spec$dml_args))
}

expect_native_refit <- function(m, native) {
  expect_identical(m$coefs, native$coefs)
  for (i in seq_along(native$results$main$treat)) {
    r <- m$results$main$treat[[i]]
    stripped <- r; stripped$joint.covariance <- NULL
    expect_identical(stripped, native$results$main$treat[[i]])
    expect_identical(m$fits[[i]]$preds, native$fits[[i]]$preds)
    psi <- cbind(theta.s = r$psis$psi.theta.s, sigma2.s = r$psis$psi.sigma2.s,
                  nu2.s = r$psis$psi.nu2.s)
    expect_equal(colMeans(psi), setNames(rep(0, 3), colnames(psi)), tolerance = 1e-10)
    expect_equal(r$joint.covariance, crossprod(psi) / nrow(psi)^2, tolerance = 1e-12)
    expect_equal(unname(diag(r$joint.covariance)),
      unname(unlist(r$estimates[c("se.theta.s", "se.sigma2.s", "se.nu2.s")]))^2)
  }
  expect_null(m$auxiliary_dml_results)
  expect_null(m$refit_spec$source)
  expect_identical(m$adapter_info$method, "dml_refit")
}

test_that("dml_refit uses saved cell data and all native DML estimates and IFs", {
  dat <- refit_panel()
  for (fast in c(TRUE, FALSE)) {
    for (control in c("nevertreated", "notyettreated")) {
      fit <- refit_source(dat, fast, control = control, anticipation = 1,
                           base = "universal")
      m <- dml_from_did(fit, 3, 3, refit_args(), method = "dml_refit")
      ids <- sort(unique(dat$id[dat$g %in% if (control == "nevertreated") c(0, 3) else c(0, 3, 5)]))
      before <- dat[dat$time == 1, ]; after <- dat[dat$time == 3, ]
      before <- before[match(ids, before$id), ]; after <- after[match(ids, after$id), ]
      X <- model.matrix(~ x + z, before); rownames(X) <- NULL
      attr(X, "assign") <- NULL
      expect_identical(m$adapter_data$unit_id, ids)
      expect_equal(m$adapter_data$X, X)
      expect_equal(m$data$y, after$y - before$y)
      expect_equal(m$adapter_data$D, as.numeric(before$g == 3))
      expect_equal(m$adapter_info$cell$base_period, 1)
      expect_null(m$adapter_info$cell$if_rows)
      expect_null(m$adapter_info$cell$if_scale)
      expect_native_refit(m, refit_native(m))
      expect_output(print(m), "Method: dml_refit; tuning: clean")
    }
  }
})

test_that("source ATT, IF and source nuisance replay are unnecessary", {
  fit <- refit_source(save.if = TRUE)
  m <- dml_from_did(fit, 3, 4, refit_args(1), method = "dml_refit")
  bad <- fit; bad$att[] <- NA_real_; bad$inffunc <- "not an IF"
  local_mocked_bindings(.panel_check_source = function(...) stop("source replay called"),
                        .drdid_refit_source = function(...) stop("source refit called"))
  fresh <- dml_from_did(bad, 3, 4, refit_args(1), method = "dml_refit")
  expect_identical(fresh$results, m$results)
  expect_identical(fresh$coefs, m$coefs)
  expect_false(m$adapter_info$diagnostics$source.replay)
  expect_s3_class(.did_benchmark_refit(fresh, fresh$data$x[, 1:2]), "dml_did")
  expect_s3_class(did_pretrend_bounds(fresh, c(0, 1)), "dml_pretrend_bounds")
  for (fast in c(TRUE, FALSE)) for (est in c("reg", "ipw")) {
    alternative <- refit_source(fast = fast, est = est)
    expect_null(alternative$inffunc)
    alt <- dml_from_did(alternative, 3, 4, refit_args(1), method = "dml_refit")
    expect_equal(alt$results, m$results, tolerance = 1e-12)
  }
})

test_that("covariate benchmarks retain the full DML estimator after serialization", {
  skip_if_not_installed("rpart")
  args <- refit_args(dirty = TRUE)
  args$yreg <- list(method = "rpart", tuneGrid = data.frame(cp = c(0.01, 0.1)),
    trControl = caret::trainControl(method = "cv", number = 3, allowParallel = FALSE))
  m <- dml_from_did(refit_source(), 3, 4, args, method = "dml_refit")
  m <- unserialize(serialize(m, NULL)); m$call <- quote(stop("missing original call"))
  m$external_fit <- NULL
  reduced <- m$adapter_data$X[, c("(Intercept)", "x"), drop = FALSE]
  r <- .did_benchmark_refit(m, reduced)
  expect_native_refit(r, refit_native(r))
  expect_identical(r$adapter_info$fold_id, m$adapter_info$fold_id)
  expect_identical(r$adapter_info$rep_seeds, m$adapter_info$rep_seeds)
  expect_true(r$info$dirty.tuning)
  expect_equal(nrow(r$refit_spec$dml_args$yreg$tuneGrid), 2)
  expect_null(r$external_fit)
  b <- dml_benchmark(m, "z")
  expect_equal(b$benchmarks$z$theta.sj,
    vapply(r$results$main$treat, function(v) v$estimates$theta.s, numeric(1)))
  for (i in 1:2) expect_identical(b$benchmarks_psis$z$psi.theta.s.wo[[i]],
                                  r$results$main$treat[[i]]$psis$psi.theta.s)
  override <- .did_benchmark_refit(m, reduced, yreg = "lm")
  expect_identical(override$refit_spec$dml_args$yreg, "lm")
  expect_native_refit(override, refit_native(override))
})

test_that("pretrend DML uses common units/folds and full six-parameter covariance", {
  m <- dml_from_did(refit_source(control = "notyettreated"), 3, 4,
                    refit_args(dirty = TRUE), method = "dml_refit")
  b <- did_pretrend_bounds(m, c(0, 1), k = c(0, 0.5, 1), type = "scaled")
  pre <- b$pre_model
  expect_identical(b$model$results, m$results)
  expect_native_refit(pre, refit_native(pre))
  expect_identical(pre$adapter_data$D, m$adapter_data$D)
  expect_identical(pre$adapter_info$fold_id, m$adapter_info$fold_id)
  expect_identical(pre$adapter_info$rep_seeds, m$adapter_info$rep_seeds)
  expect_true(pre$info$dirty.tuning)
  expect_identical(pre$adapter_data$input.info$times, c(0, 1))
  dat <- as.data.frame(m$external_fit$DIDparams$data)
  expect_equal(pre$data$y,
    dat$y[pre$adapter_data$input.info$source_rows$post] -
    dat$y[pre$adapter_data$input.info$source_rows$pre])
  keys <- c("theta.s", "sigma2.s", "nu2.s")
  for (i in 1:2) {
    postr <- m$results$main$treat[[i]]; prer <- pre$results$main$treat[[i]]
    P <- sapply(keys, function(key) postr$psis[[paste0("psi.", key)]])
    Q <- sapply(keys, function(key) prer$psis[[paste0("psi.", key)]])
    V <- cov(cbind(P, Q)) * (nrow(P) - 1) / nrow(P)^2
    expect_equal(unname(b$joint.covariance[[i]]), unname(V), tolerance = 1e-12)
    e <- unlist(postr$estimates[keys]); f <- unlist(prer$estimates[keys])
    ratio <- sqrt(e[2] * e[3] / (f[2] * f[3]))
    # Independent gradient of theta.post - abs(theta.pre) * S.post / S.pre.
    w <- abs(f[1]) * ratio
    grad <- c(1, -w / (2 * e[2]), -w / (2 * e[3]),
              -sign(f[1]) * ratio, w / (2 * f[2]), w / (2 * f[3]))
    expect_equal(b$results[[3]][[i]]$se.minus,
      unname(sqrt(drop(t(grad) %*% V %*% grad))), tolerance = 1e-10)
  }
  expect_equal(pretrend_bounds(m, pre, type = "bias")$bounds,
    did_pretrend_bounds(m, c(0, 1), type = "bias")$bounds)
  broken <- pre; broken$adapter_info$fold_id[1, 1] <- 99L
  expect_error(pretrend_bounds(m, broken), "partitions differ")
  broken <- pre; broken$refit_spec$method <- "dml"
  expect_error(pretrend_bounds(m, broken), "same method")
  broken <- m; broken$adapter_data$X[1, 2] <- 99
  expect_error(did_pretrend_bounds(broken, c(0, 1)), "no longer matches")
  expect_error(did_pretrend_bounds(m, c(1, 3)), "precede treatment")
})

test_that("batch refits expose failures and work with downstream sensitivity tools", {
  fit <- refit_source(); fit$att[] <- NA_real_
  cells <- data.frame(group = c(3, 5, 999), time = c(4, 5, 999))
  b <- dml_from_did_cells(fit, cells, refit_args(1), method = "dml_refit")
  expect_identical(b$status$status, c("ok", "ok", "error"))
  expect_match(b$status$message[3], "absent")
  out <- did_cell_apply(b, did_pretrend_bounds, pre_period = c(0, 1))
  expect_identical(out$status$status, b$status$status)
  ci <- did_cell_apply(b, confidence_bounds, cf.y = 0.05, cf.d = 0.05)
  expect_identical(ci$status$status, b$status$status)
  m <- b$models[[1]]; native <- refit_native(m)
  expect_equal(confidence_bounds(m, cf.y = 0.05, cf.d = 0.05),
    confidence_bounds(native, cf.y = 0.05, cf.d = 0.05))
  expect_equal(confidence_bounds(m, cf.y = 0, cf.d = 0),
    confidence_bounds(native, cf.y = 0, cf.d = 0))
  expect_s3_class(sensemakr(m, cf.y = 0.05, cf.d = 0.05), "dml.sensemakr")
  expect_true(all(is.finite(robustness_value(m))))
  expect_error(dml_from_did_cells(fit, cells, refit_args(1), method = "dml_refit",
    on_error = "stop"), "absent")
})

test_that("refit enforces overlap, weight, model and sample requirements", {
  fit <- refit_source()
  args <- refit_args(1); args$ps.trim <- 0.49
  expect_error(dml_from_did(fit, 3, 4, args, method = "dml_refit"), "DML propensity clipping")
  batch <- dml_from_did_cells(fit, data.frame(group = 3, time = 4), args,
    method = "dml_refit")
  expect_identical(batch$status$status, "error")
  expect_match(batch$status$message, "DML propensity clipping")
  args <- refit_args(1); args$target <- "ate"
  expect_error(dml_from_did(fit, 3, 4, args, method = "dml_refit"), "fixed by the panel adapter")
  expect_error(dml_from_did_cells(fit, dml_args = list(dirty.tuning = NA),
    method = "dml_refit"), "dirty.tuning")
  dat <- refit_panel(); dat$w <- 1 + dat$id / max(dat$id)
  weighted <- refit_source(dat, weightsname = "w")
  expect_error(dml_from_did(weighted, 3, 4, refit_args(1), method = "dml_refit"), "equal sampling weights")
  for (field in c("panel", "est_method", "clustervars")) {
    bad <- fit
    bad$DIDparams[[field]] <- switch(field, panel = FALSE,
      est_method = function(...) NULL, clustervars = "g")
    expect_error(dml_from_did(bad, 3, 4, refit_args(1), method = "dml_refit"),
      switch(field, panel = "balanced panel", est_method = "custom", clustervars = "Clustering"))
  }
  bad <- fit; bad$DIDparams$data <- as.data.frame(bad$DIDparams$data)[-1, ]
  expect_error(dml_from_did(bad, 3, 4, refit_args(1), method = "dml_refit"), "balanced panel")
  expect_error(dml_from_drdid(NULL, method = "dml_refit"), "available through dml_from_did")
})

test_that("active source trimming is not inherited by full DML refitting", {
  skip_if_not_installed("did")
  set.seed(31032)
  n <- 6010L
  u <- data.frame(id = seq_len(n), x = c(rep(1, 4010), rep(0, 2000)),
    D = c(rep(1, 4000), rep(0, 10), rep(1, 1000), rep(0, 1000)),
    z = runif(n, -1, 1), u = rnorm(n))
  dat <- do.call(rbind, lapply(0:3, function(tt) transform(u, time = tt, g = 3 * D,
    y = u + tt * (0.4 * x + 0.3 * z + 0.2 * D) +
      0.8 * D * (tt >= 3) + rnorm(n))))
  fit <- refit_source(dat, save.if = TRUE)
  expect_error(dml_from_did(fit, 3, 3, refit_args(1)), "Source propensity trimming")
  m <- dml_from_did(fit, 3, 3, refit_args(1), method = "dml_refit")
  expect_equal(nrow(m$data$x), n)
  expect_false(m$adapter_info$diagnostics$source.trimming.reused)
  expect_native_refit(m, refit_native(m))
})
