did_adapter_data <- function() {
  skip_if_not_installed("did")
  skip_if_not_installed("DRDID")
  set.seed(608)
  n <- 360L
  units <- data.frame(id = sample(seq_len(n) * 7 + 100),
    g = rep(c(0, 5, 9), each = n / 3), x = runif(n, -1, 1),
    z = runif(n, -1, 1), f = factor(sample(c("a", "b"), n, TRUE),
                                   levels = c("a", "b", "unused")), u = rnorm(n))
  dat <- do.call(rbind, lapply(c(0, 2, 5, 9, 12), function(t) {
    d <- units
    d$time <- t
    d$x <- d$x + 0.04 * t * d$z
    d$y <- d$u + t / 8 + (0.3 + t / 7) * d$x + 0.7 * d$z +
      (d$g > 0 & t >= d$g) * (0.8 + t / 20) + rnorm(n)
    d$w <- 2
    d
  }))
  dat[sample.int(nrow(dat)), ]
}

did_adapter_fit <- function(data, faster_mode = TRUE, control_group = "nevertreated",
                            base_period = "varying", anticipation = 0, ...) {
  did::att_gt(yname = "y", tname = "time", idname = "id", gname = "g",
              xformla = ~ x + z, data = data, panel = TRUE, bstrap = FALSE,
              cband = FALSE, faster_mode = faster_mode,
              control_group = control_group, base_period = base_period,
              anticipation = anticipation, ...)
}

did_adapter_args <- function() {
  list(yreg = "lm", dreg = "glm", d.class = TRUE,
        cf.folds = 3, cf.seed = 17, ps.trim = 0)
}

test_that("post-treatment cells preserve ATT, ID mapping and the rescaled source IF", {
  dat <- did_adapter_data()
  models <- list()
  for (fast in c(TRUE, FALSE)) {
    for (controls in c("nevertreated", "notyettreated")) {
      fit <- did_adapter_fit(dat, faster_mode = fast, control_group = controls)
      m <- dml_from_did(fit, 5, 5, did_adapter_args())
      r <- m$results$main$treat[[1]]; cell <- m$adapter_info$cell
      j <- which(fit$group == 5 & fit$t == 5)
      allowed <- if (controls == "nevertreated") c(0, 5) else c(0, 5, 9)
      expected.ids <- sort(unique(dat$id[dat$g %in% allowed]))
      expect_s3_class(m, "dml_did")
      expect_s3_class(m, "dml")
      expect_identical(m$adapter_data$unit_id, expected.ids)
      expect_identical(r$estimates$theta.s, fit$att[j])
      expect_equal(r$psis$psi.theta.s,
        as.numeric(fit$inffunc[match(as.character(expected.ids), rownames(fit$inffunc)), j]) *
          length(expected.ids) / fit$n, tolerance = 1e-12)
      expect_equal(cell$if_rows, match(as.character(expected.ids), rownames(fit$inffunc)))
      expect_equal(cell$base_period, 2)
      expect_equal(unname(se(m)), fit$se[j], tolerance = 1e-10)
      expect_equal(m$adapter_data$D,
        as.numeric(dat$g[match(expected.ids, dat$id)] == 5))
      expect_output(print(m), "Hybrid did ATT")
      expect_error(dml_gate(m, rep(1, length(expected.ids))), "external ATT")
      models[[paste(fast, controls)]] <- m
    }
  }
  for (controls in c("nevertreated", "notyettreated")) {
    a <- models[[paste(TRUE, controls)]]; b <- models[[paste(FALSE, controls)]]
    expect_equal(a$adapter_data$X, b$adapter_data$X)
    expect_identical(a$adapter_info$fold_id, b$adapter_info$fold_id)
    expect_equal(a$results$main, b$results$main, tolerance = 1e-8)
  }
})

test_that("the reconstructed cell matches inputs received by did's estimator", {
  dat <- did_adapter_data()
  for (fast in c(TRUE, FALSE)) {
    for (base in c("varying", "universal")) {
      captured <- list()
      oracle <- function(y1, y0, D, covariates, i.weights, inffunc = TRUE) {
        ans <- DRDID::drdid_panel(y1, y0, D, covariates, i.weights,
                                  inffunc = inffunc)
        captured[[length(captured) + 1L]] <<- list(ATT = ans$ATT,
          y1 = y1, y0 = y0, D = D, X = covariates, w = i.weights)
        ans
      }
      opts <- list(yname = "y", tname = "time", idname = "id", gname = "g",
        xformla = ~ poly(x, 2) + scale(z) + f, data = dat,
        bstrap = FALSE, cband = FALSE, faster_mode = fast,
        base_period = base, anticipation = 1, control_group = "notyettreated")
      original <- do.call(did::att_gt, opts)
      recorded <- do.call(did::att_gt, c(opts, list(est_method = oracle)))
      m <- dml_from_did(original, 9, 12, did_adapter_args())
      j <- which(recorded$group == 9 & recorded$t == 12)
      hit <- which(vapply(captured, function(v) abs(v$ATT - recorded$att[j]) < 1e-10, logical(1)))
      expect_length(hit, 1)
      actual <- captured[[hit]]
      # Continuous outcomes uniquely identify oracle rows, independently of
      # the adapter's ID/cohort reconstruction.
      ix <- match(m$adapter_data$y0, actual$y0)
      expect_false(anyNA(ix))
      expect_equal(m$adapter_data$y1, as.numeric(actual$y1[ix]))
      expect_equal(m$adapter_data$D, as.numeric(actual$D[ix]))
      expect_equal(unname(m$adapter_data$X), unname(actual$X[ix, , drop = FALSE]), tolerance = 1e-12)
      expect_equal(m$adapter_data$i.weights, actual$w[ix])
      expect_equal(m$adapter_info$cell$base_period, 5)
      expect_equal(unname(coef(m)), recorded$att[j], tolerance = 1e-10)
    }
  }
})

test_that("saved cleaning, row order and serialization are honored", {
  dat <- did_adapter_data()
  bad.ids <- unique(dat$id)[1:3]
  dat$y[which(dat$id == bad.ids[1])[1]] <- NA_real_
  dat$x[which(dat$id == bad.ids[2])[1]] <- Inf
  dat <- dat[-which(dat$id == bad.ids[3])[1], ]
  dat$irrelevant <- NA_real_
  for (fast in c(TRUE, FALSE)) {
    fit <- suppressWarnings(did_adapter_fit(dat, faster_mode = fast))
    path <- tempfile(fileext = ".rds")
    saveRDS(fit, path)
    restored <- readRDS(path); unlink(path)
    restored$call <- quote(stop("original data unavailable"))
    m <- dml_from_did(restored, 5, 9, did_adapter_args())
    expect_false(any(bad.ids %in% m$adapter_data$unit_id))
    expect_equal(m$adapter_info$cell$n_global, length(unique(dat$id)) - 3)
    saved <- as.data.frame(restored$DIDparams$data)
    rows <- m$adapter_data$input.info$source_rows
    expect_equal(saved$id[rows$pre], m$adapter_data$unit_id)
    expect_equal(saved$id[rows$post], m$adapter_data$unit_id)
    expect_equal(saved$y[rows$post] - saved$y[rows$pre], m$data$y)
    # Permuting complete IF rows together with their ID labels is harmless.
    set.seed(2); ix <- sample.int(restored$n)
    restored$inffunc <- restored$inffunc[ix, , drop = FALSE]
    m2 <- dml_from_did(restored, 5, 9, did_adapter_args())
    expect_equal(m2$results$main, m$results$main, tolerance = 1e-10)
  }
})

test_that("joint covariance, sensitivity and benchmark refits are coherent", {
  fit <- did_adapter_fit(did_adapter_data())
  m <- dml_from_did(fit, 5, 9, did_adapter_args())
  r <- m$results$main$treat[[1]]; e <- r$estimates
  joint <- cbind(r$psis$psi.theta.s, r$psis$psi.sigma2.s, r$psis$psi.nu2.s)
  V <- cov(joint) * (nrow(joint) - 1) / nrow(joint)^2
  expect_equal(unname(r$joint.covariance), unname(V), tolerance = 1e-12)
  k <- sqrt(0.08 * 0.06 / (1 - 0.06))
  grad <- c(1, -k / 2 * sqrt(e$nu2.s / e$sigma2.s), -k / 2 * sqrt(e$sigma2.s / e$nu2.s))
  expected.se <- sqrt(drop(t(grad) %*% V %*% grad))
  bounds <- dml_bounds(m, cf.y = 0.08, cf.d = 0.06)
  expect_equal(bounds$results$main$treat[[1]]$estimates$se.theta.m, expected.se)
  ci <- confidence_bounds(m, cf.y = 0.08, cf.d = 0.06, max = FALSE)
  expect_equal(unname(ci["att", "lwr"]), e$theta.s - k * sqrt(e$S2) - qnorm(0.95) * expected.se)
  zero <- confidence_bounds(m, cf.y = 0, cf.d = 0)
  expect_equal(unname(zero["att", ]), unname(coef(m)) + c(-1, 1) * qnorm(0.95) * unname(se(m)))
  expect_s3_class(sensemakr(m, cf.y = 0.08, cf.d = 0.06), "dml.sensemakr")
  expect_true(all(is.finite(robustness_value(m))))
  expect_s3_class(plot(m), "ggplot")
  path <- tempfile(fileext = ".rds"); saveRDS(m, path)
  m <- readRDS(path); unlink(path)
  m$call <- quote(stop("missing original call")); m$external_fit <- NULL
  bench <- dml_benchmark(m, "z")
  s <- m$adapter_data
  reduced <- DRDID::drdid_panel(s$y1, s$y0, s$D, s$X[, c("(Intercept)", "x")],
                               i.weights = s$i.weights, inffunc = TRUE)
  expect_equal(bench$benchmarks$z$theta.sj, as.numeric(reduced$ATT))
  expect_equal(bench$benchmarks_psis$z$psi.theta.s.wo[[1]], as.numeric(reduced$att.inf.func))
  refit <- dml.sensemakr:::.did_benchmark_refit(m, s$X[, c("(Intercept)", "x")])
  expect_s3_class(refit, "dml_did")
  expect_identical(refit$adapter_info$unit_id, m$adapter_info$unit_id)
  expect_identical(refit$adapter_info$fold_id, m$adapter_info$fold_id)
  expect_equal(bench$benchmarks_psis$z$psi.sigma2.s.wo[[1]],
                refit$results$main$treat[[1]]$psis$psi.sigma2.s)
  expect_s3_class(benchmark_bounds(m, bench, kY = 0.5, kD = 0.5), "dml_benchmark_bounds")
})

test_that("batch operations retain failed cells, warnings and successful NULL results", {
  fit <- did_adapter_fit(did_adapter_data(), base_period = "universal")
  selected <- data.frame(group = c(5, 5, 5, 999), time = c(5, 2, 12, 999))
  batch <- dml_from_did_cells(fit, selected, did_adapter_args())
  expect_identical(batch$status$status, c("ok", "error", "ok", "error"))
  expect_length(batch$models, 4)
  expect_identical(names(batch$models), batch$status$cell)
  expect_null(batch$models[[2]])
  expect_match(batch$status$message[2], "post-treatment")
  expect_match(batch$status$message[4], "absent")
  expect_output(print(batch), "error")
  result <- did_cell_apply(batch, confidence_bounds, cf.y = 0.05, cf.d = 0.05)
  expect_identical(result$status$status, batch$status$status)
  expect_length(result$values, 4)
  expect_equal(result$values[[1]], confidence_bounds(batch$models[[1]], cf.y = 0.05, cf.d = 0.05))
  fun <- function(m) { warning("kept warning"); if (m$adapter_info$cell$time == 12) stop("cell failure"); NULL }
  result <- did_cell_apply(batch, fun)
  expect_length(result$values, 4)
  expect_null(result$values[[1]])
  expect_identical(result$status$status, c("ok", "error", "error", "error"))
  expect_match(result$status$warnings[1], "kept warning")
  expect_match(result$status$message[3], "cell failure")
  expect_error(did_cell_apply(batch, coef, on_error = "stop"), "post-treatment")
  expect_error(dml_from_did_cells(fit, selected, did_adapter_args(), on_error = "stop"), "post-treatment")
  all <- dml_from_did_cells(fit, dml_args = did_adapter_args())
  expect_true(all(all$status$status == "ok"))
  expect_equal(nrow(all$status), sum(fit$t >= fit$group))
  expect_equal(nrow(all$excluded), sum(fit$t < fit$group))
  expect_true(all(all$excluded$time < all$excluded$group))
  expect_s3_class(did_cell_apply(all, dml_benchmark, "z"), "did_cell_results")
  bad <- fit; j <- which(bad$group == 5 & bad$t == 5); bad$att[j] <- NA_real_
  partial <- dml_from_did_cells(bad, dml_args = did_adapter_args())
  expect_equal(sum(partial$status$status == "error"), 1)
  expect_equal(length(partial$models), nrow(all$status))
})

test_that("unsupported source paths and inconsistent inputs fail explicitly", {
  dat <- did_adapter_data(); fit <- did_adapter_fit(dat)
  run <- function(x = fit, args = did_adapter_args()) dml_from_did(x, 5, 9, args)
  bad <- fit; bad$inffunc <- NULL
  expect_error(run(bad), "compute_inffunc = TRUE")
  bad <- fit; rownames(bad$inffunc) <- NULL
  expect_error(run(bad), "unique row names")
  bad <- fit; rownames(bad$inffunc)[2] <- rownames(bad$inffunc)[1]
  expect_error(run(bad), "unique row names")
  bad <- fit; bad$inffunc <- bad$inffunc[-1, ]
  expect_error(run(bad), "dimensions")
  j <- which(fit$group == 5 & fit$t == 9)
  bad <- fit; bad$inffunc[, j] <- rev(as.numeric(bad$inffunc[, j]))
  expect_error(run(bad), "outside|ordered IF")
  bad <- fit; bad$att[j] <- bad$att[j] + 0.1
  expect_error(run(bad), "ordered IF")
  bad <- fit; bad$DIDparams$data <- as.data.frame(bad$DIDparams$data)
  bad$DIDparams$data$y[1] <- bad$DIDparams$data$y[1] + 1
  expect_error(run(bad), "tensors disagree")
  bad <- fit; bad$DIDparams$data <- as.data.frame(bad$DIDparams$data)[-1, ]
  expect_error(run(bad), "balanced panel")
  bad <- fit; bad$DIDparams$panel <- FALSE
  expect_error(run(bad), "balanced panel estimator")
  for (method in list("reg", "ipw", function(...) NULL)) {
    bad <- fit; bad$DIDparams$est_method <- method
    expect_error(run(bad), "est_method")
  }
  bad <- fit; bad$DIDparams$clustervars <- "g"
  expect_error(run(bad), "Clustering")
  bad <- fit; class(bad) <- "AGGTEobj"
  expect_error(run(bad), "aggte")
  expect_error(dml_from_did(fit, "5", 9), "numeric scalars")
  expect_error(dml_from_did(fit, 5, 2), "post-treatment")
  expect_error(dml_from_did_cells(fit, data.frame(group = c(5, 5), time = c(9, 9))), "unique")
  args <- did_adapter_args(); args$cf.reps <- 0
  expect_error(run(args = args), "cf.reps must be an integer")
  args <- did_adapter_args(); args$ps.trim <- 0.49
  expect_error(run(args = args), "DML propensity clipping")
  for (fast in c(TRUE, FALSE)) {
    no.if <- did_adapter_fit(dat, faster_mode = fast, compute_inffunc = FALSE)
    expect_error(run(no.if), "compute_inffunc = TRUE")
    varying <- did_adapter_fit(dat, faster_mode = fast, fix_weights = "varying")
    expect_error(run(varying), "repeated-cross-section")
    w <- dat; w$w <- 1 + w$id / max(w$id)
    weighted <- did_adapter_fit(w, faster_mode = fast, weightsname = "w")
    expect_error(run(weighted), "constant positive")
  }
})

test_that("equal weights, fixed weight periods and intercept-only defaults work", {
  skip_if_not_installed("ranger")
  dat <- did_adapter_data()
  for (fast in c(TRUE, FALSE)) {
    for (fix in c("first_period", "base_period")) {
      fit <- did_adapter_fit(dat, faster_mode = fast, weightsname = "w", fix_weights = fix)
      m <- dml_from_did(fit, 5, 9, did_adapter_args())
      expect_true(all(m$adapter_data$i.weights == 1))
      expect_true(m$adapter_info$diagnostics$reproduced)
    }
    fit <- did::att_gt(yname = "y", tname = "time", idname = "id", gname = "g",
      xformla = ~ 1, data = dat, bstrap = FALSE, cband = FALSE, faster_mode = fast)
    m <- dml_from_did(fit, 5, 9)
    expect_true(is.finite(as.numeric(coef(m))))
    expect_identical(colnames(m$adapter_data$X), "(Intercept)")
    expect_false(m$adapter_info$diagnostics$trim.active)
  }
})

test_that("source trimming is rejected even when did returns a finite ATT", {
  did_adapter_data() # dependency gates
  set.seed(132)
  n <- 1000L
  u <- data.frame(id = seq_len(n), x = rep(c(1, 0), c(801, 199)),
    g = c(rep(1, 800), 0, rep(1, 100), rep(0, 99)), y = rnorm(n))
  dat <- rbind(transform(u, time = 0),
               transform(u, time = 1, y = y + x + 0.7 * (g == 1) + rnorm(n)))
  for (fast in c(TRUE, FALSE)) {
    expect_warning(fit <- did::att_gt(yname = "y", tname = "time", idname = "id", gname = "g",
      xformla = ~ x, data = dat, bstrap = FALSE, cband = FALSE, faster_mode = fast),
      "No pre-treatment periods")
    expect_true(is.finite(fit$att[1]))
    expect_error(dml_from_did(fit, 1, 1, did_adapter_args()), "Source propensity trimming")
  }
})

test_that("anticipation cutoffs and no-never-treated preprocessing follow did", {
  dat <- did_adapter_data()
  for (fast in c(TRUE, FALSE)) {
    fit <- did_adapter_fit(dat, faster_mode = fast, control_group = "notyettreated", anticipation = 4)
    m <- dml_from_did(fit, 5, 5, did_adapter_args())
    expect_equal(m$adapter_info$cell$base_period, 0)
    expect_equal(m$adapter_data$unit_id, sort(unique(dat$id[dat$g %in% c(0, 5)])))
    fit <- suppressWarnings(did_adapter_fit(dat[dat$g != 0, ], faster_mode = fast,
                                             control_group = "notyettreated"))
    m <- dml_from_did(fit, 5, 5, did_adapter_args())
    expect_equal(m$adapter_data$unit_id, sort(unique(dat$id[dat$g != 0])))
    expect_true(m$adapter_info$diagnostics$reproduced)
  }
})

test_that("bootstrap sources use joint analytic IF inference and default learners run", {
  skip_if_not_installed("ranger")
  dat <- did_adapter_data()
  set.seed(31)
  fit <- did::att_gt(yname = "y", tname = "time", idname = "id", gname = "g",
    xformla = ~ x + z, data = dat, bstrap = TRUE, biters = 99, cband = TRUE)
  m <- dml_from_did(fit, 5, 9)
  cell <- m$adapter_info$cell
  expect_true(m$adapter_info$original.boot)
  expect_true(m$adapter_info$original.cband)
  expect_identical(m$adapter_info$original.se, fit$se[cell$column])
  phi <- as.numeric(fit$inffunc[, cell$column])
  expect_equal(unname(se(m)), sqrt(sum((phi - mean(phi))^2)) / fit$n, tolerance = 1e-10)
  expect_equal(unname(coef(m)), fit$att[cell$column])
  expect_true(is.finite(m$results$main$treat[[1]]$estimates$S2))
})
