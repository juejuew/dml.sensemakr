method_fixture <- function() {
  skip_if_not_installed("DRDID")
  set.seed(110108)
  n <- 240L
  X <- cbind(`(Intercept)` = 1, x = runif(n, -1, 1), z = runif(n, -1, 1))
  D <- rbinom(n, 1, plogis(0.1 + 0.3 * X[, "x"] - 0.2 * X[, "z"]))
  y0 <- rnorm(n)
  y1 <- y0 + 0.8 * D + 0.7 * X[, "x"] + 0.3 * X[, "z"] + rnorm(n)
  input <- list(y1 = y1, y0 = y0, D = D, covariates = X, unit_id = seq_len(n))
  list(input = input, source = DRDID::drdid_panel(y1, y0, D, X, inffunc = TRUE))
}

method_args <- function(dirty = FALSE, tuned = FALSE) {
  yreg <- if (tuned) list(method = "rpart", tuneGrid = data.frame(cp = c(0.01, 0.1)),
    trControl = caret::trainControl(method = "cv", number = 3, allowParallel = FALSE)) else "lm"
  list(yreg = yreg, dreg = list(method = "glm",
    trControl = caret::trainControl(method = "none", classProbs = TRUE)),
    d.class = TRUE, cf.folds = 3, cf.reps = 2, cf.seed = 17,
    ps.trim = 0, dirty.tuning = dirty, verbose = FALSE)
}

method_fit <- function(f, method = "dml", args = list()) {
  do.call(dml_from_drdid, c(list(fit = f$source, method = method, dml_args = args), f$input))
}

method_panel <- function() {
  skip_if_not_installed("did")
  set.seed(110109)
  n <- 240L
  u <- data.frame(id = seq_len(n), g = rep(c(0, 3, 4), each = n / 3),
    x = runif(n, -1, 1), z = runif(n, -1, 1), u = rnorm(n))
  do.call(rbind, lapply(1:5, function(t) transform(u, time = t,
    y = u + t * (0.4 * x + 0.2 * z + 0.5 * (g == 3)) +
      0.9 * (g > 0 & t >= g) + rnorm(n))))
}

method_did_fit <- function(dat, weightsname = NULL) {
  did::att_gt(yname = "y", tname = "time", idname = "id", gname = "g",
    xformla = ~ x + z, data = dat, weightsname = weightsname,
    panel = TRUE, bstrap = FALSE, cband = FALSE, compute_inffunc = TRUE)
}

test_that("canonical names and legacy aliases yield identical adapter estimates", {
  f <- method_fixture()
  for (pair in list(c("dml", "crossfit"), c("replicate", "parametric"))) {
    args <- if (pair[1] == "dml") method_args() else list()
    current <- method_fit(f, pair[1], args)
    legacy <- method_fit(f, pair[2], args)
    expect_identical(current$results, legacy$results)
    expect_identical(current$refit_spec$method, pair[1])
    expect_identical(legacy$refit_spec$method, pair[1])
    expect_identical(current$adapter_info$method, pair[1])
    expect_identical(legacy$adapter_info$method, pair[1])
    expect_identical(current$results$main$treat[[1]]$estimates$theta.s, as.numeric(f$source$ATT))
    expect_identical(current$results$main$treat[[1]]$psis$psi.theta.s, as.numeric(f$source$att.inf.func))
    expect_output(print(current), if (pair[1] == "dml") "Method: dml; tuning: clean" else "Replicate")
  }
  default <- do.call(dml_from_drdid, c(list(fit = f$source, dml_args = method_args()), f$input))
  expect_identical(default$refit_spec$method, "dml")
  expect_false(default$info$dirty.tuning)
})

test_that("dirty tuning validates options and never applies to replicate", {
  f <- method_fixture()
  for (bad in list(NA, 1, "TRUE", c(TRUE, FALSE), NULL, logical())) {
    args <- method_args(); args["dirty.tuning"] <- list(bad)
    expect_error(method_fit(f, "dml", args), "dirty.tuning must be TRUE or FALSE", fixed = TRUE)
  }
  for (method in c("replicate", "parametric")) {
    expect_error(method_fit(f, method, list(dirty.tuning = TRUE)), "empty dml_args")
    expect_error(method_fit(f, method, list(dirty.tuning = FALSE)), "empty dml_args")
  }
  args <- method_args(TRUE); args$target <- "ate"
  expect_error(method_fit(f, "dml", args), "fixed by the panel adapter")
  expect_error(method_fit(f, "other"), "arg.*should be one of")
})

test_that("full-sample tuning still uses fold-specific training and native DML predictions", {
  skip_if_not_installed("rpart")
  f <- method_fixture()
  args <- method_args(TRUE, TRUE); args$save.models <- TRUE
  m <- method_fit(f, "dml", args)
  direct <- do.call(dml, c(list(y = f$input$y1 - f$input$y0,
    d = f$input$D, x = f$input$covariates), m$refit_spec$dml_args))
  expect_identical(m$auxiliary_dml_results, direct$results)
  expect_true(m$info$dirty.tuning)
  expect_true(m$refit_spec$dml_args$dirty.tuning)
  expect_equal(nrow(m$refit_spec$dml_args$yreg$tuneGrid), 2L)
  expect_equal(nrow(m$yreg$yreg0$tuneGrid), 1L)
  expect_identical(m$yreg$yreg0$trControl$method, "none")
  expect_output(print(m), "tuning: dirty")
  for (r in 1:2) {
    expect_identical(m$fits[[r]]$preds, direct$fits[[r]]$preds)
    expect_identical(m$results$main$treat[[r]]$psis$psi.theta.s, as.numeric(f$source$att.inf.func))
    for (k in 1:3) {
      training <- m$adapter_info$fold_id[, r] != k
      expect_equal(nrow(m$fits[[r]]$model.d[[k]]$trainingData), sum(training))
      expect_equal(nrow(m$fits[[r]]$model.y0[[k]]$trainingData), sum(training & f$input$D == 0))
    }
  }
  restored <- unserialize(serialize(m, NULL))
  reduced <- f$input$covariates[, c("(Intercept)", "x"), drop = FALSE]
  refit <- .drdid_benchmark_refit(restored, reduced)
  expect_true(refit$info$dirty.tuning)
  expect_identical(refit$refit_spec$method, "dml")
  expect_identical(refit$adapter_info$fold_id, m$adapter_info$fold_id)
  expect_equal(nrow(refit$refit_spec$dml_args$yreg$tuneGrid), 2L)
  old <- restored; old$refit_spec$method <- "crossfit"; old$adapter_info$method <- NULL
  expect_identical(.drdid_benchmark_refit(old, reduced)$results, refit$results)
  old$refit_spec$method <- NULL
  expect_identical(.drdid_benchmark_refit(old, reduced)$results, refit$results)
  b <- dml_benchmark(restored, benchmark_covariates = "z")
  expect_equal(nrow(b$benchmarks$z), 2L)
  expect_true(all(is.finite(b$benchmarks$z$gain.Y)))
  expect_true(all(is.finite(b$benchmarks$z$gain.D)))
})

test_that("replicate saved objects with old method names remain usable", {
  f <- method_fixture()
  m <- method_fit(f, "replicate")
  old <- unserialize(serialize(m, NULL))
  old$refit_spec$method <- "parametric"; old$adapter_info$method <- NULL
  X <- f$input$covariates[, c("(Intercept)", "x"), drop = FALSE]
  expect_identical(.drdid_benchmark_refit(old, X)$results,
    .drdid_benchmark_refit(m, X)$results)
  expect_error(.drdid_benchmark_refit(old, X, dreg = "glm"), "overrides")
  # Alias normalization also permits an old pre fit paired with a new post fit.
  expect_equal(pretrend_bounds(m, old, type = "scaled")$bounds,
    pretrend_bounds(m, m, type = "scaled")$bounds)
  expect_identical(pretrend_bounds(m, old)$method, "replicate")
})

test_that("did batches and pre-trend refits inherit canonical methods and dirty tuning", {
  skip_if_not_installed("rpart")
  fit <- method_did_fit(method_panel())
  args <- method_args(TRUE, TRUE)
  cells <- data.frame(group = c(3, 4), time = c(4, 5))
  batch <- dml_from_did_cells(fit, cells, method = "dml", dml_args = args, on_error = "stop")
  expect_identical(batch$status$status, c("ok", "ok"))
  m <- batch$models[[1L]]
  expect_identical(m$adapter_info$method, "dml")
  expect_output(print(m), "tuning: dirty")
  out <- did_cell_apply(batch, did_pretrend_bounds, pre_period = c(1, 2), type = "scaled")
  expect_identical(out$status$status, c("ok", "ok"))
  for (i in 1:2) {
    b <- out$values[[i]]
    expect_identical(b$method, "dml")
    expect_true(b$pre_model$info$dirty.tuning)
    expect_true(b$pre_model$refit_spec$dml_args$dirty.tuning)
    expect_identical(b$pre_model$adapter_info$fold_id, batch$models[[i]]$adapter_info$fold_id)
    expect_equal(nrow(b$pre_model$refit_spec$dml_args$yreg$tuneGrid), 2L)
    expect_true(all(is.finite(b$bounds$lwr)))
  }
  X <- m$adapter_data$X[, 1:2, drop = FALSE]
  refit <- .did_benchmark_refit(m, X)
  expect_true(refit$info$dirty.tuning)
  expect_identical(refit$adapter_info$fold_id, m$adapter_info$fold_id)
  old <- unserialize(serialize(m, NULL)); old$refit_spec$method <- "crossfit"
  old$adapter_info$method <- NULL
  expect_equal(did_pretrend_bounds(old, c(1, 2), type = "scaled")$bounds, out$values[[1]]$bounds)
  expect_identical(.did_benchmark_refit(old, X)$results, refit$results)
  canonical <- dml_from_did_cells(fit, cells, method = "replicate", on_error = "stop")
  legacy <- dml_from_did_cells(fit, cells, method = "parametric", on_error = "stop")
  for (i in 1:2) expect_identical(canonical$models[[i]]$results, legacy$models[[i]]$results)
  mr <- canonical$models[[1L]]
  old <- unserialize(serialize(mr, NULL)); old$refit_spec$method <- "parametric"
  old$adapter_info$method <- NULL
  expect_equal(did_pretrend_bounds(old, c(1, 2))$bounds, did_pretrend_bounds(mr, c(1, 2))$bounds)
  expect_identical(.did_benchmark_refit(old, X)$results, .did_benchmark_refit(mr, X)$results)
  expect_error(dml_from_did_cells(fit, cells, method = "dml",
    dml_args = list(dirty.tuning = NA)), "dirty.tuning")
})

test_that("unequal-weight fallback passes dirty tuning through to native DML", {
  dat <- method_panel(); dat$w <- 1 + dat$id / max(dat$id)
  fit <- method_did_fit(dat, "w")
  args <- method_args(TRUE)
  expect_warning(m <- dml_from_did(fit, 3, 4, method = "dml", dml_args = args), "Unequal")
  s <- m$adapter_data
  native <- do.call(dml, c(list(y = s$y1 - s$y0, d = s$D, x = s$X), m$adapter_info$dml_args))
  expect_true(m$info$dirty.tuning)
  expect_identical(m$adapter_info$method, "dml")
  expect_identical(m$results, native$results)
  expect_error(dml_from_did(fit, 3, 4, method = "replicate"), "equal sampling weights")
})
