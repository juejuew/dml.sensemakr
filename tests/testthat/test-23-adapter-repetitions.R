repeated_adapter_fixture <- function(method = "trad", reps = 3L) {
  skip_if_not_installed("DRDID")
  set.seed(421)
  n <- 360L
  X <- cbind(`(Intercept)` = 1, x1 = runif(n, -1, 1), x2 = runif(n, -1, 1))
  D <- rbinom(n, 1, plogis(0.1 + 0.5 * X[, "x1"] + 0.6 * X[, "x2"]))
  y0 <- X[, "x1"] + rnorm(n)
  y1 <- y0 + 0.7 * D + X[, "x1"] + 1.5 * X[, "x2"] + rnorm(n)
  args <- list(yreg = "lm", dreg = "glm", d.class = TRUE,
                cf.folds = 3, cf.reps = reps, cf.seed = 17, ps.trim = 0)
  if (method == "did") {
    skip_if_not_installed("did")
    G <- ifelse(D == 1, 2, ifelse(seq_len(n) %% 2 == 0, 3, 0))
    dat <- do.call(rbind, lapply(0:3, function(t) data.frame(
      id = seq_len(n) * 11, time = t, g = G,
      y = y0 + t / 2 * (y1 - y0) + rnorm(n, sd = 0.2),
      x1 = X[, "x1"], x2 = X[, "x2"])))
    source <- did::att_gt(yname = "y", tname = "time", idname = "id", gname = "g",
      xformla = ~ x1 + x2, data = dat, bstrap = FALSE, cband = FALSE)
    model <- dml_from_did(source, 2, 2, dml_args = args)
    cell <- model$adapter_info$cell
    theta <- source$att[cell$column]
    phi <- as.numeric(source$inffunc[cell$if_rows, cell$column]) * cell$if_scale
  } else {
    fun <- if (method == "imp") DRDID::drdid_imp_panel else DRDID::drdid_panel
    source <- fun(y1, y0, D, X, inffunc = TRUE)
    model <- dml_from_drdid(source, y1 = y1, y0 = y0, D = D,
                            covariates = X, dml_args = args)
    theta <- as.numeric(source$ATT); phi <- as.numeric(source$att.inf.func)
  }
  list(source = source, model = model, theta = theta, phi = phi, args = args)
}

repeated_direct_dml <- function(model, X = model$adapter_data$X, args = model$refit_spec$dml_args) {
  s <- model$adapter_data
  do.call(dml, c(list(y = s$y1 - s$y0, d = s$D, x = X), args))
}

test_that("cf.reps is validated without discarding or rounding the request", {
  args <- dml.sensemakr:::.drdid_dml_args
  expect_equal(args(list(), 10)$cf.reps, 1)
  expect_equal(args(list(cf.reps = 3), 10)$cf.reps, 3)
  for (bad in list(NULL, 0, -1, 1.5, NA_real_, NaN, Inf, TRUE, "3", c(1, 2), 1000001))
    expect_error(args(list(cf.reps = bad), 10), "cf.reps must be an integer")
})

test_that("did and both DRDID methods retain every repetition and the external ATT IF", {
  for (method in c("did", "trad", "imp")) {
    f <- repeated_adapter_fixture(method)
    m <- f$model; s <- m$adapter_data; n <- length(s$D)
    pure <- repeated_direct_dml(m)
    expect_length(m$results$main$treat, 3)
    expect_length(m$fits, 3)
    expect_length(m$adapter_info$dml.att, 3)
    expect_equal(m$info$cf.reps, 3)
    expect_equal(dim(m$adapter_info$fold_id), c(n, 3))
    expect_identical(m$fits, pure$fits)
    expect_equal(m$auxiliary_dml_results, pure$results)
    set.seed(f$args$cf.seed)
    expect_identical(m$adapter_info$rep_seeds, sample.int(1e6, 3))
    expect_false(identical(m$adapter_info$fold_id[, 1], m$adapter_info$fold_id[, 2]))
    for (i in 1:3) {
      r <- m$results$main$treat[[i]]; original <- pure$results$main$treat[[i]]
      expect_identical(r$estimates$theta.s, f$theta)
      expect_identical(r$psis$psi.theta.s, f$phi)
      expect_identical(r$psis$psi.sigma2.s, original$psis$psi.sigma2.s)
      expect_identical(r$psis$psi.nu2.s, original$psis$psi.nu2.s)
      expect_identical(m$adapter_info$dml.att[i], original$estimates$theta.s)
      joint <- cbind(f$phi, original$psis$psi.sigma2.s, original$psis$psi.nu2.s)
      expected.V <- cov(joint) * (n - 1) / n^2
      expect_equal(unname(r$joint.covariance), unname(expected.V), tolerance = 1e-12)
      # Cross-fitting returns the training-arm proportion for each held-out
      # unit; this independently checks metadata against the executed folds.
      folds <- m$adapter_info$fold_id[, i]
      expected.phat <- vapply(seq_len(n), function(j) mean(s$D[folds != folds[j]]), numeric(1))
      expect_equal(m$fits[[i]]$preds$phat, expected.phat)
    }
    external.se <- sqrt(sum((f$phi - mean(f$phi))^2)) / n
    for (combine in c("mean", "median")) {
      expect_equal(unname(coef(m, combine.method = combine)), f$theta)
      expect_equal(unname(se(m, combine.method = combine)), external.se, tolerance = 1e-12)
      expect_equal(unname(confidence_bounds(m, cf.y = 0, cf.d = 0,
                                            combine.method = combine)["att", ]),
                    f$theta + c(-1, 1) * qnorm(0.95) * external.se, tolerance = 1e-10)
    }
  }
})

test_that("repeated sensitivity inference matches independent endpoint IF calculations", {
  f <- repeated_adapter_fixture("did", reps = 4L)
  m <- f$model; n <- length(f$phi)
  k <- sqrt(0.08 * 0.06 / (1 - 0.06))
  ends <- lapply(m$results$main$treat, function(r) {
    e <- r$estimates
    pS <- e$nu2.s * r$psis$psi.sigma2.s + e$sigma2.s * r$psis$psi.nu2.s
    pm <- f$phi - k / (2 * sqrt(e$S2)) * pS
    pp <- f$phi + k / (2 * sqrt(e$S2)) * pS
    c(lwr = f$theta - k * sqrt(e$S2), upr = f$theta + k * sqrt(e$S2),
      se.lwr = sqrt(sum((pm - mean(pm))^2)) / n,
      se.upr = sqrt(sum((pp - mean(pp))^2)) / n)
  })
  ends <- do.call(rbind, ends)
  expect_gt(diff(range(ends[, "lwr"])), 1e-6)
  b <- dml_bounds(m, cf.y = 0.08, cf.d = 0.06)
  expect_length(b$results$main$treat, 4)
  for (i in 1:4) {
    expect_equal(b$results$main$treat[[i]]$estimates$se.theta.m, ends[i, "se.lwr"], ignore_attr = TRUE)
    expect_equal(b$results$main$treat[[i]]$estimates$se.theta.p, ends[i, "se.upr"], ignore_attr = TRUE)
  }
  for (combine in c("mean", "median")) {
    agg <- if (combine == "mean") mean else median
    lower <- agg(ends[, "lwr"]); upper <- agg(ends[, "upr"])
    expected <- c(lower - qnorm(0.95) * sqrt(agg(ends[, "se.lwr"]^2 + (ends[, "lwr"] - lower)^2)),
                  upper + qnorm(0.95) * sqrt(agg(ends[, "se.upr"]^2 + (ends[, "upr"] - upper)^2)))
    ci <- confidence_bounds(m, cf.y = 0.08, cf.d = 0.06, combine.method = combine, max = FALSE)
    expect_equal(unname(ci["att", ]), expected, tolerance = 1e-10)
    envelope <- confidence_bounds(m, cf.y = 0.08, cf.d = 0.06, combine.method = combine)
    expect_lte(envelope["att", "lwr"], ci["att", "lwr"] + 1e-10)
    expect_gte(envelope["att", "upr"], ci["att", "upr"] - 1e-10)
    expect_true(all(is.finite(robustness_value(m, combine.method = combine))))
    sens <- sensemakr(m, cf.y = 0.08, cf.d = 0.06, combine.method = combine)
    expect_s3_class(sens, "dml.sensemakr")
    expect_equal(sens$conf.bounds, envelope)
    expect_identical(sens$info$combine.method, combine)
    expect_equal(sens$bounds$theta.minus, lower)
    expect_equal(sens$bounds$theta.plus, upper)
    expect_equal(sens$sensitivity_stats[, "rv"],
                 unname(robustness_value(m, alpha = 1, combine.method = combine)))
    expect_equal(sens$sensitivity_stats[, "rva"],
                 unname(robustness_value(m, alpha = 0.05, combine.method = combine)))
  }
})

test_that("sensemakr benchmark tables and plots use the selected repetition aggregation", {
  f <- repeated_adapter_fixture("did", reps = 4L)
  for (combine in c("mean", "median")) {
    s <- sensemakr(f$model, benchmark_covariates = "x2", cf.y = 0.05, cf.d = 0.04,
                   combine.method = combine, alpha = 0.10)
    expected <- as.data.frame(benchmark_bounds(f$model, s$bench.bounds,
                              level = 0.90, combine.method = combine))
    gains <- summary(s$bench.bounds, combine.method = combine)$benchmarks
    expect_equal(s$bench.table$lwr, expected$lwr)
    expect_equal(s$bench.table$upr, expected$upr)
    expect_equal(s$bench.table$lwr.fixed, expected$lwr.fixed)
    expect_equal(s$bench.table$upr.fixed, expected$upr.fixed)
    expect_equal(s$bench.table$cf.y, unname(gains[, "gain.Y"]))
    expect_equal(s$bench.table$cf.d, unname(gains[, "gain.D"]))
  }
  expect_error(sensemakr(f$model, combine.method = "unknown"), "arg.*one of")
  testthat::local_mocked_bindings(ovb_contour_plot = function(...) list(...),
                                 .package = "dml.sensemakr")
  s$info$combine.method <- "mean"
  expect_identical(plot(s)$combine.method, "mean")
  expect_identical(plot(s, combine.method = "median")$combine.method, "median")
  s$info$combine.method <- NULL
  expect_identical(plot(s)$combine.method, "median")
})

test_that("repeated benchmarks preserve paired folds, IFs and source uncertainty", {
  for (method in c("did", "trad", "imp")) {
    f <- repeated_adapter_fixture(method)
    path <- tempfile(fileext = ".rds"); saveRDS(f$model, path)
    m <- readRDS(path); unlink(path)
    m$external_fit <- NULL; m$call <- quote(stop("original call unavailable"))
    s <- m$adapter_data; X <- s$X[, c("(Intercept)", "x1")]
    fun <- if (method == "did") dml.sensemakr:::.did_benchmark_refit else
      dml.sensemakr:::.drdid_benchmark_refit
    refit <- fun(m, X)
    expect_identical(refit$adapter_info$fold_id, m$adapter_info$fold_id)
    expect_identical(refit$adapter_info$rep_seeds, m$adapter_info$rep_seeds)
    expect_identical(refit$adapter_info$unit_id, m$adapter_info$unit_id)
    source.fun <- if (method == "imp") DRDID::drdid_imp_panel else DRDID::drdid_panel
    independent <- source.fun(s$y1, s$y0, s$D, X, inffunc = TRUE)
    pure <- repeated_direct_dml(m, X)
    bench <- dml_benchmark(m, "x2")
    expect_equal(nrow(bench$benchmarks$x2), 3)
    expect_equal(bench$benchmarks$x2$theta.sj, rep(as.numeric(independent$ATT), 3))
    for (i in 1:3) {
      expect_equal(bench$benchmarks_psis$x2$psi.theta.s.wo[[i]], as.numeric(independent$att.inf.func))
      expect_equal(bench$benchmarks_psis$x2$psi.sigma2.s.wo[[i]], pure$results$main$treat[[i]]$psis$psi.sigma2.s)
      expect_equal(bench$benchmarks_psis$x2$psi.nu2.s.wo[[i]], pure$results$main$treat[[i]]$psis$psi.nu2.s)
    }
    diff.phi <- f$phi - as.numeric(independent$att.inf.func)
    diff.se <- sqrt(mean(diff.phi^2) / length(diff.phi))
    for (combine in c("mean", "median")) {
      tab <- summary(bench, combine.method = combine)$benchmarks
      expect_equal(unname(tab["x2", "delta"]), f$theta - as.numeric(independent$ATT))
      expect_equal(unname(tab["x2", "se.delta"]), diff.se)
      expect_s3_class(benchmark_bounds(m, bench, combine.method = combine), "dml_benchmark_bounds")
    }
    # A changed refit seed must never silently pair IFs from different folds.
    m$refit_spec$dml_args$cf.seed <- 18
    expect_error(fun(m, X), "cross-fitting partitions differ")
  }
})

test_that("batch conversion and callbacks retain every cell's repetitions", {
  f <- repeated_adapter_fixture("did")
  batch <- dml_from_did_cells(f$source,
    cells = data.frame(group = c(2, 2), time = c(2, 3)), dml_args = f$args)
  expect_identical(batch$status$status, rep("ok", 2))
  expect_equal(vapply(batch$models, function(m) length(m$results$main$treat), integer(1)),
                c(`g=2,t=2` = 3L, `g=2,t=3` = 3L))
  result <- did_cell_apply(batch, confidence_bounds, cf.y = 0.05, cf.d = 0.05,
                            combine.method = "mean")
  expect_identical(result$status$status, rep("ok", 2))
  expect_equal(result$values[[1]], confidence_bounds(batch$models[[1]],
    cf.y = 0.05, cf.d = 0.05, combine.method = "mean"))
  bench <- did_cell_apply(batch, dml_benchmark, benchmark_covariates = "x2")
  expect_identical(bench$status$status, rep("ok", 2))
  expect_equal(vapply(bench$values, function(b) nrow(b$benchmarks$x2), integer(1)),
                c(`g=2,t=2` = 3L, `g=2,t=3` = 3L))
})

test_that("later invalid repetitions cannot be silently discarded", {
  f <- repeated_adapter_fixture("trad")
  s <- f$model$adapter_data
  raw <- repeated_direct_dml(f$model)
  run <- function() dml.sensemakr:::.panel_hybrid_aux(s, f$theta, f$phi, f$args)
  # Exercise failures in the final, rather than the first, backend result.
  backend <- raw
  fake.dml <- function(...) backend
  formals(fake.dml) <- formals(dml)
  testthat::local_mocked_bindings(dml = fake.dml, .package = "dml.sensemakr")
  backend$results$main$treat[[3]]$trim.summary$trimmed_num$all <- 1L
  expect_error(run(), "Repetition 3: DML propensity clipping")
  backend <- raw; backend$fits[[3]]$preds$dhat[1] <- NA_real_
  expect_error(run(), "Repetition 3: DML propensity clipping")
  backend <- raw; backend$results$main$treat[[3]]$estimates$S2 <- NA_real_
  expect_error(run(), "Repetition 3: DML produced")
  backend <- raw; backend$results$main$treat[[3]]$psis$psi.nu2.s[1] <- NA_real_
  expect_error(run(), "aligned finite unit-level vectors")
  backend <- raw; backend$results$main$treat <- backend$results$main$treat[1:2]
  expect_error(run(), "every requested")
})

test_that("all partitions are preflighted before any nuisance fit", {
  skip_if_not_installed("DRDID")
  # Two treated units: they must fall in different folds for each training set
  # to contain a treated observation. Find a seed valid only in repetition 1.
  D <- c(1, 1, rep(0, 8)); n <- length(D)
  candidate <- which(vapply(1:100, function(seed) {
    set.seed(seed); seeds <- sample.int(1e6, 2)
    good <- vapply(seeds, function(s) {
      set.seed(s); f <- rep(1:2, length.out = n)[sample.int(n)]
      f[1] != f[2]
    }, logical(1))
    good[1] && !good[2]
  }, logical(1)))[1]
  expect_false(is.na(candidate))
  y0 <- seq_len(n); y1 <- y0 + seq_len(n)^2 / 10 + D
  source <- DRDID::drdid_panel(y1, y0, D, covariates = NULL, inffunc = TRUE)
  args <- list(yreg = "lm", dreg = "glm", d.class = TRUE,
                cf.folds = 2, cf.seed = candidate, ps.trim = 0)
  one <- dml_from_drdid(source, y1 = y1, y0 = y0, D = D, dml_args = args)
  expect_length(one$results$main$treat, 1)
  args$cf.reps <- 2
  fake.dml <- function(...) stop("nuisance fitting started")
  formals(fake.dml) <- formals(dml)
  testthat::local_mocked_bindings(dml = fake.dml,
                                 .package = "dml.sensemakr")
  expect_error(dml_from_drdid(source, y1 = y1, y0 = y0, D = D, dml_args = args),
                "Repetition 2, training fold")
})

test_that("one repetition stays compatible and default learners support repeats", {
  f <- repeated_adapter_fixture("did", reps = 1L)
  args <- f$args; args$cf.reps <- NULL
  default <- dml_from_did(f$source, 2, 2, dml_args = args)
  expect_identical(default$results, f$model$results)
  expect_identical(default$adapter_info$fold_id, f$model$adapter_info$fold_id)
  expect_null(dim(default$adapter_info$fold_id))
  expect_length(default$adapter_info$dml.att, 1)
  skip_if_not_installed("ranger")
  repeated <- dml_from_did(f$source, 2, 2, dml_args = list(cf.reps = 2))
  expect_length(repeated$results$main$treat, 2)
  expect_true(all(vapply(repeated$results$main$treat,
                         function(r) is.finite(r$estimates$S2) && r$estimates$S2 > 0,
                         logical(1))))
  s <- f$model$adapter_data
  source <- DRDID::drdid_imp_panel(s$y1, s$y0, s$D, covariates = NULL, inffunc = TRUE)
  constant <- dml_from_drdid(source, y1 = s$y1, y0 = s$y0, D = s$D,
                             dml_args = list(cf.reps = 2))
  expect_length(constant$results$main$treat, 2)
  expect_identical(constant$info$dreg$preProcess, "ignore")
})
