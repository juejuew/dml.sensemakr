# Deterministic, small parametric fits exercise the actual external estimator.
adapter_fixture <- function(method = "imp") {
  skip_if_not_installed("DRDID")
  set.seed(421)
  n <- 300
  X <- cbind(`(Intercept)` = 1, x1 = runif(n, -1, 1), x2 = runif(n, -1, 1))
  D <- rbinom(n, 1, plogis(0.1 + 0.5 * X[, "x1"] + 0.6 * X[, "x2"]))
  y0 <- X[, "x1"] + rnorm(n)
  y1 <- y0 + 0.7 * D + X[, "x1"] + 1.5 * X[, "x2"] + rnorm(n)
  fun <- if (method == "imp") DRDID::drdid_imp_panel else DRDID::drdid_panel
  source <- fun(y1, y0, D, X, inffunc = TRUE)
  args <- list(yreg = "lm", dreg = "glm", d.class = TRUE,
               cf.folds = 3, cf.seed = 17, ps.trim = 0)
  input <- list(y1 = y1, y0 = y0, D = D, covariates = X, dml_args = args)
  list(source = source, input = input,
       model = do.call(dml_from_drdid, c(list(fit = source), input)))
}

test_that("both panel methods preserve the external ATT and full unscaled IF", {
  for (method in c("imp", "trad")) {
    f <- adapter_fixture(method)
    m <- f$model; result <- m$results$main$treat[[1]]
    expect_s3_class(m, "dml_drdid")
    expect_s3_class(m, "dml")
    expect_identical(result$estimates$theta.s, as.numeric(f$source$ATT))
    expect_identical(result$psis$psi.theta.s, as.numeric(f$source$att.inf.func))
    expect_equal(unname(coef(m)), as.numeric(f$source$ATT))
    expect_equal(unname(se(m)), as.numeric(f$source$se), tolerance = 1e-8)
    expect_equal(result$estimates$theta.1s - result$estimates$theta.0s,
                 f$source$ATT)
    expect_true(m$adapter_info$diagnostics$reproduced)
    expect_length(m$adapter_info$diagnostics$propensity, length(f$input$D))
    expect_output(print(m), "Hybrid panel DiD")
    expect_error(dml_gate(m, groups = list(all = rep(TRUE, length(f$input$D)))),
                 "external")
  }
})

test_that("joint IF covariance agrees with the delta method and downstream bounds", {
  f <- adapter_fixture(); m <- f$model
  r <- m$results$main$treat[[1]]; e <- r$estimates
  phi <- cbind(as.numeric(f$source$att.inf.func), r$psis$psi.sigma2.s,
               r$psis$psi.nu2.s)
  n <- nrow(phi)
  V <- cov(phi) * (n - 1) / n^2
  expect_equal(unname(r$joint.covariance), unname(V), tolerance = 1e-12)
  k <- sqrt(0.08 * 0.06 / (1 - 0.06))
  gradient <- c(1, -k / 2 * sqrt(e$nu2.s / e$sigma2.s),
                 -k / 2 * sqrt(e$sigma2.s / e$nu2.s))
  expected.se <- sqrt(drop(t(gradient) %*% V %*% gradient))
  b <- dml_bounds(m, cf.y = 0.08, cf.d = 0.06)
  expect_s3_class(plot(m), "ggplot")
  expect_s3_class(plot(b), "ggplot")
  expect_equal(b$results$main$treat[[1]]$estimates$se.theta.m,
               expected.se, tolerance = 1e-8)
  ci <- confidence_bounds(m, cf.y = 0.08, cf.d = 0.06, max = FALSE)
  expect_equal(unname(ci["att", "lwr"]),
               e$theta.s - k * sqrt(e$S2) - qnorm(0.95) * expected.se,
               tolerance = 1e-8)
  zero <- confidence_bounds(m, cf.y = 0, cf.d = 0.06)
  expect_equal(unname(zero["att", ]),
               e$theta.s + c(-1, 1) * qnorm(0.95) * f$source$se,
               tolerance = 1e-8)
  envelope <- confidence_bounds(m, cf.y = 0.08, cf.d = 0.06)
  expect_lte(envelope["att", "lwr"], ci["att", "lwr"] + 1e-10)
  expect_gte(envelope["att", "upr"], ci["att", "upr"] - 1e-10)
  rv <- robustness_value(m)
  expect_true(all(is.finite(rv) & rv >= 0 & rv <= 1))
  expect_s3_class(sensemakr(m, cf.y = 0.08, cf.d = 0.06), "dml.sensemakr")
})

test_that("long input uses DRDID's canonical unit order and records original rows", {
  f <- adapter_fixture(); v <- f$input; n <- length(v$D)
  long <- data.frame(id = rep(seq_len(n), 2), time = rep(0:1, each = n),
                      D = rep(v$D, 2), y = c(v$y0, v$y1),
                      x1 = rep(v$covariates[, "x1"], 2),
                      x2 = rep(v$covariates[, "x2"], 2))
  source <- DRDID::drdid(yname = "y", tname = "time", idname = "id",
                         dname = "D", xformla = ~ x1 + x2, data = long,
                         panel = TRUE, inffunc = TRUE)
  shuffled <- long[sample.int(nrow(long)), ]
  m <- dml_from_drdid(source, data = shuffled, dml_args = v$dml_args)
  expect_equal(as.numeric(coef(m)), as.numeric(source$ATT))
  expect_equal(m$adapter_info$unit_id, seq_len(n))
  rows <- m$adapter_data$input.info$source_rows
  expect_equal(shuffled$id[rows$pre], seq_len(n))
  expect_equal(shuffled$id[rows$post], seq_len(n))
  expect_equal(shuffled$time[rows$pre], rep(0, n))
  expect_equal(shuffled$time[rows$post], rep(1, n))
  expect_equal(m$adapter_info$fold_id, f$model$adapter_info$fold_id)
  expect_error(dml_from_drdid(source, data = long[-1, ], dml_args = v$dml_args),
                 "balanced")
  long$y[1] <- NA_real_
  expect_error(dml_from_drdid(source, data = long, dml_args = v$dml_args),
                 "Missing")
})

test_that("unsafe or incompatible inputs are rejected explicitly", {
  f <- adapter_fixture()
  run <- function(source = f$source, input = f$input)
    do.call(dml_from_drdid, c(list(fit = source), input))
  bad <- f$source; bad$att.inf.func <- NULL
  expect_error(run(bad), "inffunc = TRUE")
  bad <- f$source; bad$argu$panel <- FALSE
  expect_error(run(bad), "repeated cross-sections")
  bad <- f$source; bad$argu$type <- "reg"
  expect_error(run(bad), "doubly robust")
  bad <- f$source; bad$ps.flag <- 2
  expect_error(run(bad), "fallback")
  input <- f$input; input$i.weights <- seq_along(input$D)
  expect_error(run(input = input), "constant positive")
  input <- f$input; input$y1 <- rev(input$y1)
  expect_error(run(input = input), "ordered IF")
  bad <- f$source; bad$att.inf.func <- rev(bad$att.inf.func)
  expect_error(run(bad), "ordered IF")
  input <- f$input; input$dml_args$cf.reps <- 0
  expect_error(run(input = input), "cf.reps must be an integer")
  input <- f$input; input$dml_args$ps.trim <- 0.49
  expect_error(run(input = input), "DML propensity clipping")
  input <- f$input; input$dml_args$cf.seed <- 1.5
  expect_error(run(input = input), "integer cf.seed")
  v <- f$input
  trimmed <- DRDID::drdid_imp_panel(v$y1, v$y0, v$D, v$covariates,
                                    trim.level = 0.05, inffunc = TRUE)
  expect_error(run(trimmed), "Source propensity trimming")
})

test_that("benchmarks refit both estimators after serialization", {
  f <- adapter_fixture()
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path))
  saveRDS(f$model, path)
  m <- readRDS(path)
  # The adapter refit must not evaluate its original call or use its environment.
  m$call <- quote(stop("unavailable original call"))
  m$call.env <- NULL
  bench <- dml_benchmark(m, "x2")
  v <- f$input
  reduced <- DRDID::drdid_imp_panel(v$y1, v$y0, v$D,
                                    v$covariates[, c("(Intercept)", "x1")],
                                    inffunc = TRUE)
  expect_equal(bench$benchmarks$x2$theta.sj, as.numeric(reduced$ATT))
  expect_equal(bench$benchmarks_psis$x2$psi.theta.s.wo[[1]],
               as.numeric(reduced$att.inf.func))
  expect_equal(bench$benchmarks_psis$x2$psi.theta.s[[1]],
               as.numeric(f$source$att.inf.func))
  refit <- dml.sensemakr:::.drdid_benchmark_refit(m,
              v$covariates[, c("(Intercept)", "x1"), drop = FALSE])
  expect_identical(refit$adapter_info$unit_id, m$adapter_info$unit_id)
  expect_identical(refit$adapter_info$fold_id, m$adapter_info$fold_id)
  expect_equal(bench$benchmarks_psis$x2$psi.sigma2.s.wo[[1]],
               refit$results$main$treat[[1]]$psis$psi.sigma2.s)
})

test_that("source bootstrap SEs do not replace joint analytic IF inference", {
  f <- adapter_fixture(); v <- f$input
  source <- DRDID::drdid_imp_panel(v$y1, v$y0, v$D, v$covariates,
                                   boot = TRUE, boot.type = "multiplier",
                                   nboot = 99, inffunc = TRUE)
  m <- do.call(dml_from_drdid, c(list(fit = source), v))
  expect_identical(m$adapter_info$original.se, source$se)
  expect_true(m$adapter_info$original.boot)
  expect_equal(as.numeric(se(m)), as.numeric(f$source$se), tolerance = 1e-8)
  expect_identical(m$results$main$treat[[1]]$psis$psi.theta.s,
                   as.numeric(source$att.inf.func))
})

test_that("intercept-only sources and reduced models are supported", {
  f <- adapter_fixture(); v <- f$input
  v$covariates <- NULL
  source <- DRDID::drdid_imp_panel(v$y1, v$y0, v$D, covariates = NULL,
                                   inffunc = TRUE)
  m <- do.call(dml_from_drdid, c(list(fit = source), v))
  expect_equal(as.numeric(coef(m)), as.numeric(source$ATT))
  expect_true(is.finite(m$results$main$treat[[1]]$estimates$S2))
  expect_error(dml_benchmark(m, "(Intercept)"), "retain at least one")
})

test_that("sensemakr uses the requested alpha in hybrid confidence tables", {
  f <- adapter_fixture(); m <- f$model
  s <- sensemakr(m, cf.y = 0.08, cf.d = 0.06,
                 benchmark_covariates = "x2", alpha = 0.1)
  ci <- confidence_bounds(m, cf.y = 0.08, cf.d = 0.06, level = 0.9)
  expect_equal(s$bounds$lwr, unname(ci["att", "lwr"]))
  expect_equal(s$bounds$upr, unname(ci["att", "upr"]))
  b <- as.data.frame(benchmark_bounds(m, s$bench.bounds, level = 0.9))
  expect_equal(s$bench.table$lwr, b$lwr)
  expect_equal(s$bench.table$upr, b$upr)
})

test_that("default DML learners work with an external panel fit", {
  skip_if_not_installed("ranger")
  f <- adapter_fixture(); v <- f$input; v$dml_args <- list()
  m <- do.call(dml_from_drdid, c(list(fit = f$source), v))
  expect_equal(as.numeric(coef(m)), as.numeric(f$source$ATT))
  expect_equal(m$info$cf.folds, 5L)
  expect_true(is.finite(m$results$main$treat[[1]]$estimates$S2))
})

test_that("benchmark labels cannot silently overwrite another group", {
  f <- adapter_fixture()
  for (groups in list(list(same = "x1", same = "x2"),
                      list(x1 = "x2", "x1"), c("x1", "x1"))) {
    expect_error(dml_benchmark(f$model, groups), "group names must be unique")
  }
})

test_that("raw long weights are checked before DRDID normalizes them", {
  f <- adapter_fixture(); v <- f$input; n <- length(v$D)
  long <- data.frame(id = rep(seq_len(n), 2), time = rep(0:1, each = n),
                    D = rep(v$D, 2), y = c(v$y0, v$y1),
                    x1 = rep(v$covariates[, "x1"], 2),
                    x2 = rep(v$covariates[, "x2"], 2), w = -2)
  # DRDID itself accepts these weights after dividing them by their mean.
  source <- DRDID::drdid(yname = "y", tname = "time", idname = "id",
                        dname = "D", xformla = ~ x1 + x2, data = long,
                        weightsname = "w", panel = TRUE, inffunc = TRUE)
  run <- function(data = long, weightsname = "w")
    dml_from_drdid(source, data = data, weightsname = weightsname,
                  dml_args = v$dml_args)
  expect_error(run(), "constant positive")
  for (bad in list(rep(0, 2 * n), rep(Inf, 2 * n), rep(NA_real_, 2 * n),
                   rep("2", 2 * n), rep(c(1, 2), each = n))) {
    long$w <- bad
    expect_error(run(long), "constant positive")
  }
  expect_error(run(weightsname = c("w", "x1")), "Supply weightsname explicitly")
  # A constant positive scale must produce exactly the unweighted analysis.
  long$w <- 2.5
  m <- run(long)
  expect_equal(m$adapter_data$i.weights, rep(1, n))
  expect_equal(m$results$main, f$model$results$main)
  for (bad in list(rep(-2, n), rep(0, n), rep(Inf, n), rep(NA_real_, n),
                   rep("2", n), seq_len(n))) {
    expect_error(do.call(dml_from_drdid,
                        c(list(fit = f$source, i.weights = bad), v)),
                 "constant positive")
  }
})

test_that("intercept-only learners resolve shared defaults and explicit overrides", {
  skip_if_not_installed("ranger")
  f <- adapter_fixture(); v <- f$input; v$covariates <- NULL
  source <- DRDID::drdid_imp_panel(v$y1, v$y0, v$D, covariates = NULL,
                                  inffunc = TRUE)
  cases <- list(
    list(args = list(), y = "ranger", d = "ranger"),
    list(args = list(reg = "lm"), y = "lm", d = "lm"),
    list(args = list(reg = "glm", yreg = "lm"), y = "lm", d = "glm"),
    list(args = list(reg = "lm", dreg = "glm", d.class = TRUE),
         y = "lm", d = "glm"),
    list(args = list(reg = list(method = "glm", preProcess = "ignore"),
                     yreg = list(yreg0 = "lm", yreg1 = "lm"), d.class = TRUE),
         y = "lm", d = "glm")
  )
  for (case in cases) {
    v$dml_args <- case$args
    expect_warning(m <- do.call(dml_from_drdid, c(list(fit = source), v)), NA)
    expect_identical(attr(m$info$yreg$yreg0$method, "name"), case$y)
    expect_identical(attr(m$info$dreg$method, "name"), case$d)
    expect_identical(m$info$yreg$yreg0$preProcess, "ignore")
    expect_identical(m$info$dreg$preProcess, "ignore")
    expect_equal(as.numeric(coef(m)), as.numeric(source$ATT))
    expect_true(is.finite(m$results$main$treat[[1]]$estimates$S2))
  }
})

test_that("long and wide inputs agree with factors and interactions for both methods", {
  f <- adapter_fixture(); v <- f$input; n <- length(v$D)
  units <- data.frame(id = 100 + 7 * seq_len(n),
                      x1 = v$covariates[, "x1"], x2 = v$covariates[, "x2"])
  units$region <- cut(units$x2, c(-Inf, -1/3, 1/3, Inf),
                      labels = c("west", "central", "east"))
  X <- model.matrix(~ x1 * region + x2, units)
  long <- units[rep(seq_len(n), 2), ]
  long$time <- rep(c(2000, 2005), each = n)
  long$D <- rep(v$D, 2); long$y <- c(v$y0, v$y1); long$w <- 2.5
  long <- long[sample.int(2 * n), ]
  for (method in c("imp", "trad")) {
    source.long <- DRDID::drdid(yname = "y", tname = "time", idname = "id",
      dname = "D", xformla = ~ x1 * region + x2, data = long,
      weightsname = "w", panel = TRUE, estMethod = method, inffunc = TRUE)
    fun <- if (method == "imp") DRDID::drdid_imp_panel else DRDID::drdid_panel
    source.wide <- fun(v$y1, v$y0, v$D, X, i.weights = rep(2.5, n), inffunc = TRUE)
    ml <- dml_from_drdid(source.long, data = long, dml_args = v$dml_args)
    mw <- dml_from_drdid(source.wide, y1 = v$y1, y0 = v$y0, D = v$D,
      covariates = X, unit_id = units$id, i.weights = rep(2.5, n), dml_args = v$dml_args)
    expect_equal(as.numeric(source.long$ATT), as.numeric(source.wide$ATT))
    expect_equal(as.numeric(source.long$att.inf.func), as.numeric(source.wide$att.inf.func))
    expect_identical(ml$adapter_info$unit_id, units$id)
    expect_identical(ml$adapter_info$fold_id, mw$adapter_info$fold_id)
    expect_identical(colnames(ml$adapter_data$X), colnames(X))
    expect_equal(c(ml$adapter_data$X), c(X))
    for (key in c("y0", "y1", "D", "i.weights"))
      expect_equal(ml$adapter_data[[key]], mw$adapter_data[[key]])
    expect_equal(ml$results$main, mw$results$main, tolerance = 1e-10)
    expect_equal(ml$fits[[1]]$preds, mw$fits[[1]]$preds, tolerance = 1e-10)
    expect_equal(confidence_bounds(ml, cf.y = .08, cf.d = .06),
                 confidence_bounds(mw, cf.y = .08, cf.d = .06))
  }
})

test_that("multiple grouped benchmarks and learner overrides match independent fits", {
  skip_if_not_installed("ranger")
  groups <- list(first = "x1", both = c("x1", "x2"))
  ylearner <- list(method = "ranger", num.trees = 100, num.threads = 1,
    preProcess = "ignore", trControl = list(method = "none", allowParallel = FALSE),
    tuneGrid = data.frame(mtry = 1, splitrule = "variance", min.node.size = 10))
  dlearner <- list(method = "glm", family = binomial(link = "probit"),
                   preProcess = "ignore")
  for (method in c("imp", "trad")) {
    f <- adapter_fixture(method); m <- f$model; v <- f$input; n <- length(v$D)
    bench <- dml_benchmark(m, groups, yreg = ylearner, dreg = dlearner)
    expect_identical(names(bench$benchmarks), names(groups))
    expect_identical(names(bench$benchmarks_psis), names(groups))
    tab <- summary(bench)$benchmarks
    expect_identical(rownames(tab), names(groups))
    full <- m$results$main$treat[[1]]
    for (label in names(groups)) {
      X <- v$covariates[, !colnames(v$covariates) %in% groups[[label]], drop = FALSE]
      fun <- if (method == "imp") DRDID::drdid_imp_panel else DRDID::drdid_panel
      source <- fun(v$y1, v$y0, v$D, X, trim.level = f$source$argu$trim.level,
                     inffunc = TRUE)
      # Independent public DML fit: do not use the adapter's refit helper.
      auxiliary <- dml(v$y1 - v$y0, v$D, X, model = "npm", target = "att",
        conditional = TRUE, centered = FALSE, cf.reps = 1, dirty.tuning = FALSE,
        yreg = ylearner, dreg = dlearner, d.class = TRUE,
        cf.folds = 3, cf.seed = 17, ps.trim = 0, verbose = FALSE)
      reduced <- auxiliary$results$main$treat[[1]]
      b <- bench$benchmarks[[label]]; p <- bench$benchmarks_psis[[label]]
      expect_equal(b$theta.sj, as.numeric(source$ATT))
      expect_equal(p$psi.theta.s.wo[[1]], as.numeric(source$att.inf.func))
      expect_equal(p$psi.sigma2.s.wo[[1]], reduced$psis$psi.sigma2.s)
      expect_equal(p$psi.nu2.s.wo[[1]], reduced$psis$psi.nu2.s)
      expect_equal(p$psi.theta.s[[1]], as.numeric(f$source$att.inf.func))
      expect_equal(p$psi.sigma2.s[[1]], full$psis$psi.sigma2.s)
      expect_equal(p$psi.nu2.s[[1]], full$psis$psi.nu2.s)
      sf <- full$estimates$sigma2.s; sr <- reduced$estimates$sigma2.s
      nf <- full$estimates$nu2.s; nr <- reduced$estimates$nu2.s
      expect_equal(b$gain.Y, max(0, sr / sf - 1))
      expect_equal(b$gain.D, max(0, nf / nr - 1))
      phi <- cbind(f$source$att.inf.func, full$psis$psi.sigma2.s,
        full$psis$psi.nu2.s, source$att.inf.func,
        reduced$psis$psi.sigma2.s, reduced$psis$psi.nu2.s)
      gradients <- cbind(delta = c(1, 0, 0, -1, 0, 0),
        gain.Y = c(0, -sr / sf^2, 0, 0, 1 / sf, 0),
        gain.D = c(0, 0, 1 / nr, 0, 0, -nf / nr^2))
      joint <- cov(phi) * (n - 1) / n^2
      expected.se <- sqrt(diag(t(gradients) %*% joint %*% gradients))
      expect_equal(unname(tab[label, paste0("se.", names(expected.se))]),
                   unname(expected.se), tolerance = 1e-8)
      expect_equal(p$psi.GY[[1]], drop(phi %*% gradients[, "gain.Y"]))
      expect_equal(p$psi.GD[[1]], drop(phi %*% gradients[, "gain.D"]))
    }
    expect_true(all(is.finite(as.matrix(tab))))
    expect_true(all(is.finite(as.matrix(as.data.frame(benchmark_bounds(m, bench))))))
    # Benchmark refits must not mutate the full model's learner settings.
    expect_identical(m$refit_spec$dml_args$yreg, "lm")
    expect_identical(m$refit_spec$dml_args$dreg, "glm")
  }
})
