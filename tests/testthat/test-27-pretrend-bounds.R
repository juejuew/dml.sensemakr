pretrend_fixture <- function(n = 90L, method = "parametric") {
  skip_if_not_installed("DRDID")
  set.seed(10626)
  x <- runif(n, -1, 1); z <- runif(n, -1, 1)
  D <- rbinom(n, 1, plogis(0.2 + 0.35 * x - 0.2 * z))
  Xpre <- cbind(`(Intercept)` = 1, x = x, z = z)
  Xpost <- cbind(`(Intercept)` = 1, x = x + 0.25 * z^2, z = z)
  eps <- rnorm(n); y0 <- rnorm(n)
  y1 <- y0 + 0.7 * D + 0.4 * x - 0.2 * z + eps
  y2 <- y1 + 1.1 * D + 0.6 * Xpost[, "x"] + 0.2 * z +
    0.4 * eps + 0.5 * rnorm(n)
  input <- list(y0 = y0, y1 = y1, y2 = y2, D = D,
                Xpre = Xpre, Xpost = Xpost, id = seq_len(n) * 7 + 10)
  args <- if (method == "parametric") list() else list(yreg = "lm", dreg = "glm",
    d.class = TRUE, cf.folds = 3, cf.reps = 2, cf.seed = 23, ps.trim = 0)
  make <- function(pre = FALSE, ix = seq_len(n), sign = 1) {
    lo <- if (pre) input$y0 else input$y1
    hi <- if (pre) input$y1 else input$y2
    X <- if (pre) Xpre else Xpost
    if (pre) hi <- lo + sign * (hi - lo)
    source <- DRDID::drdid_panel(hi[ix], lo[ix], D[ix], X[ix, , drop = FALSE], inffunc = TRUE)
    dml_from_drdid(source, y1 = hi[ix], y0 = lo[ix], D = D[ix],
      covariates = X[ix, , drop = FALSE], unit_id = input$id[ix],
      method = method, dml_args = args)
  }
  list(input = input, post = make(), pre = make(TRUE), make = make)
}

# Independent base-R oracle: no adapter moment, bounds or IF helper is used.
pretrend_mass_oracle <- function(input, w, k = 1.3) {
  fit <- function(y, X) {
    d <- input$D
    ps <- glm.fit(X, d, weights = w, family = quasibinomial(),
      control = glm.control(epsilon = 1e-13, maxit = 100))
    beta <- lm.wfit(X[d == 0, , drop = FALSE], y[d == 0], w[d == 0])$coefficients
    r <- y - as.numeric(X %*% beta)
    odds <- ps$fitted.values / (1 - ps$fitted.values)
    theta <- weighted.mean(r[d == 1], w[d == 1]) -
      weighted.mean(r[d == 0], (w * odds)[d == 0])
    p <- weighted.mean(d, w)
    h <- odds * (1 - p) / p
    a <- weighted.mean(r[d == 0]^2, w[d == 0])
    b <- 2 * weighted.mean(h[d == 1], w[d == 1]) -
      weighted.mean(h[d == 0]^2, w[d == 0])
    c(theta = theta, S = sqrt(a * b))
  }
  post <- fit(input$y2 - input$y1, input$Xpost)
  pre <- fit(input$y1 - input$y0, input$Xpre)
  plain <- k * abs(pre["theta"])
  scaled <- plain * post["S"] / pre["S"]
  unname(c(post["theta"] - plain, post["theta"] + plain,
           post["theta"] - scaled, post["theta"] + scaled))
}

test_that("both pre-trend formulas and their joint IFs match independent refitting", {
  f <- pretrend_fixture(72L)
  n <- length(f$input$D); w <- rep(1, n); step <- 1e-4
  a <- pretrend_bounds(f$post, f$pre, k = 1.3, type = "bias")
  b <- pretrend_bounds(f$post, f$pre, k = 1.3, type = "scaled")
  expected <- pretrend_mass_oracle(f$input, w)
  observed <- c(a$bounds$theta.minus, a$bounds$theta.plus,
                b$bounds$theta.minus, b$bounds$theta.plus)
  expect_equal(observed, expected, tolerance = 1e-7)
  numerical <- t(vapply(seq_len(n), function(i) {
    wp <- wm <- w; wp[i] <- wp[i] + step; wm[i] <- wm[i] - step
    n * (pretrend_mass_oracle(f$input, wp) - pretrend_mass_oracle(f$input, wm)) / (2 * step)
  }, numeric(4)))
  aa <- a$results[[1L]][[1L]]$psis; bb <- b$results[[1L]][[1L]]$psis
  analytic <- cbind(aa$psi.theta.minus, aa$psi.theta.plus,
                    bb$psi.theta.minus, bb$psi.theta.plus)
  cat("\nPre-trend IF maximum absolute errors (bias lower/upper, scaled lower/upper):\n")
  print(apply(abs(analytic - numerical), 2L, max), digits = 5)
  expect_equal(unname(analytic), unname(numerical), tolerance = 3e-5)
  expected.se <- apply(numerical, 2L, function(x) sqrt(sum((x - mean(x))^2)) / n)
  expect_equal(c(a$bounds$se.minus, a$bounds$se.plus,
                 b$bounds$se.minus, b$bounds$se.plus), expected.se, tolerance = 1e-7)
  V <- a$joint.covariance[[1L]]
  grad <- c(1, 0, 0, -1.3 * sign(a$bounds$theta.pre), 0, 0)
  expect_equal(a$bounds$se.minus^2, unname(t(grad) %*% V %*% grad)[1L], tolerance = 1e-12)
  expect_gt(abs(a$bounds$se.minus^2 - (V[1, 1] + 1.3^2 * V[4, 4])), 1e-5)
})

test_that("pre-trend sign, zero multipliers, row ordering and serialization work", {
  f <- pretrend_fixture()
  a <- pretrend_bounds(f$post, f$pre, k = c(0, 0.5, 2), level = 0.975)
  neg <- pretrend_bounds(f$post, f$make(TRUE, sign = -1), k = c(0, 0.5, 2), level = 0.975)
  expect_equal(a$bounds$theta.minus, neg$bounds$theta.minus, tolerance = 1e-10)
  expect_equal(a$bounds$theta.plus, neg$bounds$theta.plus, tolerance = 1e-10)
  expect_equal(a$bounds$se.minus, neg$bounds$se.minus, tolerance = 1e-10)
  expect_equal(a$bounds$se.plus, neg$bounds$se.plus, tolerance = 1e-10)
  r <- f$post$results$main$treat[[1L]]
  expect_equal(a$bounds$theta.minus[1L], r$estimates$theta.s)
  expect_equal(a$bounds$se.minus[1L], r$estimates$se.theta.s, tolerance = 1e-12)
  expect_equal(c(a$bounds$lwr[1L], a$bounds$upr[1L]),
    r$estimates$theta.s + c(-1, 1) * qnorm(0.975) * r$estimates$se.theta.s)
  permuted <- f$make(TRUE, ix = rev(seq_along(f$input$D)))
  direct <- pretrend_bounds(f$post, f$pre, type = "scaled")
  shuffled <- pretrend_bounds(f$post, permuted, type = "scaled")
  expect_equal(direct$bounds, shuffled$bounds, tolerance = 1e-8)
  expect_equal(direct$joint.covariance, shuffled$joint.covariance, tolerance = 1e-8)
  restored <- unserialize(serialize(f$post, NULL)); restored$external_fit <- NULL
  expect_equal(pretrend_bounds(restored, f$pre, type = "scaled")$bounds, direct$bounds)
  expect_output(print(direct), "uncertainty propagated")
})

test_that("undefined delta inference is explicit and invalid pairs are rejected", {
  f <- pretrend_fixture()
  zero <- f$pre; zero$results$main$treat[[1L]]$estimates$theta.s <- 0
  expect_warning(z <- pretrend_bounds(f$post, zero, k = c(0, 1), type = "scaled"),
                 "non-differentiable")
  expect_identical(z$bounds$status, c("ok", "zero_pretrend"))
  expect_true(is.finite(z$bounds$se.minus[1L]))
  expect_true(is.na(z$bounds$se.minus[2L]))
  expect_true(is.na(z$bounds$lwr[2L]))
  expect_equal(z$bounds$theta.minus, rep(unname(coef(f$post)), 2))
  missing.ids <- f$pre; missing.ids$adapter_data$id.supplied <- FALSE
  expect_error(pretrend_bounds(f$post, missing.ids), "explicit unique unit IDs")
  other.ids <- f$pre; other.ids$adapter_data$unit_id[1L] <- 999999
  expect_error(pretrend_bounds(f$post, other.ids), "same units")
  other.d <- f$pre; other.d$adapter_data$D <- 1 - other.d$adapter_data$D
  expect_error(pretrend_bounds(f$post, other.d), "membership")
  other.method <- f$pre; other.method$refit_spec$method <- "crossfit"
  expect_error(pretrend_bounds(f$post, other.method), "same method")
  bad.scale <- f$pre; bad.scale$results$main$treat[[1L]]$estimates$sigma2.s <- 0
  expect_error(pretrend_bounds(f$post, bad.scale, type = "scaled"), "strictly positive")
  expect_s3_class(pretrend_bounds(f$post, bad.scale, type = "bias"), "dml_pretrend_bounds")
  expect_error(pretrend_bounds(f$post, f$pre, k = -1), "nonnegative")
  expect_error(pretrend_bounds(f$post, f$pre, k = Inf), "finite")
  expect_error(pretrend_bounds(f$post, f$pre, k = numeric(0)), "nonempty")
  expect_error(pretrend_bounds(f$post, f$pre, level = 1), "one-sided")
})

test_that("repeated cross-fitting combines paired bounds without dividing SEs by repetitions", {
  f <- pretrend_fixture(180L, "crossfit")
  for (cmb in c("mean", "median")) {
    b <- pretrend_bounds(f$post, f$pre, k = 1.2, type = "scaled", combine.method = cmb)
    rr <- b$results[[1L]]
    est <- vapply(rr, function(x) x$theta.minus, numeric(1))
    se <- vapply(rr, function(x) x$se.minus, numeric(1))
    combine <- if (cmb == "mean") combine.mean else combine.median
    expected <- combine(est, se)
    expect_equal(b$bounds$theta.minus, unname(expected["estimate"]))
    expect_equal(b$bounds$se.minus, unname(expected["se"]))
    expect_length(b$joint.covariance, 2L)
  }
  bad <- f$pre; bad$adapter_info$fold_id[1, 1] <- 99L
  expect_error(pretrend_bounds(f$post, bad), "partitions differ")
  bad <- f$pre; bad$results$main$treat <- bad$results$main$treat[1L]
  expect_error(pretrend_bounds(f$post, bad), "same number")
})

pretrend_did_fixture <- function() {
  skip_if_not_installed("did")
  set.seed(10627)
  n <- 240L
  units <- data.frame(id = seq_len(n) * 3, g = rep(c(0, 3, 4), each = n / 3),
    x = runif(n, -1, 1), z = runif(n, -1, 1), u = rnorm(n))
  do.call(rbind, lapply(1:5, function(t) transform(units, time = t,
    y = u + t * (0.1 + 0.4 * x + 0.2 * z + 0.7 * (g == 3)) +
      0.9 * (g > 0 & t >= g) + rnorm(n))))
}

pretrend_did_fit <- function(dat, fast = TRUE, controls = "nevertreated",
                              base = "varying", anticipation = 0) {
  did::att_gt(yname = "y", tname = "time", idname = "id", gname = "g",
    xformla = ~ x + z, data = dat, panel = TRUE, bstrap = FALSE, cband = FALSE,
    faster_mode = fast, control_group = controls, base_period = base,
    anticipation = anticipation, compute_inffunc = TRUE)
}

test_that("did placebo uses the post cell sample under both designs and base conventions", {
  dat <- pretrend_did_fixture()
  for (fast in c(TRUE, FALSE)) for (base in c("varying", "universal")) {
    fit <- pretrend_did_fit(dat, fast, "notyettreated", base)
    m <- dml_from_did(fit, 3, 4, method = "parametric")
    old <- m$results
    b <- did_pretrend_bounds(m, c(1, 2), k = c(0, 1, 2), type = "scaled")
    expect_identical(m$results, old)
    expect_identical(b$model$results, old)
    expect_identical(b$unit_id, m$adapter_data$unit_id)
    expect_identical(b$pre_model$adapter_data$D, m$adapter_data$D)
    expect_identical(b$sample_policy, "post_cell")
    expect_identical(b$pre_model$adapter_info$mode, "parametric")
    # Independent construction from the original input, not context helpers.
    ids <- m$adapter_data$unit_id
    d1 <- dat[dat$time == 1, ]; d2 <- dat[dat$time == 2, ]
    d1 <- d1[match(ids, d1$id), ]; d2 <- d2[match(ids, d2$id), ]
    X <- model.matrix(~ x + z, d1)
    src <- DRDID::drdid_panel(d2$y, d1$y, as.numeric(d1$g == 3), X, inffunc = TRUE)
    expect_equal(b$bounds$theta.pre, rep(src$ATT, 3), tolerance = 1e-10)
    expect_equal(b$pre_model$results$main$treat[[1L]]$psis$psi.theta.s,
                  as.numeric(src$att.inf.func), tolerance = 1e-10)
    expect_true(all(d1$g %in% c(0, 3)))
    expect_gt(abs(b$bounds$theta.pre[1L]), 0.1)
    if (base == "universal") {
      j <- which(fit$group == 3 & fit$t == 2)
      expect_equal(fit$att[j], 0)
      expect_false(isTRUE(all.equal(b$bounds$theta.pre[1L], fit$att[j])))
    }
    restored <- unserialize(serialize(m, NULL))
    expect_equal(did_pretrend_bounds(restored, c(1, 2), type = "scaled")$bounds,
                  b$bounds[2L, , drop = FALSE], ignore_attr = TRUE)
  }
})

test_that("did wrapper supports crossfit batches and rejects invalid dates and changed sources", {
  dat <- pretrend_did_fixture()
  fit <- pretrend_did_fit(dat)
  args <- list(yreg = "lm", dreg = "glm", d.class = TRUE, cf.folds = 3,
               cf.reps = 2, cf.seed = 17, ps.trim = 0)
  m <- dml_from_did(fit, 3, 4, dml_args = args)
  b <- did_pretrend_bounds(m, c(1, 2), type = "scaled")
  expect_identical(b$pre_model$adapter_info$fold_id, m$adapter_info$fold_id)
  expect_length(b$results[[1L]], 2L)
  batch <- dml_from_did_cells(fit, data.frame(group = c(3, 4), time = c(4, 5)),
                              method = "parametric", on_error = "stop")
  out <- did_cell_apply(batch, did_pretrend_bounds, pre_period = c(1, 2), k = 1)
  expect_identical(out$status$status, c("ok", "ok"))
  expect_s3_class(out$values[[1L]], "dml_pretrend_bounds")
  expect_error(did_pretrend_bounds(m, c(2, 1)), "increasing")
  expect_error(did_pretrend_bounds(m, c(0, 1)), "absent")
  expect_error(did_pretrend_bounds(m, c(2, 3)), "anticipation")
  bad <- m; bad$adapter_data$y0[1L] <- bad$adapter_data$y0[1L] + 1
  expect_error(did_pretrend_bounds(bad, c(1, 2)), "no longer matches")
  bad <- m; bad$external_fit <- NULL
  expect_error(did_pretrend_bounds(bad, c(1, 2)), "original MP")
  anticip <- pretrend_did_fit(dat, anticipation = 1)
  ma <- dml_from_did(anticip, 4, 5, method = "parametric")
  expect_error(did_pretrend_bounds(ma, c(2, 3)), "anticipation")
  expect_s3_class(did_pretrend_bounds(ma, c(1, 2)), "dml_pretrend_bounds")
})
