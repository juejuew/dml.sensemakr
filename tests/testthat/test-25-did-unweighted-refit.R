unweighted_did_data <- function() {
  skip_if_not_installed("did")
  skip_if_not_installed("DRDID")
  set.seed(908)
  n <- 300L
  units <- data.frame(id = seq_len(n) * 7, g = rep(c(0, 5, 9), each = n / 3),
                     x = runif(n, -1, 1), z = runif(n, -1, 1), u = rnorm(n))
  dat <- do.call(rbind, lapply(c(0, 2, 5, 9, 12), function(t) {
    d <- units
    d$time <- t
    d$x <- d$x + t * d$z / 30
    d$y <- d$u + t / 5 + (0.5 + t / 8) * d$x + d$z +
      (d$g > 0 & t >= d$g) * (0.7 + t / 20) + rnorm(n)
    d
  }))
  dat[sample.int(nrow(dat)), ]
}

unweighted_did_source <- function(data, ...) {
  did::att_gt(yname = "y", tname = "time", idname = "id", gname = "g",
              xformla = ~ x + z, data = data, panel = TRUE, bstrap = FALSE,
              cband = FALSE, ...)
}

unweighted_did_args <- function(reps = 3L, trim = 0) {
  list(yreg = "lm", dreg = "glm", d.class = TRUE, cf.folds = 3, cf.seed = 17,
         cf.reps = reps, ps.trim = trim, model = "npm", target = "att",
         conditional = TRUE, centered = FALSE, dirty.tuning = FALSE,
         y.class = FALSE, verbose = FALSE)
}

# Build the oracle from the original panel, independently of adapter extraction.
unweighted_did_oracle <- function(dat, group, time, args,
                                  controls = "nevertreated", base = 2) {
  eligible <- dat$g == group | dat$g == 0 |
    (controls == "notyettreated" & dat$g > time)
  pre <- dat[eligible & dat$time == base, ]
  post <- dat[eligible & dat$time == time, ]
  pre <- pre[order(pre$id), ]
  post <- post[match(pre$id, post$id), ]
  x <- stats::model.matrix(~ x + z, data = pre)
  rownames(x) <- NULL
  attr(x, "assign") <- NULL # did's row-subset design is a plain numeric matrix.
  do.call(dml, c(list(y = post$y - pre$y, d = as.numeric(pre$g == group), x = x), args))
}

test_that("unequal weights give the same complete fit as direct unweighted DML", {
  dat <- unweighted_did_data()
  dat$w <- 1 + dat$id / max(dat$id)
  for (fast in c(TRUE, FALSE)) {
    for (controls in c("nevertreated", "notyettreated")) {
      args <- unweighted_did_args(if (fast) 3L else 1L)
      source <- unweighted_did_source(dat, faster_mode = fast, control_group = controls,
                                 weightsname = "w")
      expect_warning(actual <- dml_from_did(source, 5, 5, args), "Unequal sampling weights")
      expected <- unweighted_did_oracle(dat, 5, 5, args, controls)
      expect_identical(class(actual), "dml")
      expect_identical(actual$adapter_info$mode, "unweighted_dml")
      expect_equal(actual$data, expected$data, tolerance = 1e-12)
      expect_equal(actual$results, expected$results, tolerance = 1e-12)
      expect_equal(actual$coefs, expected$coefs, tolerance = 1e-12)
      expect_equal(lapply(actual$fits, `[[`, "preds"),
                   lapply(expected$fits, `[[`, "preds"), tolerance = 1e-12)
      expect_equal(se(actual), se(expected), tolerance = 1e-12)
      expect_equal(length(actual$results$main$treat), args$cf.reps)
      expect_true(all(actual$adapter_data$i.weights == 1))
      expect_true(any(actual$adapter_info$source_weights != 1))
      expect_false(isTRUE(all.equal(as.numeric(coef(actual)),
        source$att[source$group == 5 & source$t == 5], tolerance = 1e-8)))
      expect_false(inherits(actual, "dml_did"))
    }
  }
})

test_that("unweighted refits retain native benchmarking after serialization", {
  dat <- unweighted_did_data()
  dat$w <- 1 + dat$id / max(dat$id)
  source <- unweighted_did_source(dat, weightsname = "w")
  args <- unweighted_did_args()
  actual <- suppressWarnings(dml_from_did(source, 5, 9, args))
  actual <- unserialize(serialize(actual, NULL))
  expected <- unweighted_did_oracle(dat, 5, 9, args)
  # Any accidental hybrid dispatch would try to use this invalid source fit.
  actual$external_fit <- NULL
  b.actual <- dml_benchmark(actual, "z")
  b.expected <- dml_benchmark(expected, "z")
  expect_equal(b.actual, b.expected, tolerance = 1e-12)
  for (combine in c("mean", "median")) {
    expect_equal(confidence_bounds(actual, cf.y = 0.05, cf.d = 0.05,
                                   combine.method = combine),
                 confidence_bounds(expected, cf.y = 0.05, cf.d = 0.05,
                                   combine.method = combine), tolerance = 1e-12)
    expect_equal(benchmark_bounds(actual, b.actual, combine.method = combine),
                 benchmark_bounds(expected, b.expected, combine.method = combine),
                 tolerance = 1e-12)
    expect_equal(robustness_value(actual, combine.method = combine),
                 robustness_value(expected, combine.method = combine), tolerance = 1e-12)
  }
})

test_that("unweighted batch refits do not require or reuse source ATT and IF", {
  dat <- unweighted_did_data()
  dat$w <- 1 + dat$id / max(dat$id)
  source <- unweighted_did_source(dat, weightsname = "w", compute_inffunc = FALSE)
  source$att[] <- NA_real_
  args <- unweighted_did_args(1L)
  cells <- data.frame(group = c(5, 9), time = c(9, 12))
  batch <- dml_from_did_cells(source, cells, args)
  expect_identical(batch$status$status, c("ok", "ok"))
  expect_true(all(grepl("Unequal sampling weights", batch$status$warnings)))
  output <- did_cell_apply(batch, sensemakr, benchmark_covariates = "z")
  expect_identical(output$status$status, c("ok", "ok"))
  for (i in seq_len(nrow(cells))) {
    direct <- unweighted_did_oracle(dat, cells$group[i], cells$time[i], args,
                                    base = c(2, 5)[i])
    expect_equal(batch$models[[i]]$results, direct$results, tolerance = 1e-12)
    expect_equal(output$values[[i]]$sensitivity_stats,
                 sensemakr(direct)$sensitivity_stats, tolerance = 1e-12)
  }
})

test_that("unweighted refits use native DML clipping and still validate saved data", {
  dat <- unweighted_did_data()
  dat$w <- 1 + dat$id / max(dat$id)
  source <- unweighted_did_source(dat, weightsname = "w")
  args <- unweighted_did_args(1L, trim = 0.49)
  expect_warning(actual <- dml_from_did(source, 5, 9, args), "Unequal sampling weights")
  direct <- unweighted_did_oracle(dat, 5, 9, args)
  expect_equal(actual$results, direct$results, tolerance = 1e-12)
  expect_gt(actual$results$main$treat[[1]]$trim.summary$trimmed_num$all, 0)
  for (bad.weight in c(0, -1, NA_real_, Inf)) {
    bad <- source
    bad$DIDparams$data <- as.data.frame(bad$DIDparams$data)
    bad$DIDparams$data$.w[1] <- bad.weight
    expect_error(dml_from_did(bad, 5, 9, args), "finite and positive")
  }
  bad <- source
  bad$DIDparams$data <- as.data.frame(bad$DIDparams$data)[-1, ]
  expect_error(dml_from_did(bad, 5, 9, args), "balanced panel")
})
