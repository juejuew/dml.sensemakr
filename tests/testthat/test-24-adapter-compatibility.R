compat_source <- function(fast = TRUE) {
  skip_if_not_installed("did")
  skip_if_not_installed("DRDID")
  env <- new.env()
  utils::data("mpdta", package = "did", envir = env)
  did::att_gt(yname = "lemp", tname = "year", idname = "countyreal",
    gname = "first.treat", xformla = ~ lpop, data = env$mpdta,
    faster_mode = fast, bstrap = FALSE, cband = FALSE, compute_inffunc = TRUE)
}

test_that("real source objects work with either design layout and three repetitions", {
  for (fast in c(TRUE, FALSE)) {
    fit <- compat_source(fast)
    m <- dml_from_did(fit, 2004, 2007, list(yreg = "lm", dreg = "glm",
      d.class = TRUE, cf.folds = 3, cf.reps = 3, cf.seed = 17, ps.trim = 0))
    j <- which(fit$group == 2004 & fit$t == 2007)
    expect_length(m$results$main$treat, 3L)
    expect_equal(unname(coef(m)), fit$att[j])
    expect_identical(m$adapter_info$reconstruction.version,
                     as.character(utils::packageVersion("did")))
    expect_identical(m$adapter_info$DRDID.version,
                     as.character(utils::packageVersion("DRDID")))
    expect_identical(m$adapter_info$compatibility$layout,
                     if (fast) "panel_tensors" else "panel_formula")
    expect_true(m$adapter_info$diagnostics$reproduced)
    for (r in m$results$main$treat) {
      ids <- match(as.character(m$adapter_data$unit_id), rownames(fit$inffunc))
      expect_equal(r$psis$psi.theta.s,
        as.numeric(fit$inffunc[ids, j]) * length(ids) / fit$n)
    }
    expect_s3_class(sensemakr(m, cf.y = .05, cf.d = .05), "dml.sensemakr")
  }
})

test_that("installed version labels do not determine acceptance of saved results", {
  fit <- compat_source()
  expected <- .did_cell_sample(.did_adapter_context(fit), 2004, 2007)
  # A new dependency label must not bypass or pre-empt the content checks.
  local_mocked_bindings(.adapter_package_version = function(package) "99.0.0")
  ctx <- .did_adapter_context(fit)
  expect_identical(ctx$version, "99.0.0")
  expect_identical(ctx$drdid.version, "99.0.0")
  actual <- .did_cell_sample(ctx, 2004, 2007)
  expect_equal(actual, expected)
  bad <- fit; bad$inffunc <- NULL
  expect_error(.did_adapter_context(bad), "source has no IF")
  bad <- fit; j <- which(bad$group == 2004 & bad$t == 2007)
  bad$att[j] <- bad$att[j] + .2
  expect_error(.did_cell_sample(.did_adapter_context(bad), 2004, 2007), "ordered IF")
})

test_that("missing metadata is reported rather than replaced with assumed defaults", {
  fit <- compat_source(FALSE)
  for (key in c("faster_mode", "control_group", "anticipation", "xformla", "panel")) {
    bad <- fit; bad$DIDparams[[key]] <- NULL
    expect_error(.did_adapter_context(bad), paste0("required DIDparams fields: ", key))
  }
  bad <- fit; bad$DIDparams$idname <- c("countyreal", "year")
  expect_error(.did_adapter_context(bad), "column names")
  bad <- fit; bad$DIDparams$data <- NULL
  expect_error(.did_adapter_context(bad), "DIDparams\\$data")
  bad <- fit; rownames(bad$inffunc) <- NULL
  expect_error(.did_adapter_context(bad), "unique row names")
  bad <- fit
  changed <- which(bad$DIDparams$data$year == 2003 & bad$DIDparams$data$first.treat == 2004)[1]
  bad$DIDparams$data$lemp[changed] <- bad$DIDparams$data$lemp[changed] + 1
  expect_error(.did_cell_sample(.did_adapter_context(bad), 2004, 2007), "ordered IF")
})

test_that("source IF scale and unit labels remain mandatory", {
  fit <- compat_source()
  j <- which(fit$group == 2004 & fit$t == 2007)
  bad <- fit; bad$inffunc[, j] <- 2 * bad$inffunc[, j]
  expect_error(.did_cell_sample(.did_adapter_context(bad), 2004, 2007), "ordered IF")
  bad <- fit; rownames(bad$inffunc) <- rev(rownames(bad$inffunc))
  expect_error(.did_cell_sample(.did_adapter_context(bad), 2004, 2007), "outside|ordered IF")
})

test_that("missing and changed DRDID interfaces produce specific errors", {
  skip_if_not_installed("DRDID")
  expect_error(.drdid_function("adapter_nonexistent_function", "x"), "required function")
  real.get <- utils::getFromNamespace
  local_mocked_bindings(getFromNamespace = function(x, ns, ...) {
    if (identical(ns, "DRDID") && identical(x, "fastglm_fit"))
      return(function(x, y, ...) NULL)
    real.get(x, ns, ...)
  }, .package = "utils")
  expect_error(.drdid_function("fastglm_fit", c("x", "y", "family", "weights", "method")),
    "Incompatible DRDID interface fastglm_fit: missing arguments family, weights, method")
})

test_that("empty propensity outputs cannot bypass the overlap check", {
  src <- .did_cell_sample(.did_adapter_context(compat_source()), 2004, 2007)
  real.get <- utils::getFromNamespace
  local_mocked_bindings(getFromNamespace = function(x, ns, ...) {
    if (identical(ns, "DRDID") && identical(x, "fastglm_fit"))
      return(function(x, y, family, weights, method) list(converged = TRUE, fitted.values = numeric()))
    real.get(x, ns, ...)
  }, .package = "utils")
  expect_error(.panel_check_source(src$theta, src$sample,
    list(method = "trad", trim.level = .995, phi = src$phi)), "one finite probability per unit")
})

test_that("changed estimator output cannot be accepted as an ATT/IF replay", {
  src <- .did_cell_sample(.did_adapter_context(compat_source()), 2004, 2007)
  real.get <- utils::getFromNamespace
  local_mocked_bindings(getFromNamespace = function(x, ns, ...) {
    if (identical(ns, "DRDID") && identical(x, "drdid_panel"))
      return(function(y1, y0, D, covariates, i.weights, boot, inffunc, trim.level)
        list(ATT = 0, att.inf.func = numeric(length(D) - 1L)))
    real.get(x, ns, ...)
  }, .package = "utils")
  expect_error(.panel_check_source(src$theta, src$sample,
    list(method = "trad", trim.level = .995, phi = src$phi)), "one finite IF value per unit")
})
