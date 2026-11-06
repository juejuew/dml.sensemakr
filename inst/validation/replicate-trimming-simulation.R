# Run after installing/loading dml.sensemakr >= 0.2.2.9009:
# source('replicate-trimming-simulation.R')
# ans <- replicate_trimming_simulation(B = 500L, output = 'trimming-validation')
# This is a numerical coverage check, not a proof of general validity.

replicate_trimming_simulation <- function(B = 500L, sizes = c(1000L, 4000L),
                                          seed = 31033L, output = NULL) {
  stopifnot(requireNamespace("dml.sensemakr", quietly = TRUE),
            requireNamespace("DRDID", quietly = TRUE),
            length(B) == 1L, is.finite(B), B >= 2L, B == as.integer(B),
            all(is.finite(sizes)), all(sizes >= 20L), all(sizes == as.integer(sizes)))
  # X is continuous across the cutoff, so the sample trim mask is not fixed.
  # c=.7 creates substantial trimming at moderate n; integration tests separately
  # exercise did's actual default c=.995. Models are correctly specified.
  cutoff <- 0.7
  prop <- function(x) plogis(0.4 + 1.6 * x)
  p <- integrate(function(x) prop(x)/2, -1, 1, rel.tol = 1e-12)$value
  a <- 0.75^2
  b <- (1-p)/p^2 * integrate(function(x) prop(x)^2/(1-prop(x))/2,
                              -1, 1, rel.tol = 1e-12)$value
  k <- sqrt(0.08 * 0.06 / (1 - 0.06))
  truth <- c(ATT = 0.8, sigma2 = a, nu2 = b, lower = 0.8 - k*sqrt(a*b))
  set.seed(seed)
  runs <- vector("list", length(sizes)*B)
  failures <- list(); j <- 0L
  for (n in sizes) for (rep in seq_len(B)) {
    j <- j + 1L
    x <- runif(n, -1, 1); X <- cbind(`(Intercept)` = 1, x = x)
    d <- rbinom(n, 1, prop(x)); y0 <- rnorm(n)
    y1 <- y0 + 0.2 + 0.5*x + 0.8*d + rnorm(n, sd = 0.75)
    run <- tryCatch({
      source <- DRDID::drdid_panel(y1, y0, d, X, trim.level = cutoff, inffunc = TRUE)
      m <- dml.sensemakr::dml_from_drdid(source, y1 = y1, y0 = y0, D = d,
                                        covariates = X, method = "replicate")
      r <- m$results$main$treat[[1]]
      bound <- dml.sensemakr::dml_bounds(m, cf.y = 0.08, cf.d = 0.06)$results$main$treat[[1]]
      estimates <- c(r$estimates$theta.s, r$estimates$sigma2.s, r$estimates$nu2.s,
                      bound$estimates$theta.m)
      standard.errors <- c(sqrt(diag(r$joint.covariance)), bound$estimates$se.theta.m)
      stopifnot(all(is.finite(estimates)), all(is.finite(standard.errors)), all(standard.errors > 0))
      data.frame(n = n, replication = rep, parameter = names(truth),
        truth = unname(truth), estimate = unname(estimates), se = unname(standard.errors),
        trimmed.controls = m$adapter_info$diagnostics$source.trimming$n.trimmed.controls,
        controls.near.cutoff = sum(d == 0 & abs(prop(x)-cutoff) < 0.02))
    }, error = function(e) e)
    if (inherits(run, "error")) {
      failures[[length(failures)+1L]] <- data.frame(n = n, replication = rep, message = conditionMessage(run))
    } else runs[[j]] <- run
  }
  raw <- do.call(rbind, runs)
  if (is.null(raw)) stop("All simulation runs failed.")
  fail <- if (length(failures)) do.call(rbind, failures) else
    data.frame(n = integer(), replication = integer(), message = character())
  summary <- do.call(rbind, lapply(split(raw, list(raw$n, raw$parameter), drop = TRUE), function(z) {
    cover <- abs(z$estimate-z$truth) <= qnorm(0.975)*z$se
    data.frame(n = z$n[1], parameter = z$parameter[1], requested = B,
      succeeded = nrow(z), failed = sum(fail$n == z$n[1]), truth = z$truth[1],
      bias = mean(z$estimate-z$truth), empirical.sd = sd(z$estimate), mean.se = mean(z$se),
      coverage.95 = mean(cover), coverage.mcse = sqrt(mean(cover)*(1-mean(cover))/nrow(z)),
      mean.trimmed.controls = mean(z$trimmed.controls),
      mean.controls.near.cutoff = mean(z$controls.near.cutoff))
  }))
  rownames(summary) <- NULL
  if (!is.null(output)) {
    dir.create(output, recursive = TRUE, showWarnings = FALSE)
    write.csv(summary, file.path(output, "coverage-summary.csv"), row.names = FALSE)
    write.csv(raw, file.path(output, "coverage-replications.csv"), row.names = FALSE)
    write.csv(fail, file.path(output, "coverage-failures.csv"), row.names = FALSE)
  }
  if (nrow(fail)) warning(nrow(fail), " runs failed; inspect failures. No successful run was silently selected.")
  list(summary = summary, failures = fail, replications = raw,
       settings = list(B = B, sizes = sizes, seed = seed, cutoff = cutoff, truth = truth))
}
