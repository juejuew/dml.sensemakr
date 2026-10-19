# Full-sample, finite-dimensional sensitivity moments for the traditional
# DRDID panel estimator. See inst/notes/parametric-inference.md for the
# estimating equations and the interpretation under model misspecification.
.panel_parametric_args <- function(args) {
  if (!is.list(args) || length(args))
    .drdid_stop("method = 'parametric' requires an empty dml_args list: source logistic/OLS models are fixed and cross-fitting is not used.")
  invisible(NULL)
}

.panel_parametric_inverse <- function(H, label) {
  if (any(!is.finite(H)) || rcond(H) < .Machine$double.eps)
    .drdid_stop("The full-sample ", label, " information matrix is singular.")
  tryCatch(chol2inv(chol(H)), error = function(e)
    .drdid_stop("The full-sample ", label, " information matrix is not positive definite."))
}

.panel_parametric_aux <- function(sample, theta, phi, spec, diagnostics) {
  if (!identical(spec$method, "trad"))
    .drdid_stop("method = 'parametric' currently supports the traditional logistic-MLE/OLS estimator only (estMethod = 'trad').")
  .drdid_check_weights(sample$i.weights, length(sample$D))
  y <- sample$y1 - sample$y0
  d <- sample$D
  X <- sample$X
  n <- length(d)
  control <- d == 0
  if (sum(control) <= ncol(X))
    .drdid_stop("Full-sample OLS requires more control units than design columns for regular sensitivity inference.")

  # Use DRDID's actual solvers and the saved design, not caret preprocessing
  # or a newly evaluated formula. This matches drdid_panel()'s nuisance fits.
  fit.fun <- .drdid_function("fastglm_fit", c("x", "y", "family", "weights", "method"))
  ps.fit <- fit.fun(x = X, y = d, family = stats::binomial(),
                    weights = sample$i.weights, method = 3)
  or.fit <- fit.fun(x = X[control, , drop = FALSE], y = y[control],
                    family = stats::gaussian(link = "identity"),
                    weights = sample$i.weights[control])
  if (!isTRUE(ps.fit$converged) ||
      length(ps.fit$coefficients) != ncol(X) ||
      length(or.fit$coefficients) != ncol(X) ||
      any(!is.finite(c(ps.fit$coefficients, or.fit$coefficients))))
    .drdid_stop("Full-sample logistic/OLS fitting failed or returned invalid coefficients.")
  e <- as.numeric(ps.fit$fitted.values)
  g <- as.numeric(X %*% or.fit$coefficients)
  if (length(e) != n || any(!is.finite(c(e, g))) ||
      any(e <= 0 | e >= 1 - 1e-6) || any(e[control] >= spec$trim.level))
    .drdid_stop("Source propensity trimming/clipping is active; parametric sensitivity inference requires the untrimmed panel estimand.")
  close <- function(a, b) isTRUE(all.equal(as.numeric(a), as.numeric(b),
                                          tolerance = 1e-7, check.attributes = FALSE))
  if (!close(e, diagnostics$propensity))
    .drdid_stop("Full-sample propensities do not reproduce the checked source nuisance fit.")
  residual <- y - g
  odds <- e / (1 - e)
  reconstructed <- mean(residual[d == 1]) -
    sum(odds[control] * residual[control]) / sum(odds[control])
  if (!close(reconstructed, theta))
    .drdid_stop("Full-sample nuisance fits do not reproduce the source ATT.")

  # Sandwich IFs for the nuisance estimating equations. These describe the
  # fitted parametric functionals even at regular pseudo-true parameters;
  # they do not make misspecified sensitivity targets causally valid.
  H.beta <- crossprod(X, (1 - d) * X) / n
  H.gamma <- crossprod(X, e * (1 - e) * X) / n
  psi.beta <- ((1 - d) * residual * X) %*%
    .panel_parametric_inverse(H.beta, "outcome")
  psi.gamma <- ((d - e) * X) %*%
    .panel_parametric_inverse(H.gamma, "propensity")
  p <- mean(d)
  psi.p <- d - p

  result <- att.npm.cond(y, d, yhat0 = g, dhat = e,
                         phat = rep(p, n), trim = 0, centered = FALSE)
  # Differentiate the actual normalized moment estimators, keeping the
  # empirical distribution fixed. The ratio for sigma2 cancels p entirely.
  grad.beta <- -2 * colMeans((1 - d) * residual * X) / mean(1 - d)
  correction.sigma <- as.numeric(psi.beta %*% grad.beta)
  l <- d / p
  lc <- (1 - d) / (1 - p)
  h <- odds * (1 - p) / p
  denominator <- mean(2 * l - lc)
  nu <- result$estimates$nu2.s
  grad.gamma <- 2 * colMeans((l * h - lc * h^2) * X) / denominator
  dl <- -l / p
  dlc <- lc / (1 - p)
  dh <- -h / (p * (1 - p))
  grad.p <- mean(2 * dl * h + 2 * l * dh - dlc * h^2 -
                   2 * lc * h * dh - nu * (2 * dl - dlc)) / denominator
  # grad.p is algebraically zero at p = mean(d) for this ratio score.
  # Keeping it explicit makes the full chain rule auditable and testable.
  correction.nu <- as.numeric(psi.gamma %*% grad.gamma) + grad.p * psi.p
  result$psis$psi.sigma2.s <- result$psis$psi.sigma2.s + correction.sigma
  result$psis$psi.nu2.s <- result$psis$psi.nu2.s + correction.nu
  result$nuisance.correction <- list(psi.sigma2.s = correction.sigma,
                                      psi.nu2.s = correction.nu)
  if (!is.finite(result$estimates$S2) || result$estimates$S2 <= 0 ||
      any(!is.finite(c(result$psis$psi.sigma2.s, result$psis$psi.nu2.s))))
    .drdid_stop("Full-sample fitting produced invalid sensitivity moments or IFs; regular sensitivity inference is unavailable.")
  # Preserve the original estimator's ATT/IF, including its normalization and
  # nuisance corrections. Only the newly computed moment IFs are supplied here.
  result <- .drdid_combine_moments(result, theta, phi, y, d)
  result$trim.summary$trimmed_obs <-
    collect.trimmed.obs(y, d, X, result$trim.summary$trimmed_indices)

  fit <- list(data = list(y = y, d = d, x = X),
    info = list(model = "npm", target = "att", conditional = TRUE,
      centered = FALSE, cf.folds = 0L, cf.reps = 1L, cross.fitting = FALSE,
      nuisance.method = "parametric", dirty.tuning = FALSE,
      yreg = list(yreg0 = list(method = "lm"), yreg1 = NULL),
      dreg = list(method = "glm")),
    fits = list(list(preds = list(yhat0 = g, yhat1 = NULL, dhat = e,
                                  phat = rep(p, n)))),
    results = list(main = list(treat = list(result))),
    coefs = list(main = list(treat = combine.cross.fits(list(result)))),
    parametric = list(coefficients = list(outcome = or.fit$coefficients,
                                          propensity = ps.fit$coefficients),
      treated.share = p, reconstructed.att = reconstructed,
      moment.gradients = list(sigma2.beta = grad.beta, nu2.gamma = grad.gamma,
                               nu2.p = grad.p),
      inference = "joint IF with full-sample nuisance estimation effects",
      assumptions = c("fixed-dimensional regular logistic/linear models",
        "correct nuisance functions for the sensitivity target",
        "overlap and finite moments", "independent panel units")))
  class(fit) <- "dml"
  list(fit = fit, args = list(), fold = integer(0), rep_seeds = integer(0))
}
