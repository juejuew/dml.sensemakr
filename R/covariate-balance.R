
# covariate balance ---------------------------------------------------------
# Helpers for the covariate balance table described in Table 2 of
# Chernozhukov, Cinelli, Newey, Sharma and Syrgkanis (2026): standardized
# differences in means and variance-ratio deviations between treated and
# control units, for a set of observed covariates.

check_balance_model <- function(model) {
  if (!inherits(model, "dml")) {
    stop("model must be an object of class 'dml' created with dml().")
  }
  if (is.null(model$data$x) || is.null(model$data$d)) {
    stop("model does not contain the covariates ('x') and/or treatment ('d') it was fitted with.")
  }
  invisible(TRUE)
}

check_balance_treatment <- function(d) {
  d.value <- unique(d)
  binary   <- all(d.value %in% c(0, 1))
  if (!binary) {
    stop("Treatment 'd' must be binary (0 = control, 1 = treated) to compute covariate balance.")
  }
  invisible(TRUE)
}

check_balance_covariates <- function(x, covariates) {

  if (is.null(covariates)) {
    return(colnames(x))
  }

  which.not <- which(!covariates %in% colnames(x))
  if (length(which.not) > 0) {
    stop("Covariates not found: ", paste(covariates[which.not], collapse = ", "), ".")
  }
  covariates
}

# per-covariate group statistics (skeleton only) -----------------------------
covariate_group_stats <- function(xj, d) {
  # TODO: compute mean and variance of xj separately for d == 0 (control)
  # and d == 1 (treated).
  list(mean.control = NA_real_,
      var.control  = NA_real_,
      mean.treated = NA_real_,
      var.treated  = NA_real_)
}

# standardized mean difference, (mu1 - mu0)/sqrt((s1^2 + s0^2)/2) -----------
std_mean_diff <- function(mean.control, mean.treated, var.control, var.treated) {
  # TODO: implement the standardized mean difference reported in Table 2 of
  # Chernozhukov et al (2026).
  NA_real_
}

# deviation of the variance ratio from one, s1^2/s0^2 - 1 -------------------
var_ratio_dev <- function(var.control, var.treated) {
  # TODO: implement the variance-ratio deviation reported in Table 2 of
  # Chernozhukov et al (2026).
  NA_real_
}

# builds the balance table for a set of covariates ---------------------------
balance_fun <- function(x, d, covariates) {

  rows <- lapply(covariates, function(covar) {
    stats <- covariate_group_stats(x[, covar], d)
    smd <- std_mean_diff(mean.control = stats$mean.control,
                         mean.treated = stats$mean.treated,
                         var.control  = stats$var.control,
                         var.treated  = stats$var.treated)
    vrd <- var_ratio_dev(var.control = stats$var.control,
                        var.treated = stats$var.treated)
    c(mean.control  = stats$mean.control,
      var.control   = stats$var.control,
      mean.treated  = stats$mean.treated,
      var.treated   = stats$var.treated,
      std.mean.diff = smd,
      var.ratio.dev = vrd)
  })

  out <- do.call(rbind, rows)
  rownames(out) <- covariates
  out
}

##' Covariate balance table for debiased machine learning
##' @description
##' Computes a covariate balance table between treated and control units, following Table 2 of Chernozhukov, Cinelli, Newey, Sharma and Syrgkanis (2026). For each covariate, the table reports the group means and variances (control and treated), the standardized difference in means, and the deviation of the variance ratio from one.
##' @param model an object of class \code{\link{dml}}. The treatment stored in \code{model} must be binary (0 = control, 1 = treated).
##' @param covariates character vector with the names of the covariates to be included in the balance table. Default (\code{NULL}) uses all covariates used to fit \code{model}.
##' @param ... arguments passed to other methods.
##' @returns An object of class \code{dml_balance} containing the covariate balance table.
##' @export
dml_balance <- function(model, covariates = NULL, ...) {

  check_balance_model(model)

  x <- model$data$x
  d <- model$data$d

  check_balance_treatment(d)
  covariates <- check_balance_covariates(x, covariates)

  out <- list()
  out$info <- list(covariates = covariates,
                   n.control  = sum(d == 0),
                   n.treated  = sum(d == 1))

  out$table <- balance_fun(x = x, d = d, covariates = covariates)

  class(out) <- "dml_balance"
  return(out)
}

##' Print and summary methods for DML covariate balance
##' @description Print and summary methods for objects of class \code{dml_balance} created with \code{\link{dml_balance}}.
##' @param x an object of class \code{\link{dml_balance}} or \code{summary_dml_balance}.
##' @param digits minimal number of significant digits.
##' @rdname summary.dml_balance
##' @export
print.dml_balance <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  print(summary(x), digits = digits, ...)
}

##' @param object an object of class \code{\link{dml_balance}}.
##' @param ... arguments passed to other methods.
##' @returns For \code{print}: the object, printed to console. For \code{summary}: the object of class \code{summary_dml_balance} holding the balance table.
##' @rdname summary.dml_balance
##' @export
summary.dml_balance <- function(object, ...) {
  out <- object
  class(out) <- "summary_dml_balance"
  out
}

##' @rdname summary.dml_balance
##' @export
print.summary_dml_balance <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\n")
  cat("Covariate Balance:", x$info$n.control, "control units,", x$info$n.treated, "treated units\n")
  cat("\n")
  print(round(x$table, digits), ...)
  cat("\nNote: columns report control/treated means and variances, the standardized mean difference, and the variance ratio deviation from one.\n")
}
