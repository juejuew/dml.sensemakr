
# covariate balance ---------------------------------------------------------
# Helpers for the covariate balance table described in Table 2 of
# Omitted Variable Bias in Difference-in-Differences Designs (2026): standardized
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

# per-covariate group statistics ---------------------------------------------
covariate_group_stats <- function(xj, d) {
  list(mean.control = mean(xj[d == 0]),
      var.control  = var(xj[d == 0]),
      mean.treated = mean(xj[d == 1]),
      var.treated  = var(xj[d == 1]))
}

# standardized mean difference, (mu1 - mu0)/sqrt((s1^2 + s0^2)/2) -----------
std_mean_diff <- function(mean.control, mean.treated, var.control, var.treated) {
  (mean.treated - mean.control) / sqrt((var.treated + var.control) / 2)
}

# deviation of the variance ratio from one, s1^2/s0^2 - 1 -------------------
var_ratio_dev <- function(var.control, var.treated) {
  var.treated / var.control - 1
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

# Formatted rendering of the balance table (Table 2 layout), in LaTeX or
# HTML -----------------------------------------------------------------------
# The LaTeX renderer uses standard $...$ math, compiled by LaTeX itself. The
# HTML renderer deliberately does NOT reuse that notation: R Markdown's
# html_document only loads MathJax via a <script> injected at render time
# that fetches https://mathjax.rstudio.com/... at runtime, which fails
# silently when the file is opened offline or directly from disk (many
# browsers block scripts on file:// pages outright). Instead, the HTML
# symbols/formulas below are built from plain HTML entities and CSS, so the
# output renders identically with no JavaScript or network dependency. The
# LaTeX Std. Mean Diff./Var. Ratio Dev. columns use siunitx's `S` column
# type so that values align on the decimal point regardless of sign; the
# document's LaTeX preamble must include \usepackage{siunitx} in addition
# to \usepackage{booktabs}.

# LaTeX math notation for the header symbols/formulas
mu0.sym      <- "$\\mu_0$"
sigma0.sym   <- "$\\sigma_0^2$"
mu1.sym      <- "$\\mu_1$"
sigma1.sym   <- "$\\sigma_1^2$"
smd.formula  <- "$\\frac{\\mu_1-\\mu_0}{\\sqrt{(\\sigma_1^2+\\sigma_0^2)/2}}$"
vrd.formula  <- "$\\frac{\\sigma_1^2}{\\sigma_0^2}-1$"

# a self-contained CSS "stacked fraction" (numerator over a rule over
# denominator), used to build the HTML formulas below without any
# JavaScript/MathJax dependency
html_fraction <- function(numerator, denominator) {
  paste0("<span style=\"display:inline-block; vertical-align:middle; text-align:center;\">",
        "<span style=\"display:block; border-bottom:1px solid black; padding:0 2px;\">",
        numerator, "</span>",
        "<span style=\"display:block; padding:0 2px;\">", denominator, "</span>",
        "</span>")
}

# HTML notation for the header symbols/formulas (plain entities + CSS, no
# MathJax needed)
mu0.html    <- "&mu;<sub>0</sub>"
sigma0.html <- "&sigma;<sub>0</sub><sup>2</sup>"
mu1.html    <- "&mu;<sub>1</sub>"
sigma1.html <- "&sigma;<sub>1</sub><sup>2</sup>"
smd.formula.html <- html_fraction(
  paste0(mu1.html, " &minus; ", mu0.html),
  paste0("&radic;((", sigma1.html, " + ", sigma0.html, ")/2)")
)
vrd.formula.html <- paste0(html_fraction(sigma1.html, sigma0.html), " &minus; 1")

# escapes the LaTeX special characters that may show up in covariate names,
# captions, or labels
escape_latex <- function(text) {
  gsub("([&%$#_{}])", "\\\\\\1", text)
}

# escapes the HTML special characters that may show up in covariate names,
# captions, or labels (the latter may end up inside an id="..." attribute,
# hence also escaping quotes)
escape_html <- function(text) {
  text <- gsub("&", "&amp;",  text, fixed = TRUE)
  text <- gsub("<", "&lt;",   text, fixed = TRUE)
  text <- gsub(">", "&gt;",   text, fixed = TRUE)
  text <- gsub("\"", "&quot;", text, fixed = TRUE)
  text
}

# shared numeric formatting for finite values: fixed decimals for ordinary
# magnitudes, thousands separator for large ones (used by both renderers)
format_finite_num <- function(value, digits, big.threshold = 1000) {
  if (abs(value) >= big.threshold) {
    return(formatC(round(value), format = "d", big.mark = ","))
  }
  formatC(value, digits = digits, format = "f")
}

# the placeholder strings used for values that are not plain numbers -- an
# `S` column cannot parse these directly, so cells containing them must be
# wrapped in \multicolumn{1}{c}{...} (see latex_balance_rows())
latex_na    <- "--"
latex_inf   <- "$\\infty$"
latex_ninf  <- "$-\\infty$"

# formats a single numeric value for LaTeX: NA/NaN and +-Inf become
# plain-text placeholders, everything else uses format_finite_num()
format_latex_value <- function(value, digits, big.threshold = 1000) {
  if (is.na(value)) return(latex_na)
  if (is.infinite(value)) return(if (value > 0) latex_inf else latex_ninf)
  format_finite_num(value, digits, big.threshold)
}

# formats a numeric vector/matrix of values for LaTeX (see format_latex_value)
format_latex_num <- function(values, digits, big.threshold = 1000) {
  vapply(values, format_latex_value, character(1),
        digits = digits, big.threshold = big.threshold)
}

# formats a single numeric value for HTML: NA/NaN and +-Inf become HTML
# entities, everything else uses format_finite_num()
format_html_value <- function(value, digits, big.threshold = 1000) {
  if (is.na(value)) return("--")
  if (is.infinite(value)) return(if (value > 0) "&infin;" else "&minus;&infin;")
  format_finite_num(value, digits, big.threshold)
}

# formats a numeric vector/matrix of values for HTML (see format_html_value)
format_html_num <- function(values, digits, big.threshold = 1000) {
  vapply(values, format_html_value, character(1),
        digits = digits, big.threshold = big.threshold)
}

# wraps non-numeric placeholders (NA/Inf) so they render correctly inside a
# siunitx `S` column; plain numbers are left untouched
wrap_if_placeholder <- function(values) {
  placeholder <- values %in% c(latex_na, latex_inf, latex_ninf)
  ifelse(placeholder, paste0("\\multicolumn{1}{c}{", values, "}"), values)
}

# the siunitx `S[table-format=...]` column spec needed to fit every finite
# value in a column, given a fixed number of decimals
siunitx_format <- function(values, digits) {
  finite.values <- values[is.finite(values)]
  int.digits <- if (length(finite.values) == 0) 1L else
    max(1L, nchar(formatC(floor(abs(finite.values)), format = "d")))
  paste0("S[table-format=-", int.digits, ".", digits, "]")
}

# LaTeX header rows: group labels, per-column rules, then the (mu, sigma^2) /
# formula labels. Std. Mean Diff./Var. Ratio Dev. cells are wrapped in
# \multicolumn{1}{c}{...} because they sit in `S` columns, which only
# accept plain numbers in the body rows below.
latex_balance_header <- function() {
  c("\\toprule",
   paste(" & \\multicolumn{2}{c}{Untreated} & \\multicolumn{2}{c}{Treated}",
         "& \\multicolumn{1}{c}{Std. Mean Diff.}",
         "& \\multicolumn{1}{c}{Var. Ratio Dev.} \\\\"),
   "\\cmidrule(lr){2-3} \\cmidrule(lr){4-5} \\cmidrule(lr){6-6} \\cmidrule(lr){7-7}",
   paste0(" & \\multicolumn{1}{c}{", mu0.sym, "} & \\multicolumn{1}{c}{", sigma0.sym, "}",
         "& \\multicolumn{1}{c}{", mu1.sym, "} & \\multicolumn{1}{c}{", sigma1.sym, "}",
         "& \\multicolumn{1}{c}{", smd.formula, "}",
         "& \\multicolumn{1}{c}{", vrd.formula, "} \\\\"),
   "\\midrule")
}

# HTML header rows: mirrors latex_balance_header(), using <th colspan> for
# the group labels and plain HTML entities/CSS (see above) for the
# symbols/formulas -- no MathJax dependency
html_balance_header <- function() {
  c("<thead>",
   "<tr>",
   "<th></th>",
   "<th colspan=\"2\" style=\"text-align:center; border-bottom:1px solid black;\">Untreated</th>",
   "<th colspan=\"2\" style=\"text-align:center; border-bottom:1px solid black;\">Treated</th>",
   "<th style=\"text-align:center; border-bottom:1px solid black;\">Std. Mean Diff.</th>",
   "<th style=\"text-align:center; border-bottom:1px solid black;\">Var. Ratio Dev.</th>",
   "</tr>",
   "<tr style=\"border-bottom:2px solid black;\">",
   "<th></th>",
   paste0("<th style=\"text-align:center;\">", mu0.html, "</th>"),
   paste0("<th style=\"text-align:center;\">", sigma0.html, "</th>"),
   paste0("<th style=\"text-align:center;\">", mu1.html, "</th>"),
   paste0("<th style=\"text-align:center;\">", sigma1.html, "</th>"),
   paste0("<th style=\"text-align:center;\">", smd.formula.html, "</th>"),
   paste0("<th style=\"text-align:center;\">", vrd.formula.html, "</th>"),
   "</tr>",
   "</thead>")
}

# resolves group_breaks into integer row indices marking the end of each
# group. Accepts either raw (1-indexed) row indices, for backward
# compatibility, or a list of covariate-name vectors defining each group in
# order -- the latter is robust to reordering or adding covariates, since it
# does not depend on hardcoded row positions. When named groups are given,
# they must account for every covariate in the table, in table order, so
# mistakes (typos, missing/reordered covariates) are caught immediately
# instead of silently drawing rules in the wrong place.
resolve_group_breaks <- function(group_breaks, covariate.names) {

  if (is.null(group_breaks) || is.numeric(group_breaks)) {
    return(group_breaks)
  }

  if (!is.list(group_breaks)) {
    stop("group_breaks must be NULL, an integer vector of row indices, ",
        "or a list of covariate-name vectors.")
  }

  flat <- unlist(group_breaks, use.names = FALSE)
  if (!identical(flat, covariate.names)) {
    stop("group_breaks must list every covariate in 'x' exactly once, ",
        "in the same order as the balance table.")
  }

  cumsum(vapply(group_breaks, length, integer(1)))
}

# one formatted LaTeX row per covariate, with an optional rule after the
# given (1-indexed) rows to separate covariate groups (e.g., region dummies
# vs. economic variables, as in Table 2)
latex_balance_rows <- function(table, digits, group_breaks = NULL) {
  covariate.names <- escape_latex(rownames(table))
  values <- format_latex_num(table, digits = digits)
  dim(values) <- dim(table)   # formatC() does not preserve matrix shape
  colnames(values) <- colnames(table)

  s.cols <- c("std.mean.diff", "var.ratio.dev")
  values[, s.cols] <- wrap_if_placeholder(values[, s.cols, drop = FALSE])

  rows <- apply(values, 1, paste, collapse = " & ")
  rows <- paste0(covariate.names, " & ", rows, " \\\\")

  if (!is.null(group_breaks)) {
    # a break after the last row would duplicate \bottomrule, so drop it
    group_breaks <- group_breaks[group_breaks < nrow(table)]
    rows[group_breaks] <- paste0(rows[group_breaks], "\n\\midrule")
  }

  rows
}

# one formatted HTML row (<tr>...</tr>) per covariate. Rather than adding a
# rule after the last row of a group (which, at the very last row, would sit
# oddly against the table's own border), a top border is added to the *first*
# row of each new group -- this also sidesteps the "duplicate rule" issue
# latex_balance_rows() has to explicitly guard against.
html_balance_rows <- function(table, digits, group_breaks = NULL) {
  covariate.names <- escape_html(rownames(table))
  values <- format_html_num(table, digits = digits)
  dim(values) <- dim(table)
  colnames(values) <- colnames(table)

  align <- c(mean.control = "right", var.control = "right",
            mean.treated = "right", var.treated = "right",
            std.mean.diff = "center", var.ratio.dev = "center")

  group.start <- if (is.null(group_breaks)) integer(0) else
    group_breaks[group_breaks < nrow(table)] + 1

  rows <- vapply(seq_len(nrow(table)), function(i) {
    border <- if (i %in% group.start) " style=\"border-top:1px solid black;\"" else ""
    cells <- paste0("<td style=\"text-align:", align, ";\">", values[i, ], "</td>",
                    collapse = "")
    paste0("<tr", border, "><td>", covariate.names[i], "</td>", cells, "</tr>")
  }, character(1))

  c("<tbody>", rows, "</tbody>")
}

# assembles the full LaTeX \begin{table}...\end{table} source
render_latex_balance <- function(x, digits = 3, caption = NULL, label = NULL,
                                 group_breaks = NULL) {

  if (is.null(caption)) {
    caption <- paste0("This table reports covariate balance between ",
                      x$info$n.treated, " treated and ", x$info$n.control,
                      " control units. From left to right, the columns report ",
                      "group means, variances, standardized mean differences, ",
                      "and deviations of the variance ratio from one.")
  }
  caption <- escape_latex(caption)
  if (!is.null(label)) label <- escape_latex(label)

  smd.spec <- siunitx_format(x$table[, "std.mean.diff"], digits)
  vrd.spec <- siunitx_format(x$table[, "var.ratio.dev"], digits)

  c("\\begin{table}[ht]",
   "\\centering",
   paste0("\\begin{tabular}{lrrrr", smd.spec, vrd.spec, "}"),
   latex_balance_header(),
   latex_balance_rows(x$table, digits = digits, group_breaks = group_breaks),
   "\\bottomrule",
   "\\end{tabular}",
   paste0("\\caption{", caption, "}"),
   if (!is.null(label)) paste0("\\label{", label, "}"),
   "\\end{table}")
}

# assembles the full HTML <table>...</table> source
render_html_balance <- function(x, digits = 3, caption = NULL, label = NULL,
                                group_breaks = NULL) {

  if (is.null(caption)) {
    caption <- paste0("This table reports covariate balance between ",
                      x$info$n.treated, " treated and ", x$info$n.control,
                      " control units. From left to right, the columns report ",
                      "group means, variances, standardized mean differences, ",
                      "and deviations of the variance ratio from one.")
  }
  caption <- escape_html(caption)

  # `label` becomes an id="..." attribute, HTML's equivalent of a LaTeX
  # \label -- it lets other content link to the table (e.g. <a href="#...">)
  id.attr <- if (!is.null(label)) paste0(" id=\"", escape_html(label), "\"") else ""

  c("<style>table.dml-balance-table td, table.dml-balance-table th { padding: 4px 10px; }</style>",
   paste0("<table class=\"dml-balance-table\"", id.attr,
         " style=\"border-collapse: collapse; margin: 0 auto;\">"),
   html_balance_header(),
   html_balance_rows(x$table, digits = digits, group_breaks = group_breaks),
   "</table>",
   paste0("<p style=\"text-align:center;\">", caption, "</p>"))
}

##' Formatted table for covariate balance
##' @description
##' Formats an object of class \code{\link{dml_balance}} as a LaTeX or HTML table, replicating the layout of Table 2 in Chernozhukov, Cinelli, Newey, Sharma and Syrgkanis (2026). Intended to be used in an R Markdown/knitr chunk with \code{results = "asis"}. For \code{format = "latex"}, the document's LaTeX preamble must include \code{\\usepackage{booktabs}} (for the table rules) and \code{\\usepackage{siunitx}} (for decimal-point alignment of the Std. Mean Diff./Var. Ratio Dev. columns). For \code{format = "html"}, the output is plain, self-contained HTML/CSS with no MathJax or other JavaScript dependency, so it renders correctly even when the file is opened offline or directly from disk.
##' @param x an object of class \code{\link{dml_balance}}.
##' @param ... arguments passed to other methods.
##' @returns A character vector with the table source, printed to the console (or knitted document) and returned invisibly.
##' @export
covariate_balance_table <- function(x, ...) {
  UseMethod("covariate_balance_table")
}

##' @param format output format, either \code{"latex"} or \code{"html"}. Default is \code{"latex"}.
##' @param digits number of decimal places to display. Default is \code{3}.
##' @param caption table caption. Default (\code{NULL}) generates a caption analogous to Table 2 of Chernozhukov et al (2026).
##' @param label a reference identifier for the table: for \code{format = "latex"}, this becomes a \code{\\label{}} for cross-referencing with \code{\\ref{}}; for \code{format = "html"}, this becomes an \code{id="..."} attribute on the \code{<table>}, so other content can link to it (e.g. \code{<a href="#label">}). Default is \code{NULL} (no label/id).
##' @param group_breaks either an integer vector with the (1-indexed) covariate rows after which a rule should be drawn, or a list of character vectors giving the covariate names in each group (which must together account for every covariate in \code{x}, in table order). Used to visually separate groups of covariates (e.g., region dummies vs. economic variables). Default (\code{NULL}) draws no such rules.
##' @rdname covariate_balance_table
##' @export
covariate_balance_table.dml_balance <- function(x, format = c("latex", "html"),
                                                digits = 3, caption = NULL, label = NULL,
                                                group_breaks = NULL, ...) {

  format <- match.arg(format)
  group_breaks <- resolve_group_breaks(group_breaks, rownames(x$table))

  out <- if (format == "latex") {
    render_latex_balance(x, digits = digits, caption = caption, label = label,
                         group_breaks = group_breaks)
  } else {
    render_html_balance(x, digits = digits, caption = caption, label = label,
                        group_breaks = group_breaks)
  }

  cat(out, sep = "\n")
  invisible(out)
}
