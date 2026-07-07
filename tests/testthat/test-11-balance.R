# Test covariate balance table (dml_balance)
#
# Includes a replication of Table 2 of Chernozhukov, Cinelli, Newey, Sharma
# and Syrgkanis (2026), which reports covariate balance between treated and
# control counties in the minimum-wage application of Callaway and Sant'Anna
# (2021). The raw county-level data is not bundled with this package, so the
# replication below uses the group means and variances exactly as reported
# in the manuscript's Table 2 and checks that our formulas reproduce the
# published standardized mean differences and variance-ratio deviations.

library(testthat)
library(dml.sensemakr)

# === covariate_group_stats ===
test_that("covariate_group_stats computes group means and variances", {
  xj <- c(1, 2, 3, 10, 20, 30)
  d  <- c(0, 0, 0, 1, 1, 1)
  stats <- dml.sensemakr:::covariate_group_stats(xj, d)
  expect_equal(stats$mean.control, mean(c(1, 2, 3)))
  expect_equal(stats$var.control, var(c(1, 2, 3)))
  expect_equal(stats$mean.treated, mean(c(10, 20, 30)))
  expect_equal(stats$var.treated, var(c(10, 20, 30)))
})

test_that("covariate_group_stats handles unbalanced group sizes", {
  xj <- c(5, 5, 5, 1, 9)
  d  <- c(0, 0, 0, 1, 1)
  stats <- dml.sensemakr:::covariate_group_stats(xj, d)
  expect_equal(stats$mean.control, 5)
  expect_equal(stats$var.control, 0)
  expect_equal(stats$mean.treated, 5)
  expect_equal(stats$var.treated, var(c(1, 9)))
})

# === std_mean_diff ===
test_that("std_mean_diff computes standardized mean difference", {
  result <- dml.sensemakr:::std_mean_diff(mean.control = 1, mean.treated = 3,
                                          var.control = 1, var.treated = 3)
  expected <- (3 - 1) / sqrt((3 + 1) / 2)
  expect_equal(result, expected)
})

test_that("std_mean_diff is zero when group means are equal", {
  expect_equal(dml.sensemakr:::std_mean_diff(1, 1, 2, 4), 0)
})

test_that("std_mean_diff is negative when treated mean is lower", {
  result <- dml.sensemakr:::std_mean_diff(mean.control = 5, mean.treated = 2,
                                          var.control = 1, var.treated = 1)
  expect_true(result < 0)
})

# === var_ratio_dev ===
test_that("var_ratio_dev computes deviation of variance ratio from one", {
  result <- dml.sensemakr:::var_ratio_dev(var.control = 2, var.treated = 4)
  expect_equal(result, 4 / 2 - 1)
})

test_that("var_ratio_dev is zero when variances are equal", {
  expect_equal(dml.sensemakr:::var_ratio_dev(5, 5), 0)
})

test_that("var_ratio_dev is negative when treated variance is smaller", {
  result <- dml.sensemakr:::var_ratio_dev(var.control = 4, var.treated = 1)
  expect_true(result < 0)
})

# === replication of Table 2 (Chernozhukov et al., 2026) ===
# Group means/variances for treated (2007 minimum-wage-increase cohort, n =
# 584) and control (never-treated, n = 1,377) counties, transcribed from
# Table 2 of the manuscript.
table2 <- data.frame(
  covariate     = c("Midwest", "South", "West", "Population (1000s)",
                    "log(Population)", "Median Inc. (1000s)",
                    "log(Median Inc.)"),
  mean.control  = c(0.336, 0.593, 0.072, 53.425, 3.018, 31.889, 3.438),
  var.control   = c(0.223, 0.241, 0.067, 23438,  1.584, 54.565, 0.048),
  mean.treated  = c(0.483, 0.301, 0.216, 84.142, 3.467, 33.080, 3.471),
  var.treated   = c(0.250, 0.212, 0.169, 32850,  1.781, 66.654, 0.053),
  std.mean.diff = c(0.303, -0.613, 0.419, 0.183, 0.346, 0.153, 0.150),
  var.ratio.dev = c(0.120, -0.128, 1.536, 0.402, 0.125, 0.222, 0.111),
  stringsAsFactors = FALSE
)
# Note: White, HS Graduates, and Poverty Rate are intentionally left out of
# this replication. Their reported variances are so small (0.003-0.026) that
# rounding them to 3 decimals in the manuscript materially changes the
# implied variance ratio, so recomputing from the rounded table entries no
# longer matches the published deviation closely enough to be a meaningful
# check.

test_that("std_mean_diff replicates Table 2 of Chernozhukov et al. (2026)", {
  computed <- mapply(dml.sensemakr:::std_mean_diff,
                     mean.control = table2$mean.control,
                     mean.treated = table2$mean.treated,
                     var.control  = table2$var.control,
                     var.treated  = table2$var.treated)
  expect_equal(unname(computed), table2$std.mean.diff, tolerance = 0.01)
})

test_that("var_ratio_dev replicates Table 2 of Chernozhukov et al. (2026)", {
  computed <- mapply(dml.sensemakr:::var_ratio_dev,
                     var.control = table2$var.control,
                     var.treated = table2$var.treated)
  expect_equal(unname(computed), table2$var.ratio.dev, tolerance = 0.02)
})

# === input checks ===
# A lightweight stand-in for a fitted dml object: dml_balance() only reads
# model$data$x and model$data$d, so we don't need to run a full dml() fit
# to exercise its input checks.
fake_x <- cbind(age = c(20, 25, 30, 40, 45, 50),
                income = c(10, 12, 15, 20, 22, 25))
fake_d_binary <- c(0, 0, 0, 1, 1, 1)
fake_d_continuous <- c(1.2, 3.4, 5.6, 7.8, 9.0, 1.1)

new_fake_dml <- function(x, d) {
  structure(list(data = list(x = x, d = d)), class = "dml")
}

test_that("dml_balance errors when model is not a dml object", {
  expect_error(dml_balance(list(data = list(x = fake_x, d = fake_d_binary))),
              "class 'dml'")
})

test_that("dml_balance errors when model is missing data", {
  bad_model <- structure(list(), class = "dml")
  expect_error(dml_balance(bad_model), "covariates")
})

test_that("dml_balance errors when treatment is not binary", {
  model <- new_fake_dml(fake_x, fake_d_continuous)
  expect_error(dml_balance(model), "binary")
})

test_that("dml_balance errors when requested covariates are not found", {
  model <- new_fake_dml(fake_x, fake_d_binary)
  expect_error(dml_balance(model, covariates = "not_a_covariate"),
              "not_a_covariate")
})

test_that("dml_balance defaults to all covariates when covariates is NULL", {
  model <- new_fake_dml(fake_x, fake_d_binary)
  bal <- dml_balance(model)
  expect_equal(bal$info$covariates, colnames(fake_x))
})

test_that("dml_balance preserves the requested covariate order, not colnames(x) order", {
  # colnames(fake_x) is c("age", "income"); group_breaks (see below) relies
  # on the table's row order matching the requested covariate order exactly
  model <- new_fake_dml(fake_x, fake_d_binary)
  bal <- dml_balance(model, covariates = c("income", "age"))
  expect_equal(rownames(bal$table), c("income", "age"))
  expect_equal(bal$info$covariates, c("income", "age"))
})

test_that("dml_balance degrades to NA/NaN rather than erroring when a treatment group is empty", {
  # d has no variation at all (still trivially passes the binary check), so
  # the treated group is empty -- mean()/var() of an empty vector return
  # NaN/NA in base R rather than erroring, and that should propagate here
  # rather than crashing.
  x <- cbind(age = c(20, 25, 30, 40, 45, 50))
  d <- rep(0, 6)
  model <- new_fake_dml(x, d)
  bal <- dml_balance(model)
  expect_true(is.na(bal$table["age", "mean.treated"]))
  expect_true(is.na(bal$table["age", "var.treated"]))
})

# === dml_balance() return structure ===
test_that("dml_balance returns correct class and structure", {
  model <- new_fake_dml(fake_x, fake_d_binary)
  bal <- dml_balance(model, covariates = c("age", "income"))

  expect_s3_class(bal, "dml_balance")
  expect_true(all(c("info", "table") %in% names(bal)))
  expect_equal(bal$info$covariates, c("age", "income"))
  expect_equal(bal$info$n.control, sum(fake_d_binary == 0))
  expect_equal(bal$info$n.treated, sum(fake_d_binary == 1))

  expect_true(is.matrix(bal$table))
  expect_equal(rownames(bal$table), c("age", "income"))
  expect_equal(colnames(bal$table),
              c("mean.control", "var.control", "mean.treated", "var.treated",
                "std.mean.diff", "var.ratio.dev"))
})

test_that("dml_balance table values match direct computation", {
  model <- new_fake_dml(fake_x, fake_d_binary)
  bal <- dml_balance(model, covariates = "age")

  age.control <- fake_x[fake_d_binary == 0, "age"]
  age.treated <- fake_x[fake_d_binary == 1, "age"]

  expect_equal(bal$table["age", "mean.control"], mean(age.control))
  expect_equal(bal$table["age", "var.control"], var(age.control))
  expect_equal(bal$table["age", "mean.treated"], mean(age.treated))
  expect_equal(bal$table["age", "var.treated"], var(age.treated))
  expect_equal(bal$table["age", "std.mean.diff"],
              (mean(age.treated) - mean(age.control)) /
                sqrt((var(age.treated) + var(age.control)) / 2))
  expect_equal(bal$table["age", "var.ratio.dev"],
              var(age.treated) / var(age.control) - 1)
})

# === dml_balance() integration with a real dml() fit ===
test_that("dml_balance works with an actual dml() fit", {
  data("pension", package = "dml.sensemakr")
  set.seed(11)
  idx <- sample(nrow(pension), 300)  # small subset for speed
  y <- pension$net_tfa[idx]
  d <- pension$e401[idx]
  x <- model.matrix(~ -1 + age + inc + educ + fsize + marr + twoearn + pira + hown,
                     data = pension[idx, ])
  fit <- dml(y, d, x, model = "plm", cf.folds = 2, cf.reps = 1, verbose = FALSE)

  bal <- dml_balance(fit, covariates = c("age", "inc"))
  expect_s3_class(bal, "dml_balance")
  expect_equal(bal$info$n.control, sum(d == 0))
  expect_equal(bal$info$n.treated, sum(d == 1))
  expect_equal(bal$table["age", "mean.control"], mean(x[d == 0, "age"]))
  expect_equal(bal$table["inc", "mean.treated"], mean(x[d == 1, "inc"]))
})

# === print and summary methods ===
test_that("summary.dml_balance returns object of correct class", {
  model <- new_fake_dml(fake_x, fake_d_binary)
  bal <- dml_balance(model)
  s <- summary(bal)
  expect_s3_class(s, "summary_dml_balance")
})

test_that("print.dml_balance does not error", {
  model <- new_fake_dml(fake_x, fake_d_binary)
  bal <- dml_balance(model)
  expect_output(print(bal), regexp = NULL)
})

test_that("print.summary_dml_balance does not error", {
  model <- new_fake_dml(fake_x, fake_d_binary)
  bal <- dml_balance(model)
  expect_output(print(summary(bal)), regexp = NULL)
})

# === format_latex_value ===
test_that("format_latex_value formats plain numbers with fixed decimals", {
  expect_equal(dml.sensemakr:::format_latex_value(0.1234, digits = 3), "0.123")
  expect_equal(dml.sensemakr:::format_latex_value(-0.1, digits = 3), "-0.100")
})

test_that("format_latex_value adds a thousands separator for large magnitudes", {
  expect_equal(dml.sensemakr:::format_latex_value(23438, digits = 3), "23,438")
})

test_that("format_latex_value handles NA/NaN and +-Inf as placeholders", {
  expect_equal(dml.sensemakr:::format_latex_value(NA_real_, digits = 3), "--")
  expect_equal(dml.sensemakr:::format_latex_value(NaN, digits = 3), "--")
  expect_equal(dml.sensemakr:::format_latex_value(Inf, digits = 3), "$\\infty$")
  expect_equal(dml.sensemakr:::format_latex_value(-Inf, digits = 3), "$-\\infty$")
})

# === wrap_if_placeholder ===
test_that("wrap_if_placeholder wraps only non-numeric placeholders", {
  values <- c("0.302", "--", "$\\infty$", "$-\\infty$", "-0.120")
  wrapped <- dml.sensemakr:::wrap_if_placeholder(values)
  expect_equal(wrapped[1], "0.302")
  expect_equal(wrapped[2], "\\multicolumn{1}{c}{--}")
  expect_equal(wrapped[3], "\\multicolumn{1}{c}{$\\infty$}")
  expect_equal(wrapped[4], "\\multicolumn{1}{c}{$-\\infty$}")
  expect_equal(wrapped[5], "-0.120")
})

# === escape_latex / escape_html ===
# These matter beyond cosmetics: covariate names from real data (e.g.
# model.matrix() output) very commonly contain underscores, and an
# unescaped "_" in LaTeX text mode is a compile error, not just a display
# glitch.
test_that("escape_latex escapes special characters, most importantly underscores", {
  expect_equal(dml.sensemakr:::escape_latex("income_2020"), "income\\_2020")
  expect_equal(dml.sensemakr:::escape_latex("100% & #1 {x}"), "100\\% \\& \\#1 \\{x\\}")
})

test_that("escape_html escapes HTML special characters", {
  expect_equal(dml.sensemakr:::escape_html("A & B"), "A &amp; B")
  expect_equal(dml.sensemakr:::escape_html("<tag>"), "&lt;tag&gt;")
  expect_equal(dml.sensemakr:::escape_html('say "hi"'), "say &quot;hi&quot;")
})

test_that("covariate_balance_table escapes special characters in covariate (row) names, not just captions", {
  x <- cbind(income_2020 = c(20, 25, 30, 40, 45, 50))
  d <- c(0, 0, 0, 1, 1, 1)
  model <- new_fake_dml(x, d)
  bal <- dml_balance(model)

  latex_out <- capture.output(covariate_balance_table(bal, format = "latex"))
  expect_true(any(grepl("income\\\\_2020", latex_out)))

  html_out <- capture.output(covariate_balance_table(bal, format = "html"))
  expect_true(any(grepl(">income_2020<", html_out)))
})

# === siunitx_format ===
test_that("siunitx_format scales the integer-digit width to the data", {
  expect_equal(dml.sensemakr:::siunitx_format(c(0.1, -0.5, 0.9), digits = 3),
              "S[table-format=-1.3]")
  expect_equal(dml.sensemakr:::siunitx_format(c(12.3, -1.5), digits = 2),
              "S[table-format=-2.2]")
})

test_that("siunitx_format ignores non-finite values when sizing the column", {
  expect_equal(dml.sensemakr:::siunitx_format(c(0.5, NA, Inf), digits = 3),
              "S[table-format=-1.3]")
})

test_that("siunitx_format falls back to 1 integer digit when no finite values exist", {
  expect_equal(dml.sensemakr:::siunitx_format(c(NA, Inf, -Inf), digits = 3),
              "S[table-format=-1.3]")
})

# === resolve_group_breaks ===
test_that("resolve_group_breaks passes through NULL and numeric input unchanged", {
  expect_null(dml.sensemakr:::resolve_group_breaks(NULL, c("a", "b")))
  expect_equal(dml.sensemakr:::resolve_group_breaks(c(1, 3), c("a", "b", "c", "d")), c(1, 3))
})

test_that("resolve_group_breaks converts a valid named-group list to row indices", {
  breaks <- dml.sensemakr:::resolve_group_breaks(
    list(c("age", "inc"), c("educ")), c("age", "inc", "educ"))
  expect_equal(breaks, c(2, 3))
})

test_that("resolve_group_breaks errors when the groups omit or reorder covariates", {
  expect_error(
    dml.sensemakr:::resolve_group_breaks(list(c("age"), c("inc")), c("age", "inc", "educ")),
    "exactly once"
  )
  expect_error(
    dml.sensemakr:::resolve_group_breaks(list(c("inc", "age")), c("age", "inc")),
    "exactly once"
  )
})

test_that("resolve_group_breaks errors on an invalid type", {
  expect_error(dml.sensemakr:::resolve_group_breaks("age", c("age", "inc")),
              "must be NULL")
})

# === covariate_balance_table.dml_balance (LaTeX) ===
test_that("covariate_balance_table defaults to latex format", {
  model <- new_fake_dml(fake_x, fake_d_binary)
  bal <- dml_balance(model)
  out <- capture.output(result <- covariate_balance_table(bal))
  expect_true(any(grepl("\\\\begin\\{table\\}", out)))
  expect_true(any(grepl("\\\\begin\\{tabular\\}", out)))
  expect_equal(out, result)
})

test_that("covariate_balance_table(format = 'latex') uses siunitx S columns", {
  model <- new_fake_dml(fake_x, fake_d_binary)
  bal <- dml_balance(model)
  out <- capture.output(covariate_balance_table(bal, format = "latex"))
  expect_true(any(grepl("S\\[table-format=", out)))
})

test_that("covariate_balance_table(format = 'latex') escapes special characters in a caption", {
  model <- new_fake_dml(fake_x, fake_d_binary)
  bal <- dml_balance(model)
  out <- capture.output(covariate_balance_table(bal, format = "latex",
                                                caption = "100% balanced & tested"))
  expect_true(any(grepl("100\\\\% balanced \\\\& tested", out)))
})

test_that("covariate_balance_table(format = 'latex') draws a rule at named group boundaries", {
  model <- new_fake_dml(fake_x, fake_d_binary)
  bal <- dml_balance(model)
  out <- capture.output(covariate_balance_table(bal, format = "latex",
                                                group_breaks = list("age", "income")))
  # a \midrule should appear after the "age" row and before the closing rules
  age_line <- which(grepl("^age &", out))
  expect_true(any(grepl("midrule", out[age_line + 1])))
})

test_that("covariate_balance_table(format = 'latex') does not draw a duplicate rule before \\bottomrule", {
  # group_breaks that fully partitions the table has its last boundary land
  # on the final row -- that should not produce \midrule immediately
  # followed by \bottomrule
  model <- new_fake_dml(fake_x, fake_d_binary)
  bal <- dml_balance(model)
  out <- capture.output(covariate_balance_table(bal, format = "latex",
                                                group_breaks = list("age", "income")))
  bottomrule_line <- which(grepl("\\\\bottomrule", out))
  expect_false(grepl("midrule", out[bottomrule_line - 1]))
})

test_that("covariate_balance_table(format = 'latex') renders NA/Inf gracefully instead of erroring", {
  # zero variance in the control group makes var.ratio.dev = Inf
  x <- cbind(age = c(30, 30, 30, 40, 45, 50))
  d <- c(0, 0, 0, 1, 1, 1)
  model <- new_fake_dml(x, d)
  bal <- dml_balance(model)
  expect_true(is.infinite(bal$table["age", "var.ratio.dev"]))
  out <- capture.output(covariate_balance_table(bal, format = "latex"))
  expect_true(any(grepl("\\\\infty", out)))
})

# === covariate_balance_table.dml_balance (HTML) ===
test_that("covariate_balance_table(format = 'html') returns an HTML table", {
  model <- new_fake_dml(fake_x, fake_d_binary)
  bal <- dml_balance(model)
  out <- capture.output(result <- covariate_balance_table(bal, format = "html"))
  expect_true(any(grepl("<table", out)))
  expect_true(any(grepl("</table>", out)))
  expect_equal(out, result)
})

test_that("covariate_balance_table(format = 'html') does not depend on MathJax/JS", {
  # header symbols/formulas must be plain HTML entities, not $...$ math that
  # would require MathJax (which fails to load offline or from file://)
  model <- new_fake_dml(fake_x, fake_d_binary)
  bal <- dml_balance(model)
  out <- capture.output(covariate_balance_table(bal, format = "html"))
  expect_false(any(grepl("\\$", out)))
  expect_false(any(grepl("<script", out)))
  expect_true(any(grepl("&mu;", out)))
  expect_true(any(grepl("&sigma;", out)))
})

test_that("covariate_balance_table(format = 'html') adds cell padding so values don't visually run together", {
  model <- new_fake_dml(fake_x, fake_d_binary)
  bal <- dml_balance(model)
  out <- capture.output(covariate_balance_table(bal, format = "html"))
  expect_true(any(grepl("padding", out)))
})

test_that("covariate_balance_table(format = 'html') turns label into an id attribute", {
  model <- new_fake_dml(fake_x, fake_d_binary)
  bal <- dml_balance(model)
  out <- capture.output(covariate_balance_table(bal, format = "html", label = "tab:balance"))
  expect_true(any(grepl("id=\"tab:balance\"", out)))
})

test_that("covariate_balance_table(format = 'html') omits id attribute when label is NULL", {
  model <- new_fake_dml(fake_x, fake_d_binary)
  bal <- dml_balance(model)
  out <- capture.output(covariate_balance_table(bal, format = "html"))
  expect_false(any(grepl("id=", out)))
})

test_that("covariate_balance_table(format = 'html') escapes quotes in the label", {
  model <- new_fake_dml(fake_x, fake_d_binary)
  bal <- dml_balance(model)
  out <- capture.output(covariate_balance_table(bal, format = "html", label = 'bad"label'))
  expect_true(any(grepl("id=\"bad&quot;label\"", out)))
})

test_that("covariate_balance_table(format = 'html') escapes special characters in a caption", {
  model <- new_fake_dml(fake_x, fake_d_binary)
  bal <- dml_balance(model)
  out <- capture.output(covariate_balance_table(bal, format = "html",
                                                caption = "100% balanced & tested"))
  expect_true(any(grepl("100% balanced &amp; tested", out)))
})

test_that("covariate_balance_table(format = 'html') draws a border at named group boundaries", {
  model <- new_fake_dml(fake_x, fake_d_binary)
  bal <- dml_balance(model)
  out <- capture.output(covariate_balance_table(bal, format = "html",
                                                group_breaks = list("age", "income")))
  income_line <- which(grepl(">income<", out))
  expect_true(any(grepl("border-top", out[income_line])))
})

test_that("covariate_balance_table(format = 'html') renders NA/Inf gracefully instead of erroring", {
  x <- cbind(age = c(30, 30, 30, 40, 45, 50))
  d <- c(0, 0, 0, 1, 1, 1)
  model <- new_fake_dml(x, d)
  bal <- dml_balance(model)
  out <- capture.output(covariate_balance_table(bal, format = "html"))
  expect_true(any(grepl("infin;", out)))
})

test_that("covariate_balance_table errors on an invalid format", {
  model <- new_fake_dml(fake_x, fake_d_binary)
  bal <- dml_balance(model)
  expect_error(covariate_balance_table(bal, format = "docx"))
})

# === integration: does the LaTeX output actually compile? ===
# The tests above only check for specific substrings (e.g. "\midrule",
# "S[table-format=") -- none of them confirm the emitted LaTeX is valid
# syntax that a real LaTeX engine accepts. This test wraps the table in a
# minimal standalone document and compiles it with pdflatex, skipping
# gracefully wherever no LaTeX engine is installed (e.g. most CI runners).
# Wraps a dml_balance object's LaTeX output in a minimal standalone document
# and compiles it with pdflatex, returning TRUE/FALSE for whether a PDF was
# produced (used by the two compile tests below).
compiles_with_pdflatex <- function(bal, pdflatex, ...) {
  table_tex <- capture.output(covariate_balance_table(bal, format = "latex", ...))

  doc <- c("\\documentclass{article}",
          "\\usepackage{booktabs}",
          "\\usepackage{siunitx}",
          "\\begin{document}",
          table_tex,
          "\\end{document}")

  tmp_dir <- tempfile("dml_balance_latex_test")
  dir.create(tmp_dir)
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

  tex_file <- file.path(tmp_dir, "balance_test.tex")
  writeLines(doc, tex_file)

  old_wd <- setwd(tmp_dir)
  on.exit(setwd(old_wd), add = TRUE)

  result <- suppressWarnings(system2(pdflatex,
                                     c("-interaction=nonstopmode", "-halt-on-error",
                                       basename(tex_file)),
                                     stdout = TRUE, stderr = TRUE))
  status <- attr(result, "status")
  pdf_produced <- file.exists(file.path(tmp_dir, "balance_test.pdf"))

  if (!pdf_produced) {
    cat(result, sep = "\n")  # surface the LaTeX log on failure for debugging
  }
  list(pdf_produced = pdf_produced, status = status)
}

test_that("covariate_balance_table(format = 'latex') output compiles with pdflatex", {
  pdflatex <- Sys.which("pdflatex")
  skip_if(!nzchar(pdflatex), "pdflatex not available")

  model <- new_fake_dml(fake_x, fake_d_binary)
  bal <- dml_balance(model)
  result <- compiles_with_pdflatex(bal, pdflatex, label = "tab:test")

  expect_true(result$pdf_produced)
  expect_true(is.null(result$status) || result$status == 0)
})

test_that("covariate_balance_table(format = 'latex') output with an underscore in a covariate name compiles", {
  # this is the realistic failure mode: model.matrix()-style covariate names
  # almost always contain underscores, and an unescaped "_" in LaTeX text
  # mode is a compile error, not just a display glitch -- string-matching
  # tests alone can't catch a regression in escape_latex() as reliably as
  # actually compiling the output.
  pdflatex <- Sys.which("pdflatex")
  skip_if(!nzchar(pdflatex), "pdflatex not available")

  x <- cbind(income_2020 = c(20, 25, 30, 40, 45, 50))
  d <- c(0, 0, 0, 1, 1, 1)
  model <- new_fake_dml(x, d)
  bal <- dml_balance(model)
  result <- compiles_with_pdflatex(bal, pdflatex)

  expect_true(result$pdf_produced)
  expect_true(is.null(result$status) || result$status == 0)
})
