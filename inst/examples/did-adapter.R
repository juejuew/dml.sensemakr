# Tested with did 2.5.1 / 2.5.1.902 and DRDID 1.3.0; acceptance is based on input checks.
library(dml.sensemakr)

set.seed(608)
n <- 360L
units <- data.frame(id = sample(seq_len(n) * 7 + 100),
                    g = rep(c(0, 5, 9), each = n / 3),
                    x = runif(n, -1, 1), z = runif(n, -1, 1), u = rnorm(n))
panel_data <- do.call(rbind, lapply(c(0, 2, 5, 9, 12), function(t) {
  d <- units
  d$time <- t
  d$x <- d$x + 0.04 * t * d$z
  d$y <- d$u + t / 8 + (0.3 + t / 7) * d$x + 0.7 * d$z +
    (d$g > 0 & t >= d$g) * (0.8 + t / 20) + rnorm(n)
  d
}))

source_fit <- did::att_gt(
  yname = "y", tname = "time", idname = "id", gname = "g",
  xformla = ~ x + z, data = panel_data, panel = TRUE,
  est_method = "dr", control_group = "notyettreated",
  bstrap = FALSE, cband = FALSE, compute_inffunc = TRUE
)
args <- list(yreg = "lm", dreg = "glm", d.class = TRUE,
              cf.folds = 3, cf.reps = 3, cf.seed = 17, ps.trim = 0)

# One cell. No raw-data argument is needed.
hybrid <- dml_from_did(source_fit, group = 5, time = 9, dml_args = args)
summary(hybrid)
confidence_bounds(hybrid, cf.y = 0.05, cf.d = 0.05)
sensitivity <- sensemakr(hybrid, cf.y = 0.05, cf.d = 0.05)
benchmark <- dml_benchmark(hybrid, "z")
benchmark_bounds(hybrid, benchmark)

# All post-treatment cells. Inspect status before using individual results.
batch <- dml_from_did_cells(source_fit, dml_args = args)
print(batch)
intervals <- did_cell_apply(batch, confidence_bounds, cf.y = 0.05, cf.d = 0.05)
print(intervals)
intervals$values[["g=5,t=9"]]

# Subsets preserve the requested order; errors keep a named NULL entry.
selected <- dml_from_did_cells(source_fit,
  cells = data.frame(group = c(5, 9), time = c(12, 12)), dml_args = args)
benchmarks <- did_cell_apply(selected, dml_benchmark, benchmark_covariates = "z")
benchmarks$status
