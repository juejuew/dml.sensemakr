# Tested with DRDID 1.3.0; acceptance uses input and runtime-interface checks.
library(DRDID)
library(dml.sensemakr)

set.seed(421)
n <- 300
x1 <- runif(n, -1, 1)
x2 <- runif(n, -1, 1)
D <- rbinom(n, 1, plogis(0.1 + 0.5 * x1 + 0.6 * x2))
y0 <- x1 + rnorm(n)
y1 <- y0 + 0.7 * D + x1 + 1.5 * x2 + rnorm(n)
panel_data <- data.frame(
  id = rep(seq_len(n), 2), time = rep(0:1, each = n),
  y = c(y0, y1), D = rep(D, 2), x1 = rep(x1, 2), x2 = rep(x2, 2)
)

source_fit <- DRDID::drdid(
  yname = "y", tname = "time", idname = "id", dname = "D",
  xformla = ~ x1 + x2, data = panel_data,
  panel = TRUE, estMethod = "imp", inffunc = TRUE
)

hybrid <- dml_from_drdid(
  source_fit, data = panel_data,
  dml_args = list(yreg = "lm", dreg = "glm", d.class = TRUE,
                  cf.folds = 3, cf.reps = 3, cf.seed = 17, ps.trim = 0)
)
print(hybrid)
stopifnot(isTRUE(all.equal(as.numeric(coef(hybrid)), as.numeric(source_fit$ATT))))

print(dml_bounds(hybrid, cf.y = 0.08, cf.d = 0.06))
print(confidence_bounds(hybrid, cf.y = 0.08, cf.d = 0.06))
print(robustness_value(hybrid))
print(sensemakr(hybrid, cf.y = 0.08, cf.d = 0.06))

# Both estimators are refit on the same panel units and cross-fitting partition.
benchmark <- dml_benchmark(hybrid, benchmark_covariates = "x2")
print(benchmark)
print(benchmark_bounds(hybrid, benchmark))

# Wide input: retain exactly the original source IF order and design matrix.
X <- model.matrix(~ x1 + x2)
source_wide <- DRDID::drdid_panel(y1, y0, D, covariates = X, inffunc = TRUE)
hybrid_wide <- dml_from_drdid(
  source_wide, y1 = y1, y0 = y0, D = D, covariates = X,
  unit_id = seq_len(n),
  dml_args = list(yreg = "lm", dreg = "glm", d.class = TRUE,
                  cf.folds = 3, cf.reps = 3, cf.seed = 17, ps.trim = 0)
)
print(hybrid_wide)
