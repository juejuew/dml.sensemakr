# Full-sample parametric sensitivity analysis of did::att_gt() cells.
# Requires did, DRDID and this version of dml.sensemakr.
library(dml.sensemakr)
set.seed(10615)
n <- 600L
units <- data.frame(id = seq_len(n), g = rep(c(0, 3, 4), each = n / 3),
  x = runif(n, -1, 1), z = runif(n, -1, 1), u = rnorm(n))
dat <- do.call(rbind, lapply(1:4, function(t) {
  transform(units, time = t,
    y = u + t * (0.2 + 0.4 * x + 0.3 * z) +
      0.8 * (g > 0 & t >= g) + rnorm(n))
}))
source_fit <- did::att_gt(yname = "y", tname = "time", idname = "id", gname = "g",
  xformla = ~ x + z, data = dat, panel = TRUE, est_method = "dr",
  bstrap = FALSE, cband = FALSE, compute_inffunc = TRUE)

m <- dml_from_did(source_fit, group = 3, time = 4, method = "parametric")
print(m)
cell <- m$adapter_info$cell
stopifnot(isTRUE(all.equal(unname(coef(m)), source_fit$att[cell$column])),
          isTRUE(all.equal(unname(se(m)), source_fit$se[cell$column])))

# The usual downstream interface, now using full-sample moment IFs.
confidence_bounds(m, cf.y = 0.05, cf.d = 0.05)
sensemakr(m, cf.y = 0.05, cf.d = 0.05)
b <- dml_benchmark(m, "z")
benchmark_bounds(m, b, kY = 0.5, kD = 0.5)

batch <- dml_from_did_cells(source_fit, method = "parametric")
print(batch)
did_cell_apply(batch, confidence_bounds, cf.y = 0.05, cf.d = 0.05)

# Inspect what was fitted and the nuisance effects included in inference.
m$parametric
m$fits[[1]]$preds
m$results$main$treat[[1]]$joint.covariance
