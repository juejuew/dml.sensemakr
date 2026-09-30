# Tested with did 2.5.1 and 2.5.1.902, DRDID 1.3.0 and adapter 0.2.2.9004.
# Install the adapter from the updated local source/package before running.
data("mpdta", package = "did")

fit <- did::att_gt(
  yname = "lemp", tname = "year", idname = "countyreal",
  gname = "first.treat", xformla = ~ lpop, data = mpdta,
  panel = TRUE, est_method = "dr", control_group = "nevertreated",
  bstrap = FALSE, cband = FALSE, compute_inffunc = TRUE
)

m <- dml.sensemakr::dml_from_did(
  fit, group = 2004, time = 2007,
  dml_args = list(
    yreg = "lm", dreg = "glm", d.class = TRUE,
    cf.folds = 3, cf.reps = 3, cf.seed = 17, ps.trim = 0
  )
)

result <- dml.sensemakr::sensemakr(m, cf.y = 0.05, cf.d = 0.05)
print(result)
