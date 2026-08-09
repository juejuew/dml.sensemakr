# quietly() muffles only a fixed list of already-catalogued, incidental
# warnings that real model fits on small/edge-case test data are known to
# produce (small cross-fitting samples pushing nu2.s negative, near-constant
# dummy columns in leave-one-out benchmark refits, expected-sampling-noise
# negative gain estimates, etc.). Any OTHER warning -- in particular a new,
# unexpected one that would indicate a real regression -- still propagates
# and prints normally. This is deliberately not a blanket suppressWarnings():
# that would also hide a genuinely new warning from the exact same call.
known_noise_patterns <- c(
  "shows no variation",
  "^NaNs produced$",
  "Unable to compute bounds.*nu\\^2 is negative",
  "treated unit\\(s\\) in this cell",
  "trimmed \\(propensity score",
  "gain\\.[YD] \\(.*is negative",
  "Dropped.*already treated",  # did::att_gt()'s own anticipation>0 notice
  "Contour plot could not be drawn"
)

quietly <- function(expr) {
  withCallingHandlers(
    expr,
    warning = function(w) {
      if (any(vapply(known_noise_patterns, grepl, x = conditionMessage(w), logical(1)))) {
        invokeRestart("muffleWarning")
      }
    }
  )
}
