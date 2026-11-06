# Adapter methods and tuning

The public methods for `dml_from_did()`, `dml_from_did_cells()` and
`dml_from_drdid()` are:

| Method | Sensitivity estimation | Tuning |
| --- | --- | --- |
| `"dml"` (default) | Cross-fitted DML moments | Clean by default; dirty optional |
| `"replicate"` | Source full-sample logistic-MLE/OLS fits and nuisance-adjusted IFs | Not applicable |

On the equal-weight adapter path both preserve the source ATT/IF and verify it
by numerical replay. Replicate describes the nuisance-fitting procedure, not
an additional bootstrap or an assurance of correct model specification.
Its fixed-dimensional model conditions are described in `parametric-inference.md`.
Replicate also supports source control residual trimming under the correct-mean
and regular-boundary conditions in `replicate-trimming.md`; sensitivity moments
continue to use all units and untrimmed propensities. The DML adapter still
rejects active source trimming. Existing unequal-weight DML fallback behavior
is unchanged.

```r
m <- dml_from_did(fit, group = 2007, time = 2007,
                  method = "replicate")

m_dml <- dml_from_did(fit, group = 2007, time = 2007,
                      method = "dml",
                      dml_args = list(dirty.tuning = TRUE))
```

`dirty.tuning = FALSE` tunes independently inside each training fold. `TRUE`
chooses hyperparameters on the full cell sample first (only controls for the
conditional ATT outcome model). Models are then trained separately on each
training fold and predict the held-out fold. Full-sample tuning does not
replace these fits with in-sample predictions, but selection uses information
from evaluation folds. Inference therefore requires suitable nuisance-rate
and tuning-selection conditions; arbitrary searches are not automatically
covered by the clean-tuning argument.

The original learner specification and `dirty.tuning` are saved in
`refit_spec$dml_args`, and the tuning choice is in `info$dirty.tuning`.
Covariate benchmarks and automatically constructed pre-period fits inherit
this choice and tune on their own data. They use the same cross-fitting
partitions/repetitions as the original fit. Selected hyperparameters may
differ between the post, pre and reduced-covariate fits.

The old names `"crossfit"` and `"parametric"` remain accepted aliases. New
objects store `"dml"`/`"replicate"` in `refit_spec$method` and
`adapter_info$method`; old saved objects are normalized when used for
benchmarking. Historical storage fields such as `parametric` and
`adapter_info$mode = "parametric"` remain for compatibility. A missing method
in an older saved adapter continues to mean DML.

The replicate method requires an empty `dml_args` list, so it rejects
`dirty.tuning` and other DML options instead of ignoring them.
