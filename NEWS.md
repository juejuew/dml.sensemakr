# dml.sensemakr 0.2.2.9004

* Replace did/DRDID version equality gates with saved-input and runtime-interface
  checks. Require the same panel scope, source ATT/IF replay and overlap checks.
* Validate missing source metadata explicitly instead of guessing defaults.
  Check the arguments and return values of the DRDID functions actually used.
* Record installed versions and the selected design layout as diagnostics.
  A serialized did fit can be adapted without did installed if it retains the
  required data, settings, design and labeled influence functions.
* Validate both did 2.5.1 and 2.5.1.902 with DRDID 1.3.0, including three-repetition
  demos, cell reconstruction, corruption rejection, batching and benchmarking.

# dml.sensemakr 0.2.2.9003

* Both did and DRDID adapters now accept repeated cross-fitting through
  `dml_args$cf.reps` (default 1). Preserve every repetition's DML moments and
  combine them with the original ATT and IF using a separate joint covariance.
* Check every training partition and every repetition for invalid predictions,
  active trimming and invalid sensitivity moments. A failed repetition rejects
  the cell instead of being omitted from aggregation.
* Record all repetition seeds, fold assignments and auxiliary DML ATTs.
  Benchmark refits retain and verify all paired partitions after serialization.
* `sensemakr(..., combine.method = "mean")` now propagates the selected aggregation
  to robustness values, sensitivity intervals and benchmark tables; printing and
  default plots use the saved choice. The default remains `"median"`.
* Add repeated-estimation, inference, benchmark, batch and failure-path tests;
  retain the single-repetition behavior and all other adapter restrictions.

# dml.sensemakr 0.2.2.9002

* Add `dml_from_did()` for post-treatment group-time ATTs from `did::att_gt()`.
  Read the stored estimation sample and design, align units by IF row names,
  and undo did's global-sample IF scaling before estimating sensitivity moments.
* Support both computation modes, both control groups, anticipation and either
  base-period setting on the balanced panel estimator path. Require did
  2.5.1.902, DRDID 1.3.0, constant positive weights and inactive trimming.
* Add `dml_from_did_cells()` and `did_cell_apply()` with explicit per-cell error
  and warning records. Pre-treatment/reference cells remain listed as excluded.
* Refit both estimators on frozen cell samples and folds for benchmarks.
  Share the auxiliary DML/IF calculation with the existing DRDID adapter.
* Add executable examples and reconstruction, inference, batch, default-path
  and compatibility tests. Aggregated effects and simultaneous inference are
  outside this release's scope.

# dml.sensemakr 0.2.2.9001

* Benchmark group names must be unique, including automatically generated labels;
  duplicate names now fail before fitting instead of silently overwriting results.
* The DRDID adapter checks raw long-panel weights before normalization, enforcing
  constant positive values across units and periods, just as for wide input.
* Intercept-only adapter fits now correctly inherit the shared `reg` learner,
  honor separate `yreg`/`dreg` overrides, and avoid default range preprocessing.
* Expanded adapter regression tests cover factor/interaction designs, long/wide
  equivalence for both panel methods, grouped benchmarks and learner overrides,
  including independent checks of sensitivity moments and joint-IF standard errors.
* The Darfur test loads its data without attaching the separate `sensemakr`
  package, preventing its generics from masking those used by later tests.

# dml.sensemakr 0.2.2.9000

* Added `dml_from_drdid()` for hybrid conditional ATT sensitivity analysis from
  DRDID 1.3.0 panel results. The adapter preserves the external ATT and full IF,
  estimates sensitivity moments with cross-fitted DML, validates the common
  sample and supports joint-IF inference and hybrid benchmark refits. Initial
  support requires equal weights and no active propensity trimming/clipping.
* `sensemakr()` now passes its requested `alpha` to both manual and benchmark
  confidence-bound tables, keeping them consistent with robustness values.

# dml.sensemakr 0.2.2

* `confidence_bounds()` now returns the widest confidence envelope over all
  confounding no stronger than the supplied sensitivity values by default.
  Set `max = FALSE` to recover the previous fixed-corner calculation.
* Robustness values now use the exact confidence-envelope definition for mean
  and median aggregation, including interior extrema.

# dml.sensemakr 0.2.0

* Initial CRAN submission.
* Updated reference to Chernozhukov, Cinelli, Newey, Sharma, and Syrgkanis (2026), Review of Economics and Statistics.
* Bug fixes in `robustness_value()` for `dml.bounds` objects.
* Improved documentation throughout.
