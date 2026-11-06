# dml.sensemakr 0.2.2.9012

* Add `target = "atu"` to both did adapters with `method = "dml_refit"`.
  Keep ATT as the default and reject ATU for source-preserving methods.
* Fit the treated outcome arm and use the native conditional ATU estimates,
  all IFs, sensitivity scale and joint covariance; retain the original sign
  treatment minus no treatment. Validate the appropriate outcome-training arm.
* Preserve ATU through covariate benchmarks, repetitions, serialization, batch
  fitting, bounds and summaries. Record the selected control target population.
* State treatment-onset parallel trends and long-model sensitivity assumptions;
  ordinary ATT parallel trends alone is insufficient. Explicitly reject ATU
  pre-trend benchmarking pending a separate transport restriction.
* Add independent cell/score oracles, ATT/ATU symmetry and regression tests,
  a theory note and a compact R Markdown demo.

# dml.sensemakr 0.2.2.9011

* Add `method = "dml_refit"` to the did adapters: extract the saved cell data
  and design, then estimate ATT, sigma2.s, nu2.s and all their influence
  functions in one native DML fit. Original ATT/IF values and source control
  trimming are not reused. Source influence functions need not be saved.
* Support built-in dr/reg/ipw MP objects with equal weights and a balanced
  panel. Preserve the existing dml/replicate methods and legacy aliases.
* Integrate batch fitting, covariate benchmarks and pre-trend benchmarks with
  matched folds, repetitions, tuning settings and joint IF inference.
* Reject active DML clipping, invalid overlap or nonregular sensitivity moments;
  add a compact R Markdown demo and native-DML equivalence tests.

# dml.sensemakr 0.2.2.9010

* Explicitly require the control OLS solver to report `converged = TRUE` for
  replicate sensitivity inference. A false, missing or indeterminate flag stops
  the fit even if coefficients are finite. Existing coefficient, prediction,
  matrix and source-replay checks remain in place. Valid fits are unchanged.

# dml.sensemakr 0.2.2.9009

* Support traditional DRDID control residual trimming with `method = "replicate"`,
  including did cells, covariate benchmarks and pre-trend benchmarks. Preserve
  the source ATT/IF while computing sensitivity moments on all source units
  using untrimmed propensities. Recompute each benchmark's mask after refitting.
* Independently reconstruct the source ATT IF, including nuisance effects;
  record control trimming indices, IDs, threshold and retained weight.
* State and derive the correct-control-mean and regular-boundary conditions for
  active trimming. No generic misspecification-robust inference is claimed on
  this path. Source clipping and active trimming in the DML adapter stay rejected.
* Add derivative checks, integration tests, a continuous-boundary simulation
  script and a compact trimming demo.

# dml.sensemakr 0.2.2.9008

* Rename panel-adapter methods to `"dml"` (default) and `"replicate"`.
  Keep `"crossfit"` and `"parametric"` as compatibility aliases, including
  when reading old saved objects for covariate or pre-trend benchmarking.
* Allow `dml_args$dirty.tuning = TRUE` in the DML adapters. Clean tuning remains
  the adapter default. Full-sample tuning selects hyperparameters before
  fold-specific training; benchmarking retains the setting and retunes each fit.
* Record canonical method names and report the tuning mode. The replicate
  method still uses the source logistic-MLE/OLS fits and rejects DML options.
  Statistical assumptions, source ATT/IF preservation and overlap checks remain
  unchanged; full-sample tuning requires suitable selection conditions.

# dml.sensemakr 0.2.2.9007

* Add `pretrend_bounds()` for bias-magnitude and scale-adjusted pre-trend
  benchmarking with ID-aligned joint influence functions and benchmark uncertainty.
* Add `did_pretrend_bounds()` to construct an explicit forward placebo contrast
  on a saved post-treatment cell's units and treatment groups, retaining its
  parametric or repeated cross-fitting method and original post ATT/IF.
* Keep the nonzero-population-placebo condition explicit. Exactly zero estimated
  placebos return point bounds but no delta-method inference at positive k.

# dml.sensemakr 0.2.2.9006

* Add `method = "parametric"` to did single-cell/batch and traditional DRDID
  panel adapters. Keep `method = "crossfit"` as the default.
* Preserve the external ATT/IF and reconstruct full-sample logistic MLE and
  control-arm OLS with the source solvers and saved design. Compute new
  conditional sensitivity moments without cross-fitting or random seeds.
* Include nuisance estimating-equation effects in the moment IFs and joint
  covariance. Document the derivation, sufficient assumptions and limits of
  sensitivity interpretation under misspecification, including benchmarks.
* Retain the parametric method when benchmarking after serialization. Reject
  unequal weights, active source trimming, improved DRDID fits and DML/learner
  overrides on this new path. Existing cross-fitted paths retain their behavior.
* Add source-reproduction, numerical IF differentiation, downstream inference,
  serialized benchmark, batch and rejection tests and a runnable did example.

# dml.sensemakr 0.2.2.9005

* For did objects with unequal positive sampling weights, automatically refit
  the selected cell with equal weights using native DML. Warn that the original
  weighted ATT and IF are not used, and record the mode and original weights.
* Return an ordinary `dml` object on this path so reporting, repeated fitting,
  serialization and benchmarking use the same calculations as direct DML.
  Equal-weight inputs retain the existing hybrid behavior.

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
