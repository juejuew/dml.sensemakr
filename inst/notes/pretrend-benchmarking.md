# Pre-trend benchmarking for panel adapters

This implementation follows Section 4.2, Section 5.2 and Appendix C.2 of
Wang, Sant'Anna, Chernozhukov and Cinelli, *Omitted Variable Bias in
Difference-in-Differences Designs* (2026). In the supplied manuscript these
appear on pages 18-20, 24 and 53. This is a benchmark based on one explicitly
chosen pre-treatment comparison, not a general multiple-pre-period procedure.

## Two restrictions

Let theta be the post-treatment short ATT and q the pre-treatment placebo ATT.
Let S = sqrt(a*b), where a = sigma2 and b = nu2 are the uncentered conditional
sensitivity moments. Define Spre analogously. For each user-specified k >= 0:

* `type = "bias"`: lower/upper = theta +/- k * abs(q).
* `type = "scaled"`: lower/upper = theta +/- k * (S/Spre) * abs(q).

These encode different transport assumptions. k is a chosen multiplier and is
not estimated or validated by a pre-trend test. The second restriction permits
the observable sensitivity scale to change between periods. It requires
strictly positive finite scale components in both fits. Period lengths are not
automatically rescaled; k refers to outcome changes over the selected intervals.

The true causal effect in the pre comparison must be zero, including no
anticipation. Validity also requires the identifying/transport assumptions,
overlap, moments and nuisance conditions appropriate to each fit. For Option 2,
use the regular fixed-dimensional logistic/OLS conditions described in
`parametric-inference.md`. The two fits' predictions are allowed to differ.

## Joint inference

Every unit contributes the six-vector

    (psi.theta, psi.a, psi.b, psi.q, psi.apre, psi.bpre).

Vectors are aligned by explicit unit ID; the covariance of the six estimates is
the centered cross-product divided by n^2. For R = S/Spre,

    psi.R = R/2 * (psi.a/a + psi.b/b - psi.apre/apre - psi.bpre/bpre).

For the bias method set R = 1 and psi.R = 0. The IF of the bias width is

    psi.B = k * (R * sign(q) * psi.q + abs(q) * psi.R).

The bound IFs are psi.theta +/- psi.B. This is algebraically Appendix C.2's
formula, written without division by q. The formulas use the original ATT IFs
and the adapter's already-corrected moment IFs; they do not add another nuisance
correction or assume independent pre/post estimates.

For repeated cross-fitting, matching partitions and repetition order are
required. Compute each repetition's bound and its IF SE, then apply the existing
package mean/median combination rule (including between-repetition spread,
without dividing the SE by the square root of the repetition count). No failed
repetition is dropped. Results are pointwise in k and in the selected cell.

`level` is the one-sided confidence level of each endpoint, consistent with
`benchmark_bounds()`. The default is 0.95; `level = 0.975` yields at least 95%
asymptotic coverage of both population endpoints by the union bound, under the
stated regularity conditions. Original source bootstrap critical values and
simultaneous bands are not reused.

## The zero-pre-trend boundary

Appendix C.2 assumes the population q is nonzero. The absolute value is
nondifferentiable at zero. Writing `sign(0) = 0` would incorrectly erase the
benchmark's uncertainty; this implementation does not do so. If the estimated
q is exactly zero, positive-k rows retain their plug-in point bounds but have
NA SEs/confidence bounds, a `zero_pretrend` status, and a warning. The k = 0 row
still recovers the post ATT/analytic IF interval. A nonzero estimate does not
establish the population condition; delta inference is not uniformly justified
near zero. No threshold, p-value screening or nonregular alternative is claimed.

## did convenience interface

    m <- dml_from_did(fit, group, time, method = "parametric")
    b <- did_pretrend_bounds(m, pre_period = c(base, end), k = c(0, 1, 2))
    b$bounds
    pretrend_bounds(m, b$pre_model, k = c(0, 1, 2), type = "scaled")

The post ATT and IF remain unchanged. The original MP object and its saved panel
must still be present. The pre fit is a new DRDID placebo on precisely the post
cell's units and treatment-group labels, using the saved design at `base`.
Both dates must precede the cohort's anticipation-adjusted treatment threshold;
`end` must be no later than the post cell's baseline.

This construction handles both fast and slow did layouts and both source base
period conventions. It deliberately does not copy a stored pre cell: under
not-yet-treated controls that cell can use different units, and under a universal
baseline it can be normalized to zero or use a reversed time contrast. The new
object records `sample_policy = "post_cell"` and its requested pre periods.

Cross-fitting settings are inherited from m, giving matched partitions. The
function can be used as FUN in `did_cell_apply()` for batch cells. Saving and
restoring the objects preserves their method. The original post-only adapter
continues rejecting pre-treatment cells; there is no silent change to its API.

## Explicit pairs and current limits

`pretrend_bounds(post_model, pre_model, ...)` also works with explicitly supplied
DRDID adapters. Both fits must be equal-weight conditional ATT adapters with
the same units, treatment labels and method. Supply `unit_id` for wide DRDID
data. Different row order is allowed and aligned; different samples are rejected.
The caller must verify pre-treatment timing for explicit wide-data fits.

This version does not add unequal-weight inference, higher-level clustering,
different-sample transport, repeated cross-sections, aggregation, multiple-pretrend
selection or maxima, automatic k selection, or simultaneous confidence bands.

## Validation

The new test file compares both bounds' analytic IFs against independent
empirical-mass perturbations: each unit is reweighted, both models are refit
with base-R glm.fit/lm.wfit, and both scalar bound endpoints are recomputed.
It also checks cross-period covariance, sign reversal, k = 0, unit reordering,
serialization, zero/nonpositive scales, invalid dates and IDs, repeated folds,
both did layouts/base conventions, not-yet-treated sample freezing, anticipation
and batch behavior. These are implementation checks, not a coverage simulation.
