# Replicate with source control trimming

## Scope

From version 0.2.2.9009, `method = "replicate"` accepts active control trimming
in the traditional DRDID panel estimator, including the default 0.995 rule used
by the supported `did::att_gt(est_method = "dr")` path. No extra option is needed.
The source ATT and ordered IF are preserved after independent reconstruction.
The public `"parametric"` alias has the same behavior.

This extension covers the source's **control residual weights**. It does not
remove observations from the logistic fit, control OLS, or sensitivity moments;
it does not trim treated units. Numerical probability clipping at `1 - 1e-6`,
zero/one probabilities, unequal sampling weights, improved (`imp`) estimation,
and arbitrary custom trimming procedures are outside this extension. The DML
adapter still rejects active source trimming. Existing independent-unit, panel,
data-alignment and source-replay requirements are unchanged.

## Population target and sufficient conditions

Let Y be the outcome change, p = P(D=1),
g(X) = E[Y | D=0,X], and e(X) = P(D=1 | X). The short DiD parameter is
theta_s = E[Y - g(X) | D=1]. Its relation to the causal long parameter and the
OVB bounds still requires the original identifying and sensitivity assumptions.

For a fixed threshold c, define

\[
q_c(D,X)=(1-D)1\{e(X)<c\}\frac{e(X)}{1-e(X)},\qquad A_c=E[q_c]>0.
\]

The probability limit of the source estimator with correct control regression is

\[
T_c=\frac{E[D\{Y-g(X)\}]}p-
    \frac{E[q_c\{Y-g(X)\}]}{A_c}=\theta_s.
\]

The second numerator is zero by conditioning on X and D=0. Thus active source
trimming need not change this target. The equality does **not** hold in general
for an incorrect control regression, even if the propensity model is correct.
This extension does not claim doubly robust inference under active trimming.

One explicit set of sufficient conditions for the implementation is:

1. IID panel units; fixed-dimensional, correctly specified linear control mean
   and logistic propensity models; nonsingular population information matrices.
2. Bounded design covariates and conditional outcome moments of order 4+delta
   for some delta > 0, uniformly over X and D. Propensities are bounded away
   from zero and one, with the upper bound strictly below the source numerical
   clipping cap. These are sufficient conditions, not asserted to be necessary.
3. A fixed, non-data-selected c with A_c > 0. Near c,
   P(|e(X)-c| <= t) <= C*t for small t. A support gap is also sufficient.
   In particular, there must not be an atom at c. A positive observed gap from
   c does not establish a population gap or verify this condition.
4. The original sensitivity/identification assumptions; positive regular
   sensitivity scale for delta-method bounds. Reduced covariate and pre-period
   fits must each satisfy their own model and boundary conditions.

The code checks numerical prerequisites and records trimming diagnostics. It
cannot establish correct population models, population overlap, or the boundary
condition. Active trimming is not a remedy for weak overlap in the full-target
sensitivity moments, which continue to use the original, untrimmed odds.

## Why the preserved source IF is valid under these conditions

Write r = Y-X'beta, with beta the correct control coefficient, and let

\[
H_\beta=E[(1-D)XX'],\quad
\psi_\beta=H_\beta^{-1}(1-D)Xr,\quad
\mu_T=E[DX]/p,\quad \mu_c=E[q_cX]/A_c.
\]

An influence function for T_c in the stated correctly specified parametric
model is

\[
\phi_c=\frac Dp(r-\theta_s)-\frac{q_c}{A_c}r
       +(\mu_c-\mu_T)'\psi_\beta.
\]

Trimming generally makes the OLS correction nonzero; treating fitted nuisance
values as known would omit it. The population derivative with respect to the
propensity parameter is zero: both the smooth odds derivative and the boundary
contribution multiply control residuals whose conditional expectation is zero.

More explicitly, root-n estimation of a fixed-dimensional logistic parameter
changes the mask only near the threshold. Under the boundary and moment
conditions, the L2 norm of that residual-weight difference tends to zero.
The finite-dimensional threshold class is stochastically equicontinuous, so
its centered empirical contribution is o_p(n^-1/2). Its expectation is exactly
zero at the correct g. Multiplying the mask change by the O_p(n^-1/2) error
in beta also gives o_p(n^-1/2). Hence no first-order boundary term remains.
The ordinary smooth beta effect above must still be retained.

The actual DRDID IF uses empirical ratio normalizations and retains a finite-
sample propensity correction as well. Let hats indicate fitted values, let
eta_T and eta_c be the respective weighted residual means, and set

\[
G_\beta=P_n[\widehat q_cX]/\widehat A_c-P_n[DX]/\widehat p,
\quad
G_\gamma=-P_n[\widehat q_c(r-\widehat\eta_c)X]/\widehat A_c.
\]

The adapter independently checks the source vector against

\[
\widehat\phi_c=\frac D{\widehat p}(r-\widehat\eta_T)
 -\frac{\widehat q_c}{\widehat A_c}(r-\widehat\eta_c)
 +G_\beta'\widehat\psi_\beta+G_\gamma'\widehat\psi_\gamma.
\]

Under the stated conditions this is consistent for the population IF above.
If g is misspecified, the missing boundary term need not vanish and T_c need
not equal theta_s. The generic pseudo-true sandwich interpretation in the
untrimmed note must therefore **not** be extended to active trimming. Neither
matching the source IF nor finite-sample numerical differentiation proves this
population validity.

## Sensitivity and benchmarking

The estimators of sigma2, nu2 and their full-sample nuisance corrections are
unchanged from `parametric-inference.md`. They use every source unit and the
raw fitted e, without capping or zeroing it. Conditional on the same nuisance
fits, changing the source trim threshold changes ATT/IF but not these moments.
Their joint covariance uses the aligned vector

\[
(\widehat\phi_c,\widehat\psi_{\sigma^2},\widehat\psi_{\nu^2})'.
\]

The usual delta method then gives bound standard errors. At zero confounding
strength, these reduce to the source's analytic ATT intervals.

Covariate benchmarks keep all units and the source threshold, then refit models
and recompute the mask. They must not reuse the full model's mask. Pre-trend
benchmarks likewise retain all post-cell units and compute their own nuisance
fits and masks on the requested pre-period design. Masks may differ across fits;
joint covariance is aligned by the original unit IDs, not by retained controls.
The existing nonzero-placebo condition for positive-k pre-trend delta inference
continues to apply.

Inspect `m$adapter_info$diagnostics$source.trimming` for the threshold, control
indices/IDs, count, retained weight and observed distance to the cutoff. The
field `boundary.conditions.verified = FALSE` is intentional: numerical
diagnostics do not verify population assumptions. The native `trim.summary`
still describes sensitivity-score clipping, not source control trimming.

## Validation and sources

`test-29-replicate-trimming.R` independently checks weighted numerical
derivatives of ATT, moments and bounds, preservation of the actual source IF,
unchanged full-sample moments, source-mask changes across benchmarks, both did
layouts, batches and pre-trend covariance. These are computational checks.
`inst/validation/replicate-trimming-simulation.R` additionally checks empirical
coverage under correctly specified models with a continuous propensity
distribution through the cutoff; simulation is not a general validity proof.

Source algorithm: DRDID `R/drdid_panel.R`,
https://github.com/pedrohcgs/DRDID/blob/master/R/drdid_panel.R.
General estimated-trimming caveats: Yang and Ding, *Asymptotic causal inference
with observational studies trimmed by the estimated propensity scores*,
https://arxiv.org/abs/1704.00666. The cancellation argument above concerns the
specific control-residual estimator and stated correct-mean model, not all
trimming estimators discussed in that paper.
