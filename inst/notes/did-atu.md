# ATU in a saved DiD cell

Version 0.2.2.9012 adds `target = "atu"` to `dml_from_did()` and
`dml_from_did_cells()`, only with `method = "dml_refit"`. ATT remains the
default. The estimator reuses the native conditional ATU DML engine; it does
not relabel, negate, preserve, or reuse the source ATT or source ATT IF.

## Target population and identification

Fix a selected cell `(g,t)` and its saved unaffected baseline `b`. Within that
cell, `D=1` denotes adoption cohort `g`, and `D=0` denotes the controls chosen
by did's saved control-group and anticipation rules. In particular, with
not-yet-treated controls the target population can change with `t`.

Let `Y_t(g)` be the potential outcome under adoption in `g`, and `Y_t(0)` the
potential outcome under no adoption through `t`. The target is

\[
\theta_U(g,t)=E[Y_t(g)-Y_t(0)\mid D=0].
\]

It asks what would happen to these selected controls if they adopted in `g`.
It is not an ATT of cohort `g`, nor an effect for every untreated unit in the
original dataset. Later adopters are controls only when their factual outcomes
through `t` are unaffected under the saved anticipation rule. Baseline outcomes
and baseline covariates must be unaffected by both the hypothetical and factual
treatment paths. Consistency, no interference and a well-defined common
treatment path are required. Higher-level clustering remains unsupported.

Define `A=Y_t(g)-Y_b(0)` and `Z=Y_t-Y_b`. For observed covariates `X`, a
sufficient identification restriction is

\[
E[A\mid X,D=1]=E[A\mid X,D=0]. \tag{PT1}
\]

This concerns potential changes under adoption in `g`, including any permitted
anticipation after `b`. It is not ordinary untreated parallel trends (PT0).
Writing `m1(X)=E[Z|X,D=1]`, consistency and PT1 give

\[
\begin{aligned}
E[m_1(X)-Z\mid D=0]
 &=E[Y_t(g)-Y_b(0)-(Y_t(0)-Y_b(0))\mid D=0]\\
 &=\theta_U(g,t).
\end{aligned}
\]

PT0 alone cannot provide this conclusion. For example, all baseline and
untreated post outcomes can be zero, observed treated post outcomes can be 1,
and the never-treated group's unobserved treated post outcomes can be any
number `a`. Observed data and PT0 are unchanged as `a` varies, but ATU equals
`a`. No reweighting algorithm or DML learner resolves this identification gap.

For sensitivity analysis, PT1 may hold given `(X,U)` instead of `X` alone.
The reported short target is `theta_U,s=E[m1,s(X)-Z|D=0]`. The long target uses
`m1(X,U)` and is causal under the corresponding long-model PT1. The sensitivity
parameters restrict how much omitted `U` changes the treated-arm outcome
regression and the treatment-assignment odds. A choice of bounds on those
parameters is substantive; numerical checks cannot verify PT1 or those bounds.
`adapter_info$identification$verified` is therefore always `FALSE`.

## Score and influence function

Let `e(X)=P(D=1|X)` and `p=P(D=1)`. For the short target, a doubly robust score
is

\[
a_U(Z,D,X)=\frac{1-D}{1-p}(m_1(X)-Z)
 +\frac{D(1-e(X))}{e(X)(1-p)}(Z-m_1(X)).
\]

At the true nuisance functions its expectation equals the short target. The
population effect IF is

\[
\psi_{\theta_U}=a_U-\frac{1-D}{1-p}\theta_{U,s}.
\]

The sign is treatment minus no treatment. Both outcome prediction for controls
and residual correction for treated units enter with this sign convention.
For cross-fitting, all nuisance values, including `p`, are obtained on training
folds. The actual finite-sample estimator uses
`mean(a_U)/mean((1-D)/(1-phat))`; its IF divides the corresponding centered
score by the same denominator. Thus fold-specific normalization is preserved,
not replaced by a full-sample approximation.

## Sensitivity moments and symmetry

The conditional ATU scale uses

\[
\sigma_{1,s}^2=E[(Z-m_1(X))^2\mid D=1],\qquad
h_1(X)=\frac{p}{1-p}\frac{1-e(X)}{e(X)},\qquad
\nu_{1,s}^2=E[h_1(X)^2\mid D=1].
\]

Here `h1` is the ratio of control to treated covariate densities. The uncentered
sensitivity scale is `S1^2=sigma1^2*nu1^2`. Its orthogonal moment and IF use
`l=(1-D)/(1-p)`, `lc=D/p`:

\[
\psi_{\sigma_1^2}=l_c\{(Z-m_1)^2-\sigma_{1,s}^2\},\qquad
\psi_{\nu_1^2}=2lh_1-l_ch_1^2-(2l-l_c)\nu_{1,s}^2.
\]

The native finite-sample ratios additionally divide the variance score by
`mean(lc)` and the imbalance score by `mean(2*l-lc)`. These normalizations, the
product-rule IF for `S1^2`, and all covariance signs are retained unchanged.

The validity argument is an exact role transformation of conditional ATT:
replace `D` by `1-D`, `p` by `1-p`, `e` by `1-e`, and `m0` by `m1`.
Under this transformation, the ATT effect and its IF are the negatives of the
original ATU effect and IF. The variance, imbalance and squared scale, and
their IFs, are unchanged. In particular, effect-scale covariance changes sign.
The absolute OVB bound is unchanged by this sign transformation. The existing
ATU covariate-benchmarking code uses the corresponding ATU alignment sign.

This algebra transfers the sensitivity calculation, not the factual PT0
assumption: the relabeled problem needs PT1 for causal identification. For the
long-model bound it needs PT1 given `(X,U)` and appropriate overlap there.

## Inference and supported operations

Inference requires independent sampled panel units, nonvanishing group shares,
overlap with adequate treated support for controls, finite score variances,
positive finite sensitivity scales, and consistent nuisance estimates with the
remainder-rate conditions for all three orthogonal scores. For example,
suitable weighted squared errors for the outcome and imbalance moments must
be `o_p(n^-1/2)`; effect-score nuisance products must also be negligible.
Double robustness of the effect alone does not validate the sensitivity
moments. Optional dirty tuning needs additional selection/overfitting control.
Active DML clipping, invalid propensities, deficient outcome-training folds and
nonpositive sensitivity scales stop the cell. Sampling weights must be equal.

Results are in `results$main$untr`, including a three-parameter joint IF
covariance matrix. The ATU and its IF are reestimated even if the source ATT
is missing or its IF was not saved. No new full-sample parametric correction
is added to the native DML IFs. The source `original.se`, if stored, still
refers to ATT, as indicated by `original.target = "att"`.

Fixed-strength bounds, robustness values, covariate benchmarks, benchmark
bounds, summary/printing, serialization and batch conversion use ATU. Reduced
covariate fits retain the target, units, folds, repetitions and tuning settings.
Each reduced fit needs its own nuisance regularity conditions. Batch estimates
remain cellwise; target controls can overlap across cells, so no aggregation
into an overall ATU or simultaneous inference is supplied.

`replicate` and hybrid `dml` reject ATU because their source estimates are ATT.
Both pre-trend functions reject ATU before fitting. Untreated pre-period
contrasts do not by themselves benchmark selection in treated potential
changes or treatment effects. A separate explicit ATU transport restriction
would be necessary; the existing ATT pre-trend assumption is not silently
reused. Existing ATT pre-trend functionality is unchanged.

## Use

```r
m <- dml_from_did(fit, group = 2007, time = 2007,
  method = "dml_refit", target = "atu",
  dml_args = list(yreg = "lm", dreg = "glm", d.class = TRUE, cf.seed = 1))
summary(m)
dml_bounds(m, cf.y = 0.05, cf.d = 0.05)
```

The call estimates the short ATU functional. A causal interpretation and the
sensitivity bounds rely on the assumptions above; the API does not establish
them from the observed data.
