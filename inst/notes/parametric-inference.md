# Full-sample parametric panel sensitivity mode

## Scope and interface

Use `dml_from_did(fit, group, time, method = "parametric")`,
`dml_from_did_cells(fit, method = "parametric")`, or
`dml_from_drdid(fit, ..., method = "parametric")` for a traditional panel fit.
The default remains `method = "crossfit"`. The new method requires an empty
`dml_args`, equal positive sampling weights, no active source trimming or
clipping, and the existing adapter's balanced-panel/source-data/IF checks.
It supports the traditional logistic-MLE/OLS estimator only. It does not call
`dml()` with one fold, select learners, or use random fold assignments.

The source ATT and ordered influence function are retained. For did, the
global IF is first rescaled to the selected cell, exactly as in the existing
adapter. The analytic ATT SE is recovered from this IF; a bootstrap SE or
simultaneous critical value from the source is not reproduced. Nuisance models
are reconstructed using `DRDID`'s `fastglm_fit` solvers and the saved design,
including intercepts and transformations. Replay checks use numerical tolerance
1e-7. They do not certify model correctness or bitwise identity across solvers
and package versions.

The sensitivity moments below are new estimates, not values stored in the
original did/DRDID result. `fits[[1]]$preds` stores their full-sample predictions;
`parametric` stores coefficients, gradients and assumptions. There is one result,
`info$cf.folds` is zero, and fold/seed metadata are empty.

## Estimators and sufficient assumptions

All notation refers to the frozen two-period cell. Let Y be the outcome change,
D its treatment-group indicator, and X its saved design (including any intercept).
With n independent units, write Pn for the empirical mean. Fit

\[
g_\beta(X)=X^\top\beta,\quad
e_\gamma(X)=\operatorname{expit}(X^\top\gamma),\quad
\widehat p=P_nD.
\]

OLS uses only D=0; logistic MLE uses the entire cell. Define

\[
r_i=Y_i-X_i^\top\widehat\beta,\quad
l_i=D_i/\widehat p,\quad c_i=(1-D_i)/(1-\widehat p),\quad
h_i=\frac{\widehat e_i}{1-\widehat e_i}\frac{1-\widehat p}{\widehat p}.
\]

The estimates use the same uncentered conditional moment parameterization as
`att.npm.cond`:

\[
\widehat a=\widehat\sigma^2=\frac{P_n[c r^2]}{P_n c},\qquad
\widehat b=\widehat\nu^2=\frac{P_n[2lh-ch^2]}{P_n[2l-c]},\qquad
\widehat S^2=\widehat a\widehat b.
\]

Sufficient conditions for their usual sensitivity interpretation include the
original identifying assumptions; correct conditional mean and propensity
models; fixed model dimension; regular, nonsingular nuisance estimation;
overlap and sufficient moments. Under such conditions, ordinary smooth
estimating-equation theory applies without cross-fitting. The code checks
finite fits, nonsingular matrices and positive sensitivity scale, but it cannot
verify these population assumptions. No arbitrary dimension-to-sample cutoff
is represented as a statistical guarantee.

Under misspecification, regular sandwich IFs can describe the corresponding
pseudo-true parametric functionals. They do not turn those functionals into the
nonparametric sensitivity targets. ATT double robustness alone does not suffice
for sensitivity validity. A reduced benchmark model has its own conditions:
dropping a covariate can make even a previously correct logistic/linear model
misspecified. Benchmark output is a model-refit comparison, with a causal
calibration interpretation only when the relevant conditions also hold for
the reduced fit.

## Nuisance and moment influence functions

Use the empirical Jacobians

\[
H_\beta=P_n[(1-D)XX^\top],\qquad
H_\gamma=P_n[e(1-e)XX^\top].
\]

The unit IFs of the nuisance estimates are

\[
\psi_{\beta,i}=H_\beta^{-1}(1-D_i)X_i r_i,\quad
\psi_{\gamma,i}=H_\gamma^{-1}X_i(D_i-e_i),\quad
\psi_{p,i}=D_i-\widehat p.
\]

Start with the fixed-nuisance ratio IFs

\[
\psi^0_{a,i}=\frac{c_i(r_i^2-\widehat a)}{P_n c},\qquad
\psi^0_{b,i}=\frac{2l_i h_i-c_i h_i^2-(2l_i-c_i)\widehat b}{P_n[2l-c]}.
\]

Differentiate each empirical ratio holding the empirical distribution fixed:

\[
A_\beta=-\frac{2P_n[(1-D)rX]}{P_n(1-D)},\qquad
B_\gamma=\frac{2P_n[(lh-ch^2)X]}{P_n[2l-c]}.
\]

For the p derivative, write Q=2lh-ch^2, V=2l-c and

\[
l_p=-l/p,\quad c_p=c/(1-p),\quad h_p=-h/[p(1-p)],
\]
\[
B_p=\frac{P_n[2l_p h+2l h_p-c_p h^2-2chh_p-\widehat b(2l_p-c_p)]}{P_n V}.
\]

The full moment IFs are

\[
\psi_{a,i}=\psi^0_{a,i}+A_\beta^\top\psi_{\beta,i},\qquad
\psi_{b,i}=\psi^0_{b,i}+B_\gamma^\top\psi_{\gamma,i}+B_p\psi_{p,i}.
\]

The OLS normal equations make A_beta zero up to solver error. The particular
ratio normalization makes B_p zero algebraically at p=PnD, even in the empirical
sample. Both derivatives are kept explicit in code for auditability. The
propensity correction B_gamma need not vanish in finite samples and should not
be silently discarded. Under correct models, its population value vanishes by
orthogonality. Retaining all derivatives also supports regular pseudo-true
functional inference without claiming misspecification-robust causal bounds.

## Joint covariance, bounds and benchmarks

Let phi_i be the source ATT IF on the cell scale. Assemble

\[
\Psi_i=(\phi_i,\psi_{a,i},\psi_{b,i})^\top,\qquad
\widehat V=\frac1{n^2}\sum_i(\Psi_i-\bar\Psi)(\Psi_i-\bar\Psi)^\top.
\]

The source ATT IF already contains its own normalization and nuisance effects;
no second ATT nuisance correction is added. The sensitivity scale IF is
`b * psi.a + a * psi.b`. Existing bound and benchmark routines consume these
jointly aligned IFs. At fixed confounding strength k, the lower bound is
`theta - k * sqrt(a*b)` with gradient
`(1, -k/2*sqrt(b/a), -k/2*sqrt(a/b))`. Zero-strength intervals reduce to the
source cell's analytic normal intervals. Statements remain cellwise.

Benchmark refits freeze units, outcomes and the surviving design columns, then
refit both DRDID and these full-sample moments. They retain the same method
after serialization and do not evaluate the original call. Learner overrides
are rejected instead of silently changing the method.

## Independent validation

`tests/testthat/test-26-adapter-parametric.R` checks source preservation,
predictions against base R OLS/logistic fits, downstream inference and benchmark
and batch behavior. Crucially, it independently perturbs each observation's
empirical mass, refits using `glm.fit`/`lm.wfit`, and numerically differentiates
ATT, sigma2, nu2, S2 and the lower sensitivity bound. The derivatives are compared
with all corresponding analytic IFs. These perturbation weights are solely a
derivative oracle; they do not extend the public adapter to unequal weights.

Numerical differentiation validates the implemented estimating equations and
chain rule. It is not a coverage theorem or a test of model correctness.

References: Sant'Anna and Zhao (2020), *Doubly Robust Difference-in-Differences
Estimators*, https://doi.org/10.1016/j.jeconom.2020.06.003; Chernozhukov et al.
(2018), *Double/debiased machine learning for treatment and structural
parameters*, https://doi.org/10.1111/ectj.12097. The source replay follows
https://github.com/pedrohcgs/DRDID/blob/master/R/drdid_panel.R.
