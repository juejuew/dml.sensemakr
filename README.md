
<!-- README.md is generated from README.Rmd. Please edit that file -->

# dml.sensemakr: Sensitivity Analysis Tools for Causal ML <img src="man/figures/dml-sensemakr-logo-sm.png" align="right" />

<!-- badges: start -->

<!-- [![CRAN status](https://www.r-pkg.org/badges/version/dml.sensemakr)](https://CRAN.R-project.org/package=dml.sensemakr) -->

<!-- [![Downloads](https://cranlogs.r-pkg.org/badges/dml.sensemakr)](https://cran.r-project.org/package=dml.sensemakr) -->

[![R-CMD-check](https://github.com/carloscinelli/dml.sensemakr/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/carloscinelli/dml.sensemakr/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/carloscinelli/dml.sensemakr/graph/badge.svg)](https://app.codecov.io/gh/carloscinelli/dml.sensemakr)
<!-- badges: end -->

`dml.sensemakr` implements a general suite of sensitivity analysis tools
for Causal Machine Learning as discussed in [Chernozhukov, V., Cinelli,
C., Newey, W., Sharma A., and Syrgkanis, V. (2026). “Long Story Short:
Omitted Variable Bias in Causal Machine Learning.” *Review of Economics
and Statistics*.](https://doi.org/10.1162/REST.a.1705)

# Development version

To install the development version on GitHub make sure you have the
package devtools installed.

``` r
# install.packages("devtools") 
devtools::install_github("carloscinelli/dml.sensemakr")
```

# CRAN

CRAN version coming soon.

# Details

For theoretical details [please see the
paper.](https://doi.org/10.1162/REST.a.1705)

For a primer on Debiased Machine Learning, [please check Chernozukohv et
al (2018).](https://academic.oup.com/ectj/article/21/1/C1/5056401)

Some presentations that may be useful:

- [Carlos’ presentation at the
  ICLR.](http://interactivecausallearning.com/Carlos_Cinelli.html)
- [Victor’s tutorial at the Chamberlain
  Seminar.](https://www.youtube.com/watch?v=PQtYqKfxH_I)

# Basic Usage

``` r
# loads package
library(dml.sensemakr)
#> See details in:
#> - Chernozhukov, V. Cinelli, C. Newey, W. Sharma, A. Syrgkanis, V. (2026). Long Story Short: Omitted Variable Bias in Causal Machine Learning. Review of Economics and Statistics.
#> - Available at: https://doi.org/10.1162/REST.a.1705

## loads data
data("pension")

# set treatment, outcome and covariates
y <- pension$net_tfa  # net total financial assets
d <- pension$e401     # 401K eligibility
x <- model.matrix(~ -1 + age + inc  + educ+ fsize + marr + twoearn + pira + hown, data = pension)

# run DML (nonparametric model)
dml.401k <- dml(y, d, x, model = "npm")
#> Debiased Machine Learning
#> 
#>  Model: Nonparametric 
#>  Target: ate 
#>  Cross-Fitting: 5 folds, 1 reps 
#>  ML Method: outcome (yreg0:ranger, yreg1:ranger), treatment (ranger)
#>  Tuning: dirty 
#> 
#> 
#> ====================================
#> Tuning parameters using all the data
#> ====================================
#> 
#> - Tuning Model for D.
#> -- Best Tune:
#>   mtry min.node.size splitrule
#> 1    2             5  variance
#> 
#> - Tuning Model for Y (non-parametric).
#> -- Best Tune:
#>   mtry min.node.size splitrule
#> 1    2             5  variance
#>   mtry min.node.size splitrule
#> 1    2             5  variance
#> 
#> 
#> ======================================
#> Repeating 5-fold cross-fitting 1 times
#> ======================================
#> 
#> -- Rep 1 -- Folds: 1  2  3  4  5

# summary of results with median method (default)
summary(dml.401k)
#> 
#> Debiased Machine Learning
#> 
#>  Model: Nonparametric 
#>  Cross-Fitting: 5 folds, 1 reps 
#>  ML Method: outcome (yreg0:ranger, yreg1:ranger, R2 = 28.038%), treatment (ranger, R2 = 11.559%)
#>  Tuning: dirty 
#> 
#> Average Treatment Effect: 
#> 
#>         Estimate Std. Error t value  P(>|t|)    
#> ate.all     8010       1163   6.889 5.63e-12 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> Note: DML estimates combined using the median method.
#> 
#> Verbal interpretation of DML procedure:
#> 
#> -- Average treatment effects were estimated using DML with 5-fold cross-fitting. In order to reduce the variance that stems from sample splitting, we repeated the procedure 1 times. Estimates are combined using the median as the final estimate, incorporating variation across experiments into the standard error as described in Chernozhukov et al. (2018). The outcome regression uses  from the R package ; the treatment regression uses Random Forest from the R package ranger.

# sensitivity analysis
sens.401k <- sensemakr(dml.401k, cf.y = 0.04, cf.d = 0.03)

# summary
summary(sens.401k)
#> ==== Original Analysis ====
#> 
#> Debiased Machine Learning
#> 
#>  Model: Nonparametric 
#>  Cross-Fitting: 5 folds, 1 reps 
#>  ML Method: outcome (yreg0:ranger, yreg1:ranger, R2 = 28.038%), treatment (ranger, R2 = 11.559%)
#>  Tuning: dirty 
#> 
#> Average Treatment Effect: 
#> 
#>         Estimate Std. Error t value  P(>|t|)    
#> ate.all     8010       1163   6.889 5.63e-12 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> Note: DML estimates combined using the median method.
#> 
#> Verbal interpretation of DML procedure:
#> 
#> -- Average treatment effects were estimated using DML with 5-fold cross-fitting. In order to reduce the variance that stems from sample splitting, we repeated the procedure 1 times. Estimates are combined using the median as the final estimate, incorporating variation across experiments into the standard error as described in Chernozhukov et al. (2018). The outcome regression uses  from the R package ; the treatment regression uses Random Forest from the R package ranger.
#> 
#> ==== Sensitivity Analysis ====
#> 
#> Null hypothesis: theta = 0 
#> Signif. level: alpha = 0.05 
#> 
#> Robustness Values:
#>         RV (%) RVa (%)
#> ate.all 6.1021  4.6027
#> 
#> Verbal interpretation of robustness values:
#> 
#> -- Robustness Value for the Bound (RV): omitted variables that explain more than RV% of the residual variation of the outcome (cf.y) and generate an additional RV% of variation on the Riesz Representer (cf.d) are sufficiently strong to make the estimated bounds include 0. Conversely, omitted variables that do not explain more than RV% of the residual variation of the outcome nor generate an additional RV% of variation on the Riesz Representer are not sufficiently strong to do so.
#> 
#> -- Robustness Value for the Confidence Bound (RVa): omitted variables that explain more than RV% of the residual variation of the outcome (cf.y) and generate an additional RV% of variation on the Riesz Representer (cf.d) are sufficiently strong to make the confidence bounds include 0, at the  significance level of alpha = 0.05. Conversely, omitted variables that do not explain more than RV% of the residual variation of the outcome nor generate an additional RV% of variation on the Riesz Representer are not sufficiently strong to do so. 
#> 
#>  The interpretation of sensitivity parameters can be further refined for each target quantity. See more below.
#> 
#> Confidence Bounds for Sensitivity Scenario:
#>               lwr       upr
#> ate.all  1572.415 14503.094
#> 
#> Confidence level: point = 95%; region = 90%.
#> Sensitivity parameters: cf.y = 0.04; cf.d = 0.03; rho2 = 1.
#> 
#> Verbal interpretation of confidence bounds:
#> 
#> -- The table shows the lower (lwr) and upper (upr) limits of the confidence bounds on the target quantity, considering omitted variables with postulated sensitivity parameters cf.y, cf.d and rho2. The confidence level "point" is the relevant coverage for most use cases, and stands for the coverage rate for the true target quantity. The confidence level "region" stands for the coverage rate of the true bounds.
#> 
#> Interpretation of sensitivity parameters:
#> 
#> -- cf.y: percentage of the residual variation of the outcome explained by latent variables.
#> -- cf.d: percentage gains in the variation of the Riesz Representer generated by latent variables:
#>    ATE: cf.d measures the percentage gains in the average precision on the treatment regression.

# contout plots
plot(sens.401k)
```

<img src="man/figures/README-basic-usage-1.png" alt="" width="100%" style="display: block; margin: auto;" />

# Acknowledgements

This project was partially supported by the Royalty Research Fund at the
University of Washington, and by the National Science Foundation.

## Panel DRDID adapter

The development adapter `dml_from_drdid()` accepts a compatible DRDID panel fit
with `inffunc = TRUE` and the corresponding estimation data. It preserves the
external ATT and full influence function, estimates conditional ATT sensitivity
moments with cross-fitted DML, and includes their joint covariance in inference.

```r
hybrid <- dml_from_drdid(
  source_fit, data = panel_data,
  yname = "y", tname = "time", idname = "id", dname = "D",
  xformla = ~ x1 + x2,
  dml_args = list(yreg = "lm", dreg = "glm", d.class = TRUE,
                  cf.folds = 3, cf.seed = 17, ps.trim = 0)
)
confidence_bounds(hybrid, cf.y = 0.08, cf.d = 0.06)
robustness_value(hybrid)
dml_benchmark(hybrid, benchmark_covariates = "x2")
```

Initial support is restricted to balanced two-period panels, equal weights,
`imp`/`trad` and no active propensity trimming/clipping. Repeated cross-fitting is available through `dml_args = list(cf.reps = 5)`.
Both estimators are refit for each benchmark. See
[the adapter guide](vignettes/drdid-adapter.Rmd) for assumptions, input alignment
and inference conventions, and the
[complete runnable example](inst/examples/drdid-adapter.R).

## did group-time ATT adapter

`dml_from_did()` adapts one post-treatment `ATT(g,t)` from `did::att_gt()`.
The following describes the default `method = "dml"` path.
It reads the saved estimation data. With equal weights, it preserves the source
ATT and correctly rescaled IF and estimates sensitivity moments with DML.
With unequal positive weights, it warns and refits the whole cell with equal
weights using native `dml()`, including the ATT and IF. The returned ordinary
`dml` object supports native reporting and benchmarking; its
`adapter_info$mode` is `"unweighted_dml"`. This is an unweighted analysis of the
retained sample, not sensitivity analysis for the original weighted ATT.

```r
hybrid <- dml_from_did(source_fit, group = 5, time = 9,
                       dml_args = list(cf.reps = 5))
confidence_bounds(hybrid, cf.y = 0.05, cf.d = 0.05)
benchmark <- dml_benchmark(hybrid, "z")

batch <- dml_from_did_cells(source_fit)
intervals <- did_cell_apply(batch, confidence_bounds, cf.y = 0.05, cf.d = 0.05)
intervals$status  # includes failed cells and their reasons
```

Compatibility is based on saved inputs, DRDID interfaces and numerical ATT/IF
reconstruction, not version equality. Tested with **did 2.5.1 / 2.5.1.902** and
**DRDID 1.3.0**. Both paths require `est_method = "dr"`, balanced panel estimation
and positive finite weights. The hybrid path also requires saved IFs and inactive
trimming; the unweighted path refits all quantities and uses native DML clipping.
Both `faster_mode` settings,
both control groups, anticipation and both base-period settings are supported.
`fix_weights = "varying"`, repeated cross-sections, coarser clustering and
`aggte()` need separate extensions. Inference is analytic and cellwise. Hybrid repetitions retain the source ATT and
IF; mean/median aggregation does not divide its variance by the number of splits.

See `vignette("did-adapter")` and the runnable
[example](inst/examples/did-adapter.R) for installation and batch usage.

### Replicate method (full-sample parametric fitting)

For equal-weight traditional panel estimates, use the source logistic-MLE/OLS
models for sensitivity moments without cross-fitting:

```r
parametric <- dml_from_did(source_fit, group = 5, time = 9,
                           method = "replicate")
confidence_bounds(parametric, cf.y = 0.05, cf.d = 0.05)
dml_benchmark(parametric, "z")
batch <- dml_from_did_cells(source_fit, method = "replicate")
```

This retains the original ATT and unit IF and includes full-sample nuisance
estimation effects in the new moment IFs and joint covariance. It does not
reproduce a source bootstrap SE. Leave `dml_args` empty: there are no folds,
seeds or learner choices. Benchmark refits retain the same parametric method.
Unequal weights, source probability clipping and the improved DRDID estimator
are unsupported on this path. The traditional `dml_from_drdid()` interface
also accepts `method = "replicate"`.

Source **control residual trimming** is supported by replicate from 0.2.2.9009.
It preserves the original ATT/IF and uses all units and untrimmed propensities
for the sensitivity moments. Valid inference requires a correct control mean
and a regular trimming boundary, as well as the existing sensitivity conditions;
it is not generally robust to misspecified outcome models. Diagnostics are in
`m$adapter_info$diagnostics$source.trimming`. See the
[trimming derivation](inst/notes/replicate-trimming.md) and
[short trimming demo](inst/examples/replicate-trimming-demo.Rmd).

The sensitivity interpretation requires suitable model correctness, overlap,
moment and regularity conditions, including for reduced benchmark models.
Reproducing the source numerically does not establish those conditions.
See the [derivation and assumptions](inst/notes/parametric-inference.md) and
the [runnable example](inst/examples/did-parametric.R).


## Panel DiD adapter methods

Use `method = "dml"` for cross-fitted DML sensitivity moments, or
`method = "replicate"` to reproduce the supported source estimator's full-sample
logistic/OLS nuisance fits. On the source-preserving path both retain the
original ATT and influence function. `"crossfit"` and `"parametric"` remain
compatibility aliases. For DML, `dml_args = list(dirty.tuning = TRUE)` enables
full-sample hyperparameter tuning followed by fold-specific fitting; the
adapter default is `FALSE`. See `inst/notes/adapter-methods.md` for details.

## Reestimate the entire did cell with DML

Since 0.2.2.9011, use `method = "dml_refit"` to estimate ATT, `sigma2.s`,
`nu2.s` and their joint influence functions with the same cross-fitted DML fit:

```r
m <- dml_from_did(fit, group = 2007, time = 2007, method = "dml_refit",
  dml_args = list(yreg = "lm", dreg = "glm", d.class = TRUE, cf.seed = 1))
confidence_bounds(m, cf.y = 0.05, cf.d = 0.05)
did_pretrend_bounds(m, pre_period = c(2005, 2006), k = c(0, 1, 2))
```

| Method (equal weights) | ATT and its IF | Sensitivity moments |
|---|---|---|
| `"dml"` (default) | Preserved from did | Cross-fitted DML |
| `"replicate"` | Preserved from did | Full-sample logistic/OLS, corrected IFs |
| `"dml_refit"` | Reestimated by DML | Same cross-fitted DML fit |

The new mode uses the saved balanced panel, cell selection and baseline design;
source ATT/IF values and source control trimming are unused. Source `dr`, `reg`
and `ipw` fits are supported, including `compute_inffunc = FALSE`. Sampling
weights must be equal; higher-level clustering, custom estimators and repeated
cross-sections remain unsupported. Active DML propensity clipping stops the fit.
Covariate benchmarks and pre-trend benchmarks also use DML with matched folds.
See `inst/examples/dml-refit-demo.Rmd` for a short runnable example. Inference
still needs identification, overlap, moment and nuisance-rate conditions;
optional dirty tuning requires suitable selection conditions.

## ATU for the selected DiD controls

Since 0.2.2.9012, `dml_refit` also supports ATU:

```r
m_atu <- dml_from_did(fit, group = 2007, time = 2007,
  method = "dml_refit", target = "atu",
  dml_args = list(yreg = "lm", dreg = "glm", d.class = TRUE, cf.seed = 1))
summary(m_atu)
dml_bounds(m_atu, cf.y = 0.05, cf.d = 0.05)
```

The target is the effect on the **selected controls** if they adopted in
`group`, evaluated at `time`. It requires conditional mean comparability of
potential changes under that adoption path; ordinary untreated parallel trends
alone does not identify ATU. Sensitivity analysis can impose this condition
given observed and omitted covariates, with explicit confounding-strength bounds.
These are identifying assumptions, not properties verified by fitting DML.

ATU uses treated-arm outcome fitting and inverse odds. All three estimates,
IFs and sensitivity moments come from DML. Bounds, robustness values,
covariate benchmarks and batch conversion support ATU. `replicate`, hybrid
`dml`, and pre-trend benchmarking reject ATU. The existing ATT behavior remains
the default. See `inst/notes/did-atu.md` and `inst/examples/atu-demo.Rmd`.
