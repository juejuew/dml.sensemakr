
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
It reads the saved estimation data, preserves the source ATT and correctly
rescaled IF, and estimates sensitivity moments with cross-fitted DML.

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
**DRDID 1.3.0**. Required scope: `est_method = "dr"`, saved IFs, balanced panel estimation,
constant positive weights and inactive trimming. Both `faster_mode` settings,
both control groups, anticipation and both base-period settings are supported.
`fix_weights = "varying"`, repeated cross-sections, coarser clustering and
`aggte()` need separate extensions. Inference is analytic and cellwise. All repetitions retain the source ATT and
IF; mean/median aggregation does not divide its variance by the number of splits.

See `vignette("did-adapter")` and the runnable
[example](inst/examples/did-adapter.R) for installation and batch usage.
