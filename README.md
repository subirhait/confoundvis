
<!-- README.md is generated from README.Rmd. Please edit that file -->

# confoundvis

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/confoundvis)](https://CRAN.R-project.org/package=confoundvis)
[![R-CMD-check](https://github.com/subirhait/confoundvis/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/subirhait/confoundvis/actions/workflows/R-CMD-check.yaml)
[![License:
GPL-3](https://img.shields.io/badge/License-GPL--3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
<!-- badges: end -->

`confoundvis` draws and reports sensitivity analyses for unmeasured
confounding. A single `confoundsens` object stores a *sensitivity path*,
the treatment estimate as a function of the strength of a hypothetical
omitted confounder, whichever framework produced it. The same robustness
curves, covariate benchmark plots, contour plots, and plain-language
reports then work for:

| Framework                            | Reference                         | Strength index                      | Path source                          |
|--------------------------------------|-----------------------------------|-------------------------------------|--------------------------------------|
| Impact threshold (ITCV)              | Frank (2000); Frank et al. (2013) | impact `r(D,U) x r(Y,U)`            | `itcv_lm()`, `from_konfound()`       |
| Partial R-squared / robustness value | Cinelli & Hazlett (2020)          | partial R-squared of the confounder | `sens_path_lm()`, `from_sensemakr()` |
| E-value                              | VanderWeele & Ding (2017)         | confounder risk ratio               | `from_evalue()`                      |

## Installation

``` r
install.packages("confoundvis")
# development version
# pak::pak("subirhait/confoundvis")
```

## Workflow

``` r
library(confoundvis)
fit <- lm(mpg ~ am + wt + hp + qsec, data = mtcars)

# 1. sensitivity paths computed from the fitted model
path <- sens_path_lm(fit, treatment = "am")    # partial R-squared
it   <- itcv_lm(fit, treatment = "am")         # ITCV

# 2. plots
plot_robustness_curve(path)
plot_robustness_curve(it$path)

# 3. benchmark against observed covariates
imp <- covariate_impacts(fit, "am")
plot_sensitivity_love(imp)
plot_sensitivity_contour(attr(imp, "threshold"), benchmarks = imp)

# 4. report
sens_report(path)
```

Results already produced by **sensemakr**, **konfound**, or **EValue**
can be converted with `from_sensemakr()`, `from_konfound()`, and
`from_evalue()`; `as_confoundsens()` also accepts a data frame of
precomputed paths, including multilevel (within/between) paths.

See `vignette("confoundvis-workflow")` for a complete example with the
public `darfur` data.

## Scope

`confoundvis` is a presentation layer. Its computations reproduce each
framework’s published formulas (tests compare them with sensemakr,
konfound, and EValue), and it inherits each framework’s assumptions. A
sensitivity display shows how strong confounding would have to be; it
cannot show whether such a confounder exists, and it cannot repair a
flawed identification strategy. `plot_reversal_cone()` and
`plot_taylor_panels()` are conceptual illustrations built on stylized
models.

## References

Cinelli, C., & Hazlett, C. (2020). Making sense of sensitivity:
Extending omitted variable bias. *JRSS-B*, 82(1), 39–67.

Frank, K. A. (2000). Impact of a confounding variable on a regression
coefficient. *Sociological Methods & Research*, 29(2), 147–194.

Frank, K. A., Maroulis, S. J., Duong, M. Q., & Kelcey, B. M. (2013).
What would it take to change an inference? *Educational Evaluation and
Policy Analysis*, 35(4), 437–460.

VanderWeele, T. J., & Ding, P. (2017). Sensitivity analysis in
observational research: Introducing the E-value. *Annals of Internal
Medicine*, 167(4), 268–274.

## Citation

``` r
citation("confoundvis")
```
