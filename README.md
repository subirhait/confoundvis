# confoundvis

<!-- badges: start -->
[![R-CMD-check](https://github.com/subirhait/confoundvis/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/subirhait/confoundvis/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

**confoundvis** is a common interface to established sensitivity analyses for
unmeasured confounding. Each framework keeps its own parameter and its own
published formula. confoundvis stores the results in one object. Robustness
curves, tables, and reports work for every framework; surfaces, contours, and
covariate benchmarks are provided only where the framework defines them. It
adds no new sensitivity parameter or measure.

| Framework | Source | confoundvis function | Agrees with |
|---|---|---|---|
| Partial R² / robustness value | Cinelli & Hazlett (2020) | `sens_path_lm()`, `from_sensemakr()` | `sensemakr` |
| Impact threshold (ITCV) | Frank (2000) | `itcv_lm()`, `from_konfound()` | `konfound` |
| Robustness of inference to replacement (RIR) | Frank et al. (2013) | `rir_lm()`, `from_konfound(index = "RIR")` | `konfound` |
| E-value | VanderWeele & Ding (2017) | `from_evalue()` | `EValue` |
| Rosenbaum bounds, matched pairs | Rosenbaum (2002) | `sens_path_pairs()` | `DOS2`, `sensitivity2x2xk`, `sensitivitymv` |
| Risk difference, binary confounder (logistic model) | Lin, Psaty & Kronmal (1998) | `sens_risk_difference()` | `marginaleffects` (no confounding), `episensr` (RD correction) |
| Instrumental variables (Anderson–Rubin) | Cinelli & Hazlett (2025) | `from_iv_sensemakr()` | `iv.sensemakr` |
| Double machine learning (PLR, IRM-ATE) | Chernozhukov et al. (2024) | `sens_dml()` for R `DoubleML` fits | Python `doubleml` |

## Common and framework-specific tools

| Tool | Partial R² | ITCV | RIR | E-value | Rosenbaum | Risk diff. | IV | DML |
|---|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
| `plot_robustness_curve()`, `robustness_points()` | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| `sens_table()` | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| `sens_report()` | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| `sens_surface()`, `plot_sensitivity_surface()` | ✓ | ✓ | – | ✓ | ✓ | ✓ | ✓ | ✓ |
| `covariate_impacts()`, `plot_sensitivity_love()` | ✓ | ✓ | – | – | – | – | – | – |
| `plot_sensitivity_contour()` (exact correlation contour) | – | ✓ | – | – | – | – | – | – |

RIR has a single sensitivity parameter, so it has no two-parameter surface.
Covariate benchmarks need a fitted linear model, so they are not available for
E-values or matched pairs. For instrumental variables, the benchmark bounds
computed by `iv.sensemakr` are shown on the surface. The R package `DoubleML`
has no sensitivity analysis; `sens_dml()` computes the one of the Python
package `doubleml` from a fit with stored predictions and agrees with it.

## Installation

```r
install.packages("confoundvis")
# development version
# install.packages("pak")
pak::pak("subirhait/confoundvis")
```

## Quick start: everything in one call

```r
library(confoundvis)
fit <- lm(mpg ~ am + wt + hp, data = mtcars)
res <- confoundvis(fit, "wt")   # partial R2, ITCV, RIR, and covariate benchmarks
res                             # table and plain-language summary
plot(res)                       # one panel per analysis
plot(res, which = "benchmarks") # measured covariates
```

`confoundvis()` also accepts a logistic regression (RIR through
'konfound'), a vector of matched-pair differences (Rosenbaum bounds), an
instrumental-variable analysis from 'iv.sensemakr', or a `DoubleMLPLR` /
`DoubleMLIRM` fit (with `store_predictions = TRUE`).

```r
# binary outcome: risk difference with a practical threshold
# (right heart catheterization data, Connors et al., 1996). ATbounds documents
# `survival` as 30-day survival, but it matches the death indicator over the
# whole follow-up (about 180 days) of the original Vanderbilt release.
data("RHC", package = "ATbounds")
rhc <- RHC
rhc$death <- 1 - rhc$survival     # death during follow-up (about 180 days)
rhc$survival <- NULL
fit <- glm(death ~ ., family = binomial, data = rhc)
rd <- sens_risk_difference(fit, "RHC", p1 = 0.4, p0 = 0.2, threshold = 0.02)
sens_report(rd)
plot_robustness_curve(rd)

# instrumental variables
data("card", package = "iv.sensemakr")
xm <- model.matrix(~ exper + expersq + black + south + smsa + reg661 + reg662 +
                     reg663 + reg664 + reg665 + reg666 + reg667 + reg668 + smsa66,
                   data = card)
iv <- iv.sensemakr::iv_fit(card$lwage, card$educ, card$nearc4, xm)
confoundvis(iv.sensemakr::sensemakr(iv, benchmark_covariates = c("black", "smsa")))

# double machine learning
library(DoubleML)
dat <- make_plr_CCDDHNR2018(n_obs = 500, return_type = "data.table")
plr <- DoubleMLPLR$new(DoubleMLData$new(dat, y_col = "y", d_cols = "d"),
                       ml_l = mlr3::lrn("regr.lm"), ml_m = mlr3::lrn("regr.lm"))
plr$fit(store_predictions = TRUE)
confoundvis(plr)
```
The functions below give finer control.

## Example

```r
library(confoundvis)
data("darfur", package = "sensemakr")
fit <- lm(peacefactor ~ directlyharmed + age + farmer_dar + herder_dar +
            pastvoted + hhsize_darfur + female + village, data = darfur)

p_r2 <- sens_path_lm(fit, treatment = "directlyharmed")   # sensemakr's RV
it   <- itcv_lm(fit, treatment = "directlyharmed")        # konfound's ITCV
rr   <- rir_lm(fit, treatment = "directlyharmed")         # konfound's RIR

sens_table(partial_r2 = p_r2, itcv = it$path, rir = rr$path)
plot_sensitivity_surface(sens_surface(p_r2))
sens_report(p_r2)
```

See `vignette("confoundvis-workflow")` for the full workflow, including
E-values and matched pairs.

## What confoundvis does not do

A sensitivity display shows how strong confounding would have to be to change
a conclusion. It cannot show whether such a confounder exists, and it inherits
the assumptions of the framework that produced it.
