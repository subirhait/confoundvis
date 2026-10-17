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

## Common and framework-specific tools

| Tool | Partial R² | ITCV | RIR | E-value | Rosenbaum |
|---|:-:|:-:|:-:|:-:|:-:|
| `plot_robustness_curve()`, `robustness_points()` | ✓ | ✓ | ✓ | ✓ | ✓ |
| `sens_table()` | ✓ | ✓ | ✓ | ✓ | ✓ |
| `sens_report()` | ✓ | ✓ | ✓ | ✓ | ✓ |
| `sens_surface()`, `plot_sensitivity_surface()` | ✓ | ✓ | – | ✓ | ✓ |
| `covariate_impacts()`, `plot_sensitivity_love()` | ✓ | ✓ | – | – | – |
| `plot_sensitivity_contour()` (exact correlation contour) | – | ✓ | – | – | – |

RIR has a single sensitivity parameter, so it has no two-parameter surface.
Covariate benchmarks need a fitted linear model, so they are not available for
E-values or matched pairs.

## Installation

```r
install.packages("confoundvis")
# development version
# install.packages("pak")
pak::pak("subirhait/confoundvis")
```

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
