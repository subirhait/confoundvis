# confoundvis 0.2.0

## Bug fixes

* Added the missing `tests/testthat.R` runner. In 0.1.0 the test files were
  shipped but never executed by `R CMD check`.
* `grid` is now declared in `Imports` (it was used via `grid::`).
* `validate_confoundsens()` now checks that `lambda` is nondecreasing
  *within* each level, so stacked multilevel paths no longer trigger a
  spurious warning.

## New features: from a fitted model to a report

* `sens_path_lm()` computes the partial R-squared sensitivity path
  (Cinelli and Hazlett, 2020) directly from an `lm` fit, including
  weighted fits.
* `itcv_lm()` computes the conditional and unconditional impact threshold
  (Frank, 2000; Frank et al., 2013) and the ITCV path.
* `covariate_impacts()` benchmarks every observed covariate term (including
  multi-column factor terms for the partial R-squared metric) and returns a
  data frame that `plot_sensitivity_love()` and `plot_sensitivity_contour()`
  accept directly.
* `robustness_points()` locates where a path reaches the null and where
  significance is lost; `sens_report()` turns these into a short
  plain-language report that states what the analysis cannot establish.

## Interoperability

* `as_confoundsens()` is now generic, with methods for data frames (as
  before), `lm` fits, and `sensemakr` objects.
* `from_sensemakr()`, `from_konfound()`, and `from_evalue()` convert results
  from the 'sensemakr', 'konfound', and 'EValue' packages. These packages are
  only suggested; the adapters read documented fields and do not call them.
* The test suite checks the computations against those packages when they are
  installed (adjusted estimates and standard errors, robustness values,
  conditional and unconditional ITCV, and E-values).

## Other changes

* `confoundsens` objects gain an optional `meta` list (framework, labels,
  null value, critical value). Objects created by 0.1.0 remain valid.
* `plot_robustness_curve()` uses axis labels from `meta` and gains
  `reference` to draw the null and critical values.
* `plot_sensitivity_love()` takes its threshold from the `"threshold"`
  attribute of `df` when not supplied, and drops rows with missing impact.
* `plot_figure2_taylor_panels()` is renamed `plot_taylor_panels()`; the old
  name is deprecated.
* `plot_reversal_cone()` documentation now states that the display uses a
  stylized bilinear bias model and that the annotated share depends on the
  plotted ranges.
* New vignette, "From a fitted model to a sensitivity report", using the
  public `darfur` data and including a limitations section.

# confoundvis 0.1.0

* First release on CRAN.
