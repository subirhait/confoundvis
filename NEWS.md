# confoundvis 0.3.0

## New features

* `sens_path_pairs()` computes a sensitivity path for matched-pair studies
  under Rosenbaum's sensitivity model. For each value of Gamma it stores
  the deviate equivalent to the upper bound on the one-sided P-value, so
  `robustness_points()`, `plot_robustness_curve()`, `sens_report()`, and
  `sens_table()` apply unchanged. Two statistics are available: Wilcoxon's
  signed rank statistic (`method = "wilcoxon"`, agreeing with
  `DOS2::senWilcox()`) and the exact bound for McNemar's statistic for
  binary outcomes (`method = "binary"`, agreeing with
  `sensitivity2x2xk::mh()`). The sensitivity value is also found by root
  finding and stored in `meta$sens_value`.
* `rir_lm()` computes the robustness of inference to replacement (RIR;
  Frank et al., 2013) from a fitted linear model: the threshold for
  significance, the percent of the estimate that would have to be due to
  bias, the number of cases that would have to be replaced with null-effect
  cases, and the corresponding path (estimate after replacing a proportion
  lambda of the cases). Results agree with `konfound::konfound(index =
  "RIR")` and reproduce the worked example of Frank et al.'s (2025)
  practical guide. A `threshold` argument allows a substantive effect
  threshold in place of statistical significance. `from_konfound()` gains
  `index = "RIR"` to import konfound's RIR output as a path.
* ITCV benchmarks following Frank et al.'s (2025) practical guide:
  `itcv_lm()` returns `benchmark_corr_product` (R_XZ * R_YZ, the impact of
  all covariates together) and `itcv_ratio_to_benchmark` (eq. 12), and
  `covariate_impacts()` adds each term's zero-order correlations, their
  product `impact_zero`, and `itcv_ratio`, the unconditional ITCV divided by
  that impact (eq. 10).
* `sens_surface()` returns the two-parameter sensitivity surface of which a
  path is one slice, and `plot_sensitivity_surface()` draws it with shaded
  zones, the null and significance boundaries, and optional benchmarks
  (sensemakr formal bounds are added automatically). Surfaces are available
  for partial R2 (Cinelli and Hazlett, 2020), ITCV (exact boundaries),
  E-values (joint bounding factor of Ding and VanderWeele, 2016), and
  Rosenbaum bounds (amplification of Rosenbaum and Silber, 2009, agreeing
  with `sensitivitymv::amplify()`). Boundaries are computed in closed form
  or by one-dimensional root finding.
* `sens_table()` collects several `confoundsens` objects into one table:
  framework, meaning of the sensitivity parameter, the framework's own
  summary statistic, null and significance crossings, and the rule used.
  Rows keep their own scales; the printed table says so.
* `plot_sensitivity_contour()` gains an `itcv` argument. Given the output of
  `itcv_lm()`, it draws, in closed form, the exact boundaries at which the
  magnitude of the confounder-adjusted partial correlation equals the
  critical value: in the original direction (solid) and in the reversed
  direction (dotted). It shades the region where the estimate is no longer
  significant, the region where it becomes significant with the opposite
  sign, and correlation pairs that cannot occur (non-positive-definite
  correlation matrix). The product curve is kept as a dashed reference.
  `threshold` may now be omitted when `itcv` is supplied, and a new
  `limits` argument sets the axis range.
* Benchmark labels in `plot_sensitivity_contour()` repel one another when
  **ggrepel** is installed (suggested), instead of overlapping.

## Bug fixes

* `itcv_lm()`: for weighted fits, the R-squared values used for the
  unconditional ITCV were computed around an unweighted mean. They now
  equal the weighted R-squared reported by `summary.lm()`. Unweighted fits
  are unaffected.
* `sens_report()` no longer points to `covariate_impacts()` for frameworks
  that it does not benchmark (E-values, Rosenbaum bounds).
* `robustness_points()` now takes the direction of the original conclusion
  from the stored threshold when there is one. Paths that start on the
  wrong side of their threshold are therefore reported as not significant
  at the first grid value instead of returning `NA`. Results for paths
  produced by earlier versions are unchanged.

## Other changes

* **DOS2** and **sensitivity2x2xk** are suggested for the cross-checks in
  the test suite.
* The workflow vignette adds the exact contour, a matched-pair example, and
  `sens_table()`.

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
