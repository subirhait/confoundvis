# confoundvis (development version)

* `sens_risk_difference(confounder_effect = "risk_difference")` no longer
  reports separate arm risks in `meta$risks`. With a confounder effect that
  may vary with the covariates, only the risk difference is determined by
  its average; the arm risks depend on how the effect varies. The risk
  difference, its standard error, and all crossings are unchanged.

# confoundvis 0.4.0

This release includes the changes prepared as 0.3.1, which was not released.

## New frameworks

* `sens_risk_difference()` gives sensitivity analysis on the absolute-risk
  scale for a binary outcome and a binary treatment analysed with logistic
  regression. A hypothetical binary confounder is described by its
  prevalence among the treated and the untreated and its association with
  the outcome (an odds ratio, or a risk difference), as in Lin, Psaty, and
  Kronmal (1998). The function returns the standardized risks and risk
  difference (ATE or ATT) for each confounder strength and the strengths at
  which the risk difference reaches a practical threshold and zero. The
  range across scenarios is reported separately from sampling error
  (delta-method standard errors). Without confounding it equals
  `marginaleffects::avg_comparisons()`, estimate and standard error; the
  risk-difference version equals the correction used by
  `episensr::confounders(type = "RD")`. `confoundvis(glm, risk = list(...))`
  adds it to the RIR. This answers the request for logistic-model support
  in the review of the 0.3.0 article. Scenarios that imply risks outside
  [0, 1] are not used, offsets are kept, the threshold may lie on either
  side of zero, and the report and plot use the level given in `alpha`.
  For the risk-difference version, the confounder-outcome risk difference
  is an average that may vary with the covariates; a value is used only if
  some such pattern keeps every person's risks in [0, 1], and the largest
  such value is the average of the person-level limits. Prevalences of 0
  or 1 are handled.
  The report says so in words when the estimate is already within the
  threshold without confounding, or when no confounder strength reaches
  the threshold or zero (the bias is bounded for given prevalences), and it
  separates a baseline confidence interval that contains the threshold from
  one entirely past it. Surfaces draw thresholds on either side of zero.
* The README quick start now leads with this workflow on simulated data.
* `plot_robustness_curve()` draws bands at the level stored in the path
  (for example 99% for a path built with `alpha = 0.01`), not always 95%.
* The vignette has new sections on risk differences (simulated data),
  instrumental variables, and DoubleML.

* `from_iv_sensemakr()` imports the sensitivity analysis of instrumental
  variable estimates from 'iv.sensemakr' (Cinelli and Hazlett, 2025). The
  path follows the equal-strength line on which 'iv.sensemakr' defines its
  robustness value and reproduces it, together with the extreme robustness
  value and the benchmark bounds. Significance is judged from the
  bias-adjusted Anderson-Rubin confidence limits, which stay valid with a
  weak instrument and can be unbounded.
* `sens_dml()` computes the omitted-variable-bias analysis of Chernozhukov
  et al. (2024) for `DoubleMLPLR` and `DoubleMLIRM` (ATE) fits from the R
  package 'DoubleML', which does not provide it. It follows the Python
  package 'doubleml' (0.11) formula by formula, including its aggregation
  over repeated sample splits, and agrees with it to machine precision for
  the bounds and confidence limits (robustness values to the tolerance of
  the Python optimizer). The fit needs `store_predictions = TRUE`.
* `sens_surface()` and `plot_sensitivity_surface()` support both: the plane
  of the two partial R2 values for instrumental variables (with the
  'iv.sensemakr' benchmark bounds) and the plane of `cf_d` and `cf_y` for
  DoubleML.
* `confoundvis()` accepts 'iv.sensemakr' results and DoubleML fits.
* `plot_robustness_curve()` draws the stored confidence limits of these
  paths instead of +/- 1.96 standard errors.
* `robustness_points()` has a new rule, "adjusted confidence limit", and
  uses a default level stored in the path (`meta$default_alpha`) when
  there is one. `sens_table()` reports RV for instrumental variables and
  RVa for DoubleML.


## Bug fixes

* `robustness_points()`, `sens_report()`, and `sens_table()` no longer
  ignore `alpha` for paths whose threshold was fixed when they were built
  (ITCV, RIR, and Rosenbaum paths). `alpha` now defaults to the stored
  level; asking for a different level is an error that says how to rebuild
  the path. Previously the stored threshold was used silently, and
  `sens_report()` could print a level that was not the one used. The level
  used is now returned in a new `alpha` column.
* `from_konfound()` gains an `alpha` argument (default 0.05, konfound's
  default) so imported paths record their level.
* `rir_lm()` with a substantive `threshold` records no significance level.

## New features

* `confoundvis()` runs everything in one call. Give it a fitted `lm` and the
  treatment name, a `glm` (RIR through 'konfound'), or a vector of
  matched-pair differences (Rosenbaum bounds). It prints a table and a
  plain-language summary, compares the result with the measured
  covariates, and has `plot()` methods for the paths, the covariate
  benchmarks, and the surface.
* `sens_profile()` runs the partial R2, ITCV, and RIR analyses for one
  coefficient in one call, adds any further paths (for example an E-value
  or Rosenbaum path), and returns an object that prints as a side-by-side
  table and plots as one figure with one panel per analysis.
* `covariate_impacts()` gains a `labels` argument for readable covariate
  names, used by `plot_sensitivity_love()`, `plot_sensitivity_contour()`,
  and `plot_sensitivity_surface()`.

## Performance

* `covariate_impacts()` computes every benchmark from one cross-product
  matrix instead of refitting the model for each term. With one million
  observations and 20 covariates it now takes about 5 seconds instead of
  more than 3 minutes. Results agree with the previous version to within
  rounding error; for an exactly aliased term the correlations are now 0
  rather than `NA`.

## Documentation

* Native functions given a `glm` fit now explain how to analyse it with
  'konfound' and `from_konfound()` (RIR for logistic regression) or with
  'EValue' and `from_evalue()`.
* `covariate_impacts()` documents how factors, ordered factors, and
  separately entered dummy variables are benchmarked.

# confoundvis 0.3.0

## Scope

* confoundvis is now a common interface to established sensitivity
  analyses only. The illustration functions `fit_local_quadratic()`,
  `plot_local_taylor()`, `plot_taylor_panels()`,
  `plot_figure2_taylor_panels()`, `simulate_taylor_demo()`,
  `plot_reversal_cone()`, and `plot_cone_comparison()` are removed. They did
  not analyse data and did not correspond to any established framework.
* The introductory vignette built on stylized paths is removed; the
  workflow vignette covers every feature with real data.
* 'grid' is no longer imported and 'gridExtra' is no longer suggested.

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
  *within* each level, so stacked subgroup paths no longer trigger a
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
