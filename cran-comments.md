## Update

This is an update from confoundvis 0.2.0 to 0.3.0.

* New: Rosenbaum bounds for matched pairs (`sens_path_pairs()`), robustness
  of inference to replacement (`rir_lm()`), two-parameter sensitivity
  surfaces (`sens_surface()`, `plot_sensitivity_surface()`), a summary
  table (`sens_table()`), and an exact ITCV contour.
* Removed exported functions: `fit_local_quadratic()`, `plot_local_taylor()`,
  `plot_taylor_panels()`, `plot_figure2_taylor_panels()`,
  `simulate_taylor_demo()`, `plot_reversal_cone()`, and
  `plot_cone_comparison()`. These were stylized illustrations that did not
  analyse data. The package has no reverse dependencies on CRAN.
* 'grid' was removed from Imports and 'gridExtra' from Suggests.
* All new Suggests packages ('DOS2', 'sensitivity2x2xk', 'sensitivitymv',
  'ggrepel') are used conditionally in tests, examples, and the vignette.

## Test environments

* local Windows 11, R 4.6.0
* win-builder R-devel and R-release
* GitHub Actions: ubuntu-latest (devel, release, oldrel-1), macos-latest,
  windows-latest

## R CMD check results

0 errors | 0 warnings | 0 notes

## Notes on possibly flagged items

* "Possibly misspelled words": Cinelli, Hazlett, VanderWeele, Rosenbaum,
  and konfound are author surnames or package names, not misspellings.
