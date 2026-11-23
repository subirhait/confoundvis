## Submission

This is an update of confoundvis (0.3.0 on CRAN) to 0.4.0.

* Fixes a bug in which the significance level was ignored for paths with a
  stored threshold (ITCV, RIR, Rosenbaum).
* Adds confoundvis(), a one-call interface, and sens_profile().
* Version 0.3.1 was prepared but not submitted; its changes are included.
* Adds from_iv_sensemakr() (instrumental-variable sensitivity, imported from
  'iv.sensemakr'), sens_dml() (omitted-variable-bias bounds for fits from
  the 'DoubleML' package), and sens_risk_difference() (risk-difference
  sensitivity for binary outcomes).

New packages in Suggests: iv.sensemakr, DoubleML, mlr3, mlr3learners,
data.table, lgr, marginaleffects. All examples and tests that use them are
conditional on their being installed.

## Test environments

* local Windows 11, R 4.6.0
* win-builder (R-devel)

## R CMD check results

0 errors | 0 warnings | 1 note

* The NOTE lists possibly misspelled words in DESCRIPTION: Chernozhukov,
  Psaty, and Kronmal are author names in the references, and confounder is
  a correct technical term.

## Reverse dependencies

There are no reverse dependencies.

