## Submission

This is an update from 0.3.0 (accepted 2026-10-17). It fixes a bug: for
ITCV, RIR, and Rosenbaum results, a significance level requested in
robustness_points(), sens_report(), and sens_table() was ignored. It also adds
confoundvis(), a one-call summary, and makes covariate_impacts() much faster.
We release it now because the bug can give wrong reported results.

## R CMD check results

0 errors | 0 warnings | 0 notes

* Local: Windows 11, R 4.6.0
* win-builder: R-devel

## Reverse dependencies

There are no reverse dependencies.
