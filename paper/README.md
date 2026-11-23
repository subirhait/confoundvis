# Replication materials

These scripts reproduce every number, table, and figure in

> Hait, S. confoundvis: an R package for sensitivity analysis of risk
> differences and established robustness measures for unmeasured confounding.
> Submitted to *BMC Medical Research Methodology*.

They use **confoundvis 0.4.1**. Install that exact version with

```r
remotes::install_github("subirhait/confoundvis@v0.4.1")
```

or from the archived release (see the Zenodo DOI in the article).

## Scripts

Run them from this folder, in order. Each writes to `output/`.

| Script | What it does | Time |
|---|---|---|
| `01_analysis.R` | All analyses: RHC case study (risk difference, scenarios, covariate comparison, established frameworks), Darfur, Card, 401(k), welders, computing time. Writes `output/tables/`, `output/figures/`, `output/key_numbers.txt`, `output/analysis_log.txt`. | ~3 min |
| `02_figures.R` | Figures 1 and 2 of the article (publication versions). | ~1 min |
| `03_validation.R` | Validation checks: closed-form solution, worked example, additive version, agreement with `marginaleffects`, simulation with a known truth; exports the DoubleML predictions. Stops if any check fails. | ~2 min |
| `04_doubleml_check.py` | Compares `sens_dml()` with the Python package `doubleml` 0.11.4 using the same cross-fitted predictions. | <1 min |

## Data

* **RHC.** Covariates from the R package `ATbounds` (data set `RHC`). The
  30-day death indicator (`dth30`) comes from the Vanderbilt release,
  <https://hbiostat.org/data/repo/rhc.csv>, which the script downloads; or
  place `rhc.csv` in this folder. The script checks that the rows of the two
  sources agree. Note: the `survival` column of `ATbounds::RHC` equals one
  minus the Vanderbilt indicator of death at any time during follow-up
  (`death`), not 30-day survival, and is not used.
* **401(k).** Downloaded by `DoubleML::fetch_401k()`, or place `sipp1991.dta`
  in this folder.
* **Darfur, Card, welders.** In the R packages `sensemakr`, `iv.sensemakr`,
  and `DOS2`.

## Software used for the article

R 4.6.0 on Windows 11; confoundvis 0.4.1, ATbounds 0.1.1, sensemakr 0.1.6,
konfound 1.0.3, EValue 4.1.4, DOS2 0.5.2, iv.sensemakr 0.1.0, DoubleML 1.0.2,
mlr3 1.8.0, mlr3learners 0.16.0, ranger 0.18.0, marginaleffects 1.0.0,
sandwich 3.1-3, ggplot2 4.0.4; Python doubleml 0.11.4. The full session
information is at the end of `output/analysis_log.txt`.
