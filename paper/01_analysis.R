# =============================================================================
# confoundvis 0.4.0 -- analyses for the manuscript (BMC Medical Research Methodology)
#
# Produces every number, table, and figure for the revised article:
#   1. Main case study: right heart catheterization (RHC) and 30-day death
#      (Connors et al., 1996), with the verified Vanderbilt `dth30` outcome.
#      - risk-difference sensitivity under several confounder scenarios
#      - calibration against measured binary covariates
#      - the established analyses on the same data, each in its own units
#        (partial R2, ITCV, RIR on the linear probability scale; logistic RIR;
#        E-values on the risk-ratio scale)
#   2. Supporting examples: Darfur (linear model), Card (instrumental
#      variables), 401(k) (double machine learning), welders (matched pairs).
#   3. Computing time for every analysis.
#
# Output: output/ in this folder
#   tables/*.csv, figures/*.pdf and *.png, key_numbers.txt, analysis_log.txt
#
# Run from this folder (paper/):  source("01_analysis.R")
# Needs internet for the Vanderbilt RHC file and the 401(k) data, or local
# copies (rhc.csv, sipp1991.dta) in this folder. About 3 minutes.
# =============================================================================

dl_dir  <- getwd()                       # local data files, if any
out_dir <- file.path(getwd(), "output")
for (d in c("tables", "figures")) dir.create(file.path(out_dir, d), recursive = TRUE,
                                             showWarnings = FALSE)

need <- c("confoundvis", "ATbounds", "sensemakr", "konfound", "EValue", "sandwich",
          "marginaleffects", "iv.sensemakr", "DoubleML", "mlr3", "mlr3learners",
          "ranger", "readstata13", "data.table", "lgr", "DOS2", "ggplot2")
miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) install.packages(miss)
suppressPackageStartupMessages(library(confoundvis))
if (utils::packageVersion("confoundvis") < "0.4.0")
  stop("This script needs confoundvis 0.4.0: ",
       "remotes::install_github('subirhait/confoundvis@v0.4.0')")

log_file <- file.path(out_dir, "analysis_log.txt")
sink(log_file, split = TRUE)
cat("confoundvis", as.character(utils::packageVersion("confoundvis")),
    "| run", format(Sys.time()), "\n")

timing <- numeric(0)
tm <- function(name, expr) {
  t0 <- proc.time()[["elapsed"]]; v <- force(expr)
  timing[name] <<- proc.time()[["elapsed"]] - t0; v
}
key <- character(0)
note <- function(...) { s <- sprintf(...); key <<- c(key, s); cat("  ", s, "\n") }
fig <- function(p, name, w = 7, h = 4.5) {
  for (ext in c("pdf", "png"))
    ggplot2::ggsave(file.path(out_dir, "figures", paste0(name, ".", ext)), p,
                    width = w, height = h, dpi = 300)
}
tab <- function(x, name) utils::write.csv(x, file.path(out_dir, "tables", paste0(name, ".csv")),
                                          row.names = FALSE)
quiet <- function(expr) suppressWarnings(suppressMessages(expr))
r3 <- function(x) round(x, 3)

# =============================================================================
cat("\n== 1. RHC and 30-day death ==\n")
data("RHC", package = "ATbounds")
rhc <- as.data.frame(RHC)
v <- NULL
if (file.exists(file.path(dl_dir, "rhc.csv"))) v <- read.csv(file.path(dl_dir, "rhc.csv"))
if (is.null(v)) v <- tryCatch(read.csv("https://hbiostat.org/data/repo/rhc.csv"),
                              error = function(e) NULL)
if (is.null(v)) stop("Vanderbilt rhc.csv not available: download it from ",
                     "https://hbiostat.org/data/repo/rhc.csv into ", dl_dir)
stopifnot(nrow(v) == nrow(rhc),
          all(rhc$RHC == as.integer(v$swang1 == "RHC")),
          all(abs(rhc$age - v$age) < 1e-6),
          all(1 - rhc$survival == as.integer(v$death == "Yes")))
cat("  Rows aligned with Vanderbilt; ATbounds `survival` = 1 - Vanderbilt `death`.\n")
rhc$death30 <- as.integer(v$dth30 == "Yes")
rhc$survival <- NULL
covs <- setdiff(names(rhc), c("RHC", "death30"))
note("n = %d; RHC = %d; 30-day deaths = %d (%.1f%%); follow-up deaths = %.1f%%",
     nrow(rhc), sum(rhc$RHC), sum(rhc$death30), 100 * mean(rhc$death30),
     100 * mean(v$death == "Yes"))
note("Crude 30-day death: RHC %.1f%%, no RHC %.1f%%",
     100 * mean(rhc$death30[rhc$RHC == 1]), 100 * mean(rhc$death30[rhc$RHC == 0]))
f_rhc <- reformulate(c("RHC", covs), response = "death30")

# ---- 1a. Main risk-difference analysis -------------------------------------
fit_b <- glm(f_rhc, family = binomial, data = rhc)
rd <- tm("RHC: sens_risk_difference (main)",
         sens_risk_difference(fit_b, "RHC", p1 = 0.4, p0 = 0.2, threshold = 0.02))
me <- marginaleffects::avg_comparisons(fit_b, variables = "RHC")
m <- rd$meta
note("Standardized 30-day risks: %.1f%% (RHC) vs %.1f%% (no RHC)",
     100 * m$risk1_obs, 100 * m$risk0_obs)
note("Risk difference %.4f (SE %.4f; 95%% CI %.4f to %.4f); marginaleffects %.4f (SE %.4f)",
     m$estimate, m$se0, m$estimate - 1.96 * m$se0, m$estimate + 1.96 * m$se0,
     me$estimate, me$std.error)
note("Scenario p1 = 0.40, p0 = 0.20: OR %.2f reaches the 0.02 threshold; OR %.2f reaches zero",
     m$lambda_threshold, m$lambda_zero)
rp <- robustness_points(rd)
note("Significance (alpha = 0.05) is lost at OR %.2f; the confidence limit reaches 0.02 at OR %.2f",
     rp$lambda_sig, m$lambda_threshold_ci)
print(sens_report(rd))
fig(plot_robustness_curve(rd, points = TRUE), "fig_rhc_riskdiff_path")
s_rd <- tm("RHC: risk-difference surface", sens_surface(rd))
fig(plot_sensitivity_surface(s_rd), "fig_rhc_riskdiff_surface", 8, 5.5)
tab(data.frame(odds_ratio = rd$lambda, risk_rhc = m$risks$risk1,
               risk_no_rhc = m$risks$risk0, risk_difference = rd$theta, se = rd$se),
    "rhc_riskdiff_path")

# ---- 1b. Scenario table -----------------------------------------------------
sc <- data.frame(p1 = c(0.3, 0.4, 0.5, 0.4, 0.6, 0.4, 0.4),
                 p0 = c(0.2, 0.2, 0.2, 0.1, 0.3, 0.2, 0.2),
                 target = c(rep("ATE", 5), "ATT", "ATE"),
                 effect = c(rep("odds_ratio", 6), "risk_difference"))
rows <- lapply(seq_len(nrow(sc)), function(i) {
  r <- quiet(sens_risk_difference(fit_b, "RHC", p1 = sc$p1[i], p0 = sc$p0[i],
                                  threshold = 0.02, target = sc$target[i],
                                  confounder_effect = sc$effect[i]))
  data.frame(p1 = sc$p1[i], p0 = sc$p0[i], target = sc$target[i],
             confounder_effect = sc$effect[i], estimate = r3(r$meta$estimate),
             to_threshold = r3(r$meta$lambda_threshold),
             to_zero = r3(r$meta$lambda_zero),
             sig_lost = r3(robustness_points(r)$lambda_sig),
             ci_to_threshold = r3(r$meta$lambda_threshold_ci),
             max_feasible = if (sc$effect[i] == "risk_difference") r3(r$meta$max_feasible) else NA)
})
scen <- do.call(rbind, rows)
cat("\n  Scenario table (strengths are odds ratios, or risk differences for the additive row):\n")
print(scen, row.names = FALSE)
tab(scen, "rhc_scenarios")

# ---- 1c. Calibration against measured binary covariates --------------------
# For each binary covariate W: its prevalence among RHC and non-RHC patients
# and its adjusted odds ratio with 30-day death in the outcome model. The
# adjusted risk difference is what a hidden confounder with the same
# prevalences and the same odds ratio would leave.
bin <- covs[vapply(rhc[covs], function(z) all(z %in% c(0, 1)), logical(1))]
b_or <- exp(coef(fit_b)[bin])
cal <- do.call(rbind, lapply(bin, function(w) {
  p1w <- mean(rhc[[w]][rhc$RHC == 1]); p0w <- mean(rhc[[w]][rhc$RHC == 0])
  orw <- b_or[[w]]
  if (is.na(orw) || abs(p1w - p0w) < 0.01) return(NULL)
  r0 <- quiet(sens_risk_difference(fit_b, "RHC", p1 = p1w, p0 = p0w, lambda = c(1, 2)))
  lam <- orw^r0$meta$direction                  # strength in the direction toward zero
  adj <- if (lam >= 1) quiet(sens_risk_difference(fit_b, "RHC", p1 = p1w, p0 = p0w,
                                                  lambda = c(1, lam)))$theta[2]
         else NA_real_
  data.frame(covariate = w, prev_rhc = r3(p1w), prev_no_rhc = r3(p0w),
             adjusted_or = r3(orw), toward_zero = lam >= 1,
             rd_if_hidden = r3(adj), rd_if_hidden_2x = if (lam >= 1) r3(quiet(
               sens_risk_difference(fit_b, "RHC", p1 = p1w, p0 = p0w,
                                    lambda = c(1, lam^2)))$theta[2]) else NA)
}))
cal <- cal[order(cal$rd_if_hidden, na.last = TRUE), ]
cat("\n  Calibration: risk difference left by a hidden confounder like each covariate\n")
print(utils::head(cal, 10), row.names = FALSE)
tab(cal, "rhc_covariate_calibration")
note("Strongest measured covariate (toward zero): %s leaves RD %.4f (2x its log-OR: %.4f)",
     cal$covariate[1], cal$rd_if_hidden[1], cal$rd_if_hidden_2x[1])

# ---- 1d. Established analyses on the same data, each in its own units ------
cat("\n  -- Linear probability model: partial R2, ITCV, RIR --\n")
fit_lp <- lm(f_rhc, data = rhc)
note("Linear probability coefficient %.5f (SE %.5f)", coef(fit_lp)[["RHC"]],
     sqrt(vcov(fit_lp)["RHC", "RHC"]))
cv_lp <- tm("RHC: confoundvis(lm)", quiet(confoundvis(fit_lp, "RHC")))
print(cv_lp)
tab(as.data.frame(cv_lp$table), "rhc_lpm_summary")
fig(plot(cv_lp), "fig_rhc_lpm_paths", 9, 4)

cat("\n  -- Logistic RIR (konfound) --\n")
cv_b <- tm("RHC: confoundvis(glm)", quiet(confoundvis(fit_b, "RHC",
               risk = list(p1 = 0.4, p0 = 0.2, threshold = 0.02))))
print(cv_b)
note("Logistic RIR: %s cases (%.1f%% of the estimate)",
     format(cv_b$paths$rir$meta$rir, big.mark = ","),
     100 * robustness_points(cv_b$paths$rir)$lambda_sig)

cat("\n  -- Risk ratio and E-values (robust HC0 SE) --\n")
fit_rr <- quiet(glm(f_rhc, data = rhc, family = poisson(link = "log")))
V <- sandwich::vcovHC(fit_rr, type = "HC0")
b <- coef(fit_rr)[["RHC"]]; se <- sqrt(V["RHC", "RHC"])
ci <- exp(b + c(-1, 1) * qnorm(0.975) * se)
ev <- EValue::evalues.RR(exp(b), lo = ci[1], hi = ci[2])
p_ev <- tm("RHC: E-value path", from_evalue(ev, lambda = seq(1, 3, by = 1e-4)))
note("Risk ratio %.3f (95%% CI %.3f to %.3f); E-value %.3f (point), %.3f (CI)",
     exp(b), ci[1], ci[2], ev["E-values", "point"], ev["E-values", "lower"])
tab(sens_table(riskdiff = rd, rir = cv_b$paths$rir, evalue = p_ev,
               partial_r2 = cv_lp$paths$partial_r2, itcv = cv_lp$paths$itcv), "rhc_all_frameworks")

# =============================================================================
cat("\n== 2. Supporting examples ==\n")
cat("\n  -- Darfur (linear model) --\n")
data("darfur", package = "sensemakr")
fit_d <- lm(peacefactor ~ directlyharmed + age + farmer_dar + herder_dar +
              pastvoted + hhsize_darfur + female + village, data = darfur)
cv_d <- tm("Darfur: confoundvis(lm)", quiet(confoundvis(fit_d, "directlyharmed")))
print(cv_d); tab(as.data.frame(cv_d$table), "darfur_summary")
fig(plot(cv_d), "fig_darfur_paths", 9, 4)

cat("\n  -- Card (instrumental variables) --\n")
data("card", package = "iv.sensemakr")
xm <- model.matrix(~ exper + expersq + black + south + smsa + reg661 + reg662 +
                     reg663 + reg664 + reg665 + reg666 + reg667 + reg668 + smsa66,
                   data = card)
ivs <- iv.sensemakr::sensemakr(iv.sensemakr::iv_fit(card$lwage, card$educ, card$nearc4, xm),
                               benchmark_covariates = c("black", "smsa"))
cv_iv <- tm("Card: confoundvis(iv)", confoundvis(ivs))
print(cv_iv); tab(as.data.frame(cv_iv$table), "card_iv_summary")
fig(plot(cv_iv, which = "surface"), "fig_card_iv_surface", 7, 5.5)
note("Card IV: estimate %.4f; RV %.5f; XRV %.6f", cv_iv$paths$iv$meta$estimate,
     cv_iv$paths$iv$meta$rv, cv_iv$paths$iv$meta$xrv)

cat("\n  -- 401(k) (double machine learning) --\n")
k401 <- NULL
if (file.exists(file.path(dl_dir, "sipp1991.dta")))
  k401 <- data.table::as.data.table(readstata13::read.dta13(file.path(dl_dir, "sipp1991.dta")))
if (is.null(k401)) k401 <- tryCatch(quiet(DoubleML::fetch_401k(return_type = "data.table")),
                                    error = function(e) NULL)
if (!is.null(k401)) {
  lgr::get_logger("mlr3")$set_threshold("warn")
  dd <- DoubleML::DoubleMLData$new(k401, y_col = "net_tfa", d_cols = "e401",
          x_cols = c("age", "inc", "educ", "fsize", "marr", "twoearn", "db", "pira", "hown"))
  set.seed(2026)
  plr <- DoubleML::DoubleMLPLR$new(dd, ml_l = mlr3::lrn("regr.ranger", num.trees = 500),
            ml_m = mlr3::lrn("regr.ranger", num.trees = 500), n_folds = 5, n_rep = 3)
  tm("401(k): DoubleML PLR fit", plr$fit(store_predictions = TRUE))
  set.seed(2026)
  irm <- DoubleML::DoubleMLIRM$new(dd, ml_g = mlr3::lrn("regr.ranger", num.trees = 500),
            ml_m = mlr3::lrn("classif.ranger", num.trees = 500), n_folds = 5, n_rep = 3,
            trimming_threshold = 0.01)
  tm("401(k): DoubleML IRM fit", irm$fit(store_predictions = TRUE))
  p_plr <- tm("401(k): sens_dml PLR", sens_dml(plr))
  p_irm <- tm("401(k): sens_dml IRM", sens_dml(irm))
  dml_tab <- sens_table(PLR = p_plr, IRM = p_irm)
  print(dml_tab); tab(as.data.frame(dml_tab), "k401_dml_summary")
  note("401(k) PLR: %.0f; RV %.4f; RVa %.4f | IRM: %.0f; RV %.4f; RVa %.4f",
       p_plr$meta$estimate, p_plr$meta$rv, p_plr$meta$rva,
       p_irm$meta$estimate, p_irm$meta$rv, p_irm$meta$rva)
  fig(plot_sensitivity_surface(tm("401(k): DML surface", sens_surface(p_plr))),
      "fig_401k_plr_surface", 7, 5.5)
} else cat("  401(k) data not available: skipped\n")

cat("\n  -- Welders (matched pairs) --\n")
data("werfel", package = "DOS2")
cv_p <- tm("Welders: confoundvis(pairs)", confoundvis(werfel$serpc_p - werfel$cerpc_p))
print(cv_p); tab(as.data.frame(cv_p$table), "welders_summary")
note("Welders: sensitivity value Gamma %.3f", cv_p$paths$rosenbaum$meta$sens_value)

# =============================================================================
cat("\n== 3. Computing time (seconds) ==\n")
tt <- data.frame(analysis = names(timing), seconds = round(unname(timing), 2))
print(tt, row.names = FALSE); tab(tt, "computing_time")

writeLines(c(sprintf("confoundvis %s, %s, %s", utils::packageVersion("confoundvis"),
                     R.version.string, format(Sys.time())), "", key),
           file.path(out_dir, "key_numbers.txt"))
cat("\nKey numbers written to", file.path(out_dir, "key_numbers.txt"), "\n")
print(sessionInfo())
sink()
cat("Done. Output in", out_dir, "\n")
