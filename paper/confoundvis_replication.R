# =============================================================================
# confoundvis: A Common Interface for Visualizing and Reporting Sensitivity
# Analyses for Unmeasured Confounding in R -- replication script (version 0.3.1)
#
# The code under each "Section x.y" heading is the code printed in the
# article, in the same order. Lines marked "# [script]" save figures, export
# numbers, or run checks; they are not printed in the article.
#
# How to run
#   RStudio : open this file and use "Source with Echo" (Ctrl+Shift+Enter).
#   Terminal: Rscript confoundvis_replication.R
# Outputs are written to the folder confoundvis_output/.
#
# Required packages: confoundvis (>= 0.3.1), ATbounds, sensemakr, konfound,
# EValue, DOS2, sensitivity2x2xk, sensitivitymv, sandwich, ggplot2,
# patchwork, ggrepel.
# =============================================================================

required <- c("confoundvis", "ATbounds", "sensemakr", "konfound", "EValue",  # [script]
              "DOS2", "sensitivity2x2xk", "sensitivitymv", "sandwich",       # [script]
              "ggplot2", "patchwork", "ggrepel")                             # [script]
miss <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]  # [script]
if (length(miss)) stop("Please install: ", paste(miss, collapse = ", "))     # [script]
stopifnot(packageVersion("confoundvis") >= "0.3.1")                          # [script]
out <- "confoundvis_output"                                                  # [script]
dir.create(out, showWarnings = FALSE)                                        # [script]
sink(file.path(out, "printed_output.txt"), split = TRUE)                     # [script]

# ---- Section 4.1: data and starting model -----------------------------------
library(confoundvis)
options(width = 80, digits = 7)
data("RHC", package = "ATbounds")              # 5,735 patients (ATbounds package)
rhc <- as.data.frame(RHC)
rhc$death <- 1 - rhc$survival                  # outcome: death within 30 days
xvars <- setdiff(names(rhc), c("survival", "death"))   # RHC and 72 covariates
fit <- lm(reformulate(xvars, response = "death"), data = rhc)  # stats::lm()
summary(fit)$coefficients["RHC", ]             # estimate, SE, t, p for RHC

# ---- Section 4 (shown in "What the package adds"): everything in one call -----
res <- confoundvis(fit, treatment = "RHC",               # one call
                   labels = c(surv2md1 = "Estimated 2-month survival"))
res                                            # table and summary
plot(res)                                      # one panel per analysis

# ---- Section 4.2: one analysis at a time ------------------------------------
grid <- seq(0, 0.12, by = 0.0001)              # confounder strengths to evaluate
p_r2 <- sens_path_lm(fit, treatment = "RHC", lambda = grid)  # partial R2 path
robustness_points(p_r2)                        # where the conclusion changes
sens_report(p_r2)                              # the same result in words

it <- itcv_lm(fit, treatment = "RHC", alpha = 0.05, lambda = grid)  # ITCV
unlist(it[c("r", "r_crit", "itcv", "itcv_unconditional")])
robustness_points(it$path)

rr <- rir_lm(fit, treatment = "RHC")           # RIR at alpha = 0.05
unlist(rr[c("estimate", "beta_threshold", "perc_bias", "rir", "n_obs")])

# ---- Section 4.3: importing sensemakr and konfound results ------------------
s <- sensemakr::sensemakr(fit, treatment = "RHC",          # sensemakr's own call
                          benchmark_covariates = "surv2md1", kd = 1:3)
p_sm <- from_sensemakr(s, lambda = grid)                   # convert to confoundsens
k_it <- konfound::konfound(fit, RHC, index = "IT", to_return = "raw_output")
p_kf <- from_konfound(k_it, dof = df.residual(fit), lambda = grid,
                      treatment = "RHC")
k_rir <- konfound::konfound(fit, RHC, index = "RIR", to_return = "raw_output")
p_kr <- from_konfound(k_rir, index = "RIR", treatment = "RHC")
# largest difference between imported and native paths (0 = identical)
c(sensemakr_theta = max(abs(p_sm$theta - p_r2$theta)),
  konfound_itcv   = max(abs(p_kf$theta - it$path$theta)),
  konfound_rir    = max(abs(p_kr$theta - rr$path$theta)),
  rir_cases       = k_rir$RIR_primary - rr$rir)

# ---- Section 4.4: importing an E-value result -------------------------------
fit_rr <- glm(formula(fit), data = rhc, family = poisson(link = "log"))  # risk ratio
V <- sandwich::vcovHC(fit_rr, type = "HC0")    # robust (HC0) variance
beta <- unname(coef(fit_rr)["RHC"])
se_rr <- sqrt(V["RHC", "RHC"])
ci_rr <- exp(beta + c(-1, 1) * qnorm(0.975) * se_rr)       # 95% interval for RR
e <- EValue::evalues.RR(exp(beta), lo = ci_rr[1], hi = ci_rr[2])
e
p_ev <- from_evalue(e, lambda = seq(1, 2, by = 0.001))     # log-RR path
robustness_points(p_ev)

# ---- Section 4.5: all analyses in one call (Figure 1) -----------------------
prof <- sens_profile(fit, treatment = "RHC", extra = list(evalue = p_ev))
prof                                           # side-by-side table
library(ggplot2)
f1 <- plot(prof) + theme_bw(base_size = 9)     # one panel per analysis
ggsave(file.path(out, "figure1.pdf"), f1, width = 7, height = 5.6)           # [script]

# ---- Section 4.6: measured-covariate benchmarks (Figures 2 and 3) -----------
lab <- c(surv2md1 = "Estimated 2-month survival", wtkilo1 = "Weight (kg)",
         dnr1_Yes = "Do-not-resuscitate order, day 1", wt0 = "Weight missing",
         hrt1 = "Heart rate", pafi1 = "PaO2/FiO2 ratio", age = "Age",
         cat1_Coma = "Primary category: coma", hema_Yes = "Hematologic diagnosis",
         card_Yes = "Cardiovascular diagnosis", resp1 = "Respiratory rate",
         cat1_MOSF_Malignancy = "Primary category: MOSF with malignancy",
         cat2_MOSF_Malignancy = "Secondary category: MOSF with malignancy",
         cat1_MOSF_Sepsis = "Primary category: MOSF with sepsis",
         gastr_Yes = "Gastrointestinal diagnosis",
         cat1_Cirrhosis = "Primary category: cirrhosis",
         neuro_Yes = "Neurological diagnosis")
ci <- covariate_impacts(fit, treatment = "RHC", metric = "itcv", labels = lab)
cr <- covariate_impacts(fit, treatment = "RHC", metric = "partial_r2", labels = lab)
head(ci[order(-abs(ci$impact)), c("label", "r_du", "r_yu", "impact")], 3)
head(cr[order(-cr$impact), c("label", "r2dz.x", "r2yz.dx", "impact")], 3)
library(patchwork)
f2 <- (plot_sensitivity_love(ci, top = 15) +
         labs(title = "A. ITCV metric", x = "Signed correlation product")) +
      (plot_sensitivity_love(cr, top = 15) +
         labs(title = "B. Partial R2 metric", x = "Bias index")) &
      theme_bw(base_size = 9)
ggsave(file.path(out, "figure2.pdf"), f2, width = 7.6, height = 4.4)         # [script]
top3 <- head(ci[order(-abs(ci$impact)), ], 3)
f3 <- plot_sensitivity_contour(itcv = it, benchmarks = top3,
                               limits = c(-0.35, 0.35)) +
      theme_bw(base_size = 9)
ggsave(file.path(out, "figure3.pdf"), f3, width = 5.6, height = 5.2)         # [script]

# ---- Section 4.7: categorical and ordinal covariates ------------------------
# Rebuild the four income categories, stored in RHC as three dummy variables
rhc$income <- factor(ifelse(rhc$income1 == 1, "$11-25k",
                     ifelse(rhc$income2 == 1, "$25-50k",
                     ifelse(rhc$income3 == 1, ">$50k", "<$11k"))),
                     levels = c("<$11k", "$11-25k", "$25-50k", ">$50k"))
rhc$income_score <- as.integer(rhc$income)     # 1 to 4, keeps the order
base <- setdiff(xvars, c("income1", "income2", "income3"))
fit_f <- lm(reformulate(c(base, "income"), response = "death"), data = rhc)
fit_s <- lm(reformulate(c(base, "income_score"), response = "death"), data = rhc)
keep <- c("covariate", "term_df", "r_du", "r_yu", "r2dz.x", "r2yz.dx", "bias_index")
rbind(cr[cr$covariate %in% c("income1", "income2", "income3"), keep],
      subset(covariate_impacts(fit_f, "RHC", metric = "partial_r2"),
             covariate == "income")[, keep],
      subset(covariate_impacts(fit_s, "RHC", metric = "partial_r2"),
             covariate == "income_score")[, keep])

# ---- Section 4.8: two-parameter surfaces (Figure 4) -------------------------
ax <- seq(0, 0.1, length.out = 201)            # axis grid for both partial R2s
s_r2 <- sens_surface(p_sm, a = ax, b = ax)     # from the sensemakr result
s_ev <- sens_surface(p_ev)                     # from the EValue result
s_r2
f4 <- plot_sensitivity_surface(s_r2, benchmarks = head(cr[order(-cr$impact), ], 2)) +
        ggtitle("A. Partial R2 (sensemakr result)") +
      plot_sensitivity_surface(s_ev) + ggtitle("B. E-value") +
      plot_layout(ncol = 1) & theme_bw(base_size = 9)
ggsave(file.path(out, "figure4.pdf"), f4, width = 6.5, height = 8)           # [script]

# ---- Section 4.9: a logistic regression through konfound --------------------
fit_lg <- glm(formula(fit), data = rhc, family = binomial)  # log-odds of death
k_lg <- konfound::konfound(fit_lg, RHC, index = "RIR", to_return = "raw_output")
p_lg <- from_konfound(k_lg, index = "RIR", treatment = "RHC")
sens_report(p_lg)

# ---- Section 5: matched pairs (Figure 5) ------------------------------------
data("werfel", package = "DOS2")               # 39 welder-control pairs
d <- werfel$serpc_p - werfel$cerpc_p           # welder minus matched control
p_gam <- sens_path_pairs(d, gamma = seq(1, 6, by = 0.01))   # one-sided 0.05
robustness_points(p_gam)
sens_report(p_gam)
p_gam2 <- sens_path_pairs(d, gamma = seq(1, 6, by = 0.01), alpha = 0.025)
robustness_points(p_gam2)                      # matches a two-sided 0.05 test
f5 <- plot_robustness_curve(p_gam, points = FALSE) + labs(title = "A. Path", x = "Gamma") +
      plot_sensitivity_surface(sens_surface(p_gam)) + ggtitle("B. Surface") +
      plot_layout(ncol = 2, widths = c(1, 1.35)) & theme_bw(base_size = 9)
ggsave(file.path(out, "figure5.pdf"), f5, width = 7.2, height = 3.3)         # [script]

sink()                                                                       # [script]

# ---- [script] Grid-step study (Table 7) -------------------------------------
gs <- do.call(rbind, lapply(c(0.01, 0.005, 0.001, 1e-4), function(h) {
  g2 <- seq(0, 0.12, by = h)
  pp <- sens_path_lm(fit, "RHC", lambda = g2)
  ii <- itcv_lm(fit, "RHC", lambda = g2)
  data.frame(step = h, points = length(g2),
             err_rv = robustness_points(pp)$lambda_null - pp$meta$rv,
             err_rvqa = robustness_points(pp)$lambda_sig - s$sensitivity_stats$rv_qa,
             err_itcv = robustness_points(ii$path)$lambda_sig - ii$itcv)
}))
print(gs)
write.csv(gs, file.path(out, "grid_step.csv"), row.names = FALSE)
g_short <- robustness_points(sens_path_lm(fit, "RHC", lambda = seq(0, 0.05, by = 1e-4)))
print(g_short)

# ---- [script] Numbers quoted in the manuscript ------------------------------
rp  <- robustness_points(p_r2)
rpe <- robustness_points(p_ev)
cf  <- summary(fit)$coefficients["RHC", ]
amp <- sensitivitymv::amplify(p_gam$meta$sens_value, c(6, 10))
tc1 <- qt(0.975, df.residual(fit) - 1)
tadj <- function(a, b) {
  B  <- cf[2] * sqrt(df.residual(fit)) * sqrt(a * b / (1 - a))
  s2 <- cf[2] * sqrt((1 - b) / (1 - a)) * sqrt(df.residual(fit) / (df.residual(fit) - 1))
  (cf[1] - B) / s2
}
b_at_001 <- uniroot(function(b) tadj(0.01, b) - tc1, c(1e-6, 0.5), tol = 1e-10)$root
ev_b_at_3 <- exp(beta) * (3 - 1) / (3 - exp(beta))
targets <- data.frame(
  quantity = c("n", "treated", "rank", "df", "estimate", "se_lm", "rv",
               "lambda_sig_r2", "r", "r_crit", "itcv", "itcv_uncond",
               "rir", "perc_bias", "beta_threshold",
               "rr", "rr_lo", "rr_hi", "evalue_point", "evalue_ci",
               "top_itcv_impact", "top_bias_index", "itcv_threshold",
               "f_threshold", "gamma", "gamma_alpha025", "amp_lambda6",
               "amp_lambda10", "surface_r2_b_at_a001", "surface_ev_b_at_a3",
               "logit_est", "logit_se", "logit_rir", "logit_perc_bias",
               "logit_beta_threshold"),
  value = c(nobs(fit), sum(rhc$RHC), fit$rank, df.residual(fit), cf[1], cf[2],
            p_r2$meta$rv, rp$lambda_sig, it$r, it$r_crit, it$itcv,
            it$itcv_unconditional, rr$rir, rr$perc_bias, rr$beta_threshold,
            exp(beta), ci_rr[1], ci_rr[2],
            e["E-values", "point"], e["E-values", "lower"],
            max(ci$impact, na.rm = TRUE), max(cr$impact),
            attr(ci, "threshold"), attr(cr, "threshold"),
            p_gam$meta$sens_value, p_gam2$meta$sens_value, amp[1], amp[2],
            b_at_001, ev_b_at_3,
            coef(fit_lg)[["RHC"]], sqrt(vcov(fit_lg)["RHC", "RHC"]),
            k_lg$RIR_primary, k_lg$perc_bias_to_change, k_lg$beta_threshold))
write.csv(targets, file.path(out, "manuscript_numbers.csv"), row.names = FALSE)

# ---- [script] Checks against closed forms and upstream packages -------------
f <- abs(cf[["t value"]]) / sqrt(df.residual(fit))
g_chk <- c(1, 2, 3, 3.5)
wil <- vapply(g_chk, function(G) DOS2::senWilcox(d, gamma = G)$pval, numeric(1))
ours_wil <- stats::pnorm(confoundvis:::.wilcox_dev(d, g_chk), lower.tail = FALSE)
db <- c(rep(1, 40), rep(-1, 22), rep(0, 100))
tab <- array(0, c(2, 2, length(db)))
for (i in seq_along(db)) {
  yt <- as.integer(db[i] == 1); yc <- as.integer(db[i] == -1)
  tab[1, 2 - yt, i] <- 1; tab[2, 2 - yc, i] <- 1
}
mh_ref <- vapply(c(1, 1.2), function(G) sensitivity2x2xk::mh(tab, Gamma = G)$pval,
                 numeric(1))
mh_ours <- sens_path_pairs(db, gamma = c(1, 1.2), method = "binary")$meta$pvalue
data("darfur", package = "sensemakr")
fd <- lm(peacefactor ~ directlyharmed + age + farmer_dar + herder_dar +
           pastvoted + hhsize_darfur + female + village, data = darfur)
sd <- sensemakr::sensemakr(fd, treatment = "directlyharmed")
pd <- sens_path_lm(fd, "directlyharmed", lambda = seq(0, 0.3, by = 1e-4))
itd <- itcv_lm(fd, "directlyharmed")
rrd <- rir_lm(fd, "directlyharmed")
kd_it <- suppressWarnings(konfound::konfound(fd, directlyharmed, index = "IT",
                                             to_return = "raw_output"))
kd_rir <- suppressWarnings(konfound::konfound(fd, directlyharmed, index = "RIR",
                                              to_return = "raw_output"))
kg <- suppressWarnings(suppressMessages(konfound::pkonfound(
  0.1, 0.03, 6320, 5, index = "RIR", to_return = "raw_output")))
pg_guide <- from_konfound(kg, index = "RIR")
rg <- confoundvis:::.rir_core(0.1, qt(0.975, 6320 - 5 - 2) * 0.03, 6320,
                              seq(0, 1, by = 0.001), 0.05, "summary", NA)
sb <- s_r2$boundary[grepl("^Significance \\(alpha", s_r2$boundary$kind), ]
sb_diag <- sb$a[which.min(abs(sb$a - sb$b))]
nb_ev <- s_ev$boundary[grepl("^Null", s_ev$boundary$kind), ]
bf <- function(a, b) a * b / (a + b - 1)
sbp <- sens_surface(p_gam)$boundary
amp_ours <- stats::approx(sbp$a, sbp$b, xout = c(6, 10))$y
alpha_err <- inherits(try(robustness_points(it$path, alpha = 0.01), silent = TRUE),
                      "try-error")
it01 <- itcv_lm(fit, "RHC", alpha = 0.01, lambda = grid)

checks <- c(
  rv_closed_form     = abs(p_r2$meta$rv - (sqrt(f^4 + 4 * f^2) - f^2) / 2) < 1e-12,
  rv_vs_sensemakr    = abs(p_r2$meta$rv - s$sensitivity_stats$rv_q) < 1e-12,
  rvqa_vs_grid       = abs(rp$lambda_sig - s$sensitivity_stats$rv_qa) < 1e-6,
  itcv_vs_konfound   = abs(it$itcv - k_it$itcvGz) < 1e-12,
  itcvu_vs_konfound  = abs(it$itcv_unconditional - k_it$itcv) < 1e-12,
  rir_vs_konfound    = rr$rir == k_rir$RIR_primary &&
                       abs(rr$perc_bias - k_rir$perc_bias_to_change) < 1e-8,
  paths_identical    = max(abs(p_sm$theta - p_r2$theta),
                           abs(p_kf$theta - it$path$theta),
                           abs(p_kr$theta - rr$path$theta)) == 0,
  profile_matches    = isTRUE(all.equal(prof$table$summary_value[1:3],
                                        c(p_r2$meta$rv, it$itcv, rr$rir))),
  onecall_matches    = isTRUE(all.equal(res$table$summary_value,
                                        c(p_r2$meta$rv, it$itcv, rr$rir))),
  evalue_null        = abs(rpe$lambda_null - e["E-values", "point"]) < 1e-3,
  evalue_ci          = abs(rpe$lambda_sig - e["E-values", "lower"]) < 1e-3,
  surface_r2_diag    = abs(sb_diag - s$sensitivity_stats$rv_qa) < 2e-3,
  surface_ev_null    = max(abs(bf(nb_ev$a, nb_ev$b) - exp(beta))) < 1e-10,
  wilcox_vs_DOS2     = max(abs(ours_wil - wil)) < 1e-12,
  binary_vs_mh       = max(abs(mh_ours - mh_ref)) < 1e-10,
  amplify_vs_smv     = max(abs(amp_ours - amp)) < 1e-3,
  logit_vs_konfound  = p_lg$meta$rir == k_lg$RIR_primary &&
                       abs(robustness_points(p_lg)$lambda_sig -
                           k_lg$perc_bias_to_change / 100) < 1e-3,
  alpha_refused      = alpha_err,
  alpha_rebuilt      = abs(robustness_points(it01$path)$lambda_sig - it01$itcv) < 1e-6,
  darfur_rv          = abs(pd$meta$rv - sd$sensitivity_stats$rv_q) < 1e-12,
  darfur_rvqa        = abs(robustness_points(pd)$lambda_sig -
                             sd$sensitivity_stats$rv_qa) < 1e-6,
  darfur_itcv        = abs(itd$itcv - kd_it$itcvGz) < 1e-12,
  darfur_rir         = rrd$rir == kd_rir$RIR_primary && rrd$rir == 677,
  guide_rir          = pg_guide$meta$rir == 2603 && rg$meta$rir == 2603 &&
                       round(rg$meta$perc_bias, 2) == 41.19)
print(checks)
writeLines(paste(names(checks), checks), file.path(out, "checks.txt"))
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
if (!all(checks)) stop("One or more checks failed: ",
                       paste(names(checks)[!checks], collapse = ", "))
message("All checks passed. Outputs saved to ", normalizePath(out))
