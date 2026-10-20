# =============================================================================
# confoundvis 0.3.0: a common interface for sensitivity analyses
# Replication script for the Observational Studies software article
#
# The code under each "Section x.y" heading is copied verbatim from the
# manuscript, in the same order. Lines marked "# [script]" are additions that
# save figures, export numbers, or run checks; they do not appear in the paper.
#
# How to run
#   RStudio: open this file and use "Source with Echo" (Ctrl+Shift+Enter).
#   Terminal: Rscript confoundvis_OBS_replication.R
#
# Outputs (folder confoundvis_output/):
#   printed_output.txt      everything printed, for comparison with the paper
#   figure1.pdf ... figure5.pdf
#   manuscript_numbers.csv, checks.txt, sessionInfo.txt
#
# Required CRAN packages: confoundvis (>= 0.3.0), ATbounds, sensemakr,
# konfound, EValue, DOS2, sensitivity2x2xk, sensitivitymv, sandwich,
# ggplot2, ggrepel, patchwork.
# =============================================================================

required <- c("confoundvis", "ATbounds", "sensemakr", "konfound", "EValue",  # [script]
              "DOS2", "sensitivity2x2xk", "sensitivitymv", "sandwich",       # [script]
              "ggplot2", "ggrepel", "patchwork")                             # [script]
missing <- required[!vapply(required, requireNamespace, logical(1),          # [script]
                            quietly = TRUE)]                                 # [script]
if (length(missing)) stop("Please install: ", paste(missing, collapse = ", "))  # [script]
stopifnot(utils::packageVersion("confoundvis") >= "0.3.0")                   # [script]
out <- "confoundvis_output"                                                  # [script]
dir.create(out, showWarnings = FALSE)                                        # [script]
sink(file.path(out, "printed_output.txt"), split = TRUE)                     # [script]

# ---- Section 4.1: data and starting model -----------------------------------
library(confoundvis)
options(width = 80, digits = 7)
data("RHC", package = "ATbounds")
rhc <- as.data.frame(RHC)
rhc$death <- 1 - rhc$survival
xvars <- setdiff(names(rhc), c("survival", "death"))
fit <- lm(reformulate(xvars, response = "death"), data = rhc)
summary(fit)$coefficients["RHC", ]
c(n = nobs(fit), treated = sum(rhc$RHC), rank = fit$rank,                    # [script]
  df = df.residual(fit))                                                     # [script]

# ---- Section 4.2: partial R2 path -------------------------------------------
grid <- seq(0, 0.12, by = 0.0001)
p_r2 <- sens_path_lm(fit, treatment = "RHC", lambda = grid)
robustness_points(p_r2)
sens_report(p_r2)

# ---- Section 4.3: impact threshold path -------------------------------------
it <- itcv_lm(fit, treatment = "RHC", alpha = 0.05, lambda = grid)
unlist(it[c("r", "r_crit", "itcv", "itcv_unconditional")])
robustness_points(it$path)

# ---- Section 4.4: robustness of inference to replacement --------------------
rr <- rir_lm(fit, treatment = "RHC")
unlist(rr[c("estimate", "beta_threshold", "perc_bias", "rir", "n_obs")])
sens_report(rr$path)

# ---- Section 4.5: importing sensemakr and konfound results ------------------
s <- sensemakr::sensemakr(fit, treatment = "RHC",
                          benchmark_covariates = "surv2md1", kd = 1:3)
p_sm <- from_sensemakr(s, lambda = grid)
k_it <- konfound::konfound(fit, RHC, index = "IT", to_return = "raw_output")
p_kf <- from_konfound(k_it, dof = df.residual(fit), lambda = grid,
                      treatment = "RHC")
k_rir <- konfound::konfound(fit, RHC, index = "RIR", to_return = "raw_output")
p_kr <- from_konfound(k_rir, index = "RIR", treatment = "RHC")

c(sensemakr_theta = max(abs(p_sm$theta - p_r2$theta)),
  konfound_itcv   = max(abs(p_kf$theta - it$path$theta)),
  konfound_rir    = max(abs(p_kr$theta - rr$path$theta)),
  rir_cases       = k_rir$RIR_primary - rr$rir)

# ---- Section 4.6: importing an E-value result -------------------------------
fit_rr <- glm(formula(fit), data = rhc, family = poisson(link = "log"))
V <- sandwich::vcovHC(fit_rr, type = "HC0")
beta <- unname(coef(fit_rr)["RHC"])
se_rr <- sqrt(V["RHC", "RHC"])
ci_rr <- exp(beta + c(-1, 1) * qnorm(0.975) * se_rr)
e <- EValue::evalues.RR(exp(beta), lo = ci_rr[1], hi = ci_rr[2])
p_ev <- from_evalue(e, lambda = seq(1, 2, by = 0.001))
e
robustness_points(p_ev)

# ---- Section 4.7: one plotting function for four frameworks (Figure 1) ------
library(ggplot2)
library(patchwork)
f1 <- plot_robustness_curve(p_r2, points = FALSE) + ggtitle("A. Partial R2") +
      plot_robustness_curve(it$path, points = FALSE) + ggtitle("B. ITCV") +
      plot_robustness_curve(rr$path, points = FALSE) + ggtitle("C. RIR") +
      plot_robustness_curve(p_ev, points = FALSE) + ggtitle("D. E-value") +
      plot_layout(ncol = 2) & theme_bw(base_size = 9)
ggsave(file.path(out, "figure1.pdf"), f1, width = 7.2, height = 6)           # [script]

# ---- Section 4.8: measured-covariate benchmarks (Figures 2 and 3) -----------
ci <- covariate_impacts(fit, treatment = "RHC", metric = "itcv")
cr <- covariate_impacts(fit, treatment = "RHC", metric = "partial_r2")
head(ci[order(-abs(ci$impact)), c("covariate", "r_du", "r_yu", "impact")], 3)
head(cr[order(-cr$impact), c("covariate", "r2dz.x", "r2yz.dx", "impact")], 3)

f2 <- (plot_sensitivity_love(ci, top = 15) +
         labs(title = "A. metric = \"itcv\"", x = "Signed correlation product")) +
      (plot_sensitivity_love(cr, top = 15) +
         labs(title = "B. metric = \"partial_r2\"", x = "Bias index")) &
      theme_bw(base_size = 9)
ggsave(file.path(out, "figure2.pdf"), f2, width = 7.2, height = 4.6)         # [script]

top3 <- head(ci[order(-abs(ci$impact)), ], 3)
f3 <- plot_sensitivity_contour(itcv = it, benchmarks = top3,
                               limits = c(-0.35, 0.35)) +
      theme_bw(base_size = 9)
ggsave(file.path(out, "figure3.pdf"), f3, width = 6, height = 4.6)           # [script]

# ---- Section 4.9: two-parameter surfaces (Figure 4) -------------------------
ax <- seq(0, 0.1, length.out = 201)
s_r2 <- sens_surface(p_sm, a = ax, b = ax)
s_ev <- sens_surface(p_ev)
s_r2
f4 <- plot_sensitivity_surface(s_r2, benchmarks = head(cr[order(-cr$impact), ], 2)) +
        ggtitle("A. Partial R2 (sensemakr result)") +
      plot_sensitivity_surface(s_ev) + ggtitle("B. E-value") +
      plot_layout(ncol = 1) & theme_bw(base_size = 9)
ggsave(file.path(out, "figure4.pdf"), f4, width = 6.5, height = 8)           # [script]

# ---- Section 4.10: several analyses in one table ----------------------------
tb <- sens_table(partial_r2 = p_r2, itcv = it$path, rir = rr$path, evalue = p_ev)
tb[, c("path", "summary_name", "summary_value", "lambda_null", "lambda_sig")]

# ---- Section 5: matched pairs (Figure 5) ------------------------------------
data("werfel", package = "DOS2")
d <- werfel$serpc_p - werfel$cerpc_p
p_gam <- sens_path_pairs(d, gamma = seq(1, 6, by = 0.01))
robustness_points(p_gam)
sens_report(p_gam)
f5 <- plot_robustness_curve(p_gam, points = FALSE) + labs(title = "A. Path", x = "Gamma") +
      plot_sensitivity_surface(sens_surface(p_gam)) + ggtitle("B. Surface") +
      plot_layout(ncol = 2, widths = c(1, 1.35)) & theme_bw(base_size = 9)
ggsave(file.path(out, "figure5.pdf"), f5, width = 7.2, height = 3.3)         # [script]

sink()                                                                       # [script]

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
               "f_threshold", "gamma", "amp_lambda6", "amp_lambda10",
               "surface_r2_b_at_a001", "surface_ev_b_at_a3"),
  value = c(nobs(fit), sum(rhc$RHC), fit$rank, df.residual(fit), cf[1], cf[2],
            p_r2$meta$rv, rp$lambda_sig, it$r, it$r_crit, it$itcv,
            it$itcv_unconditional, rr$rir, rr$perc_bias, rr$beta_threshold,
            exp(beta), ci_rr[1], ci_rr[2],
            e["E-values", "point"], e["E-values", "lower"],
            max(ci$impact, na.rm = TRUE), max(cr$impact),
            attr(ci, "threshold"), attr(cr, "threshold"),
            p_gam$meta$sens_value, amp[1], amp[2], b_at_001, ev_b_at_3))
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
# darfur (Cinelli and Hazlett, 2020) and the practical guide (Frank et al., 2025)
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
  evalue_null        = abs(rpe$lambda_null - e["E-values", "point"]) < 1e-3,
  evalue_ci          = abs(rpe$lambda_sig - e["E-values", "lower"]) < 1e-3,
  surface_r2_diag    = abs(sb_diag - s$sensitivity_stats$rv_qa) < 2e-3,
  surface_ev_null    = max(abs(bf(nb_ev$a, nb_ev$b) - exp(beta))) < 1e-10,
  wilcox_vs_DOS2     = max(abs(ours_wil - wil)) < 1e-12,
  binary_vs_mh       = max(abs(mh_ours - mh_ref)) < 1e-10,
  amplify_vs_smv     = max(abs(amp_ours - amp)) < 1e-3,
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
