# =============================================================================
# confoundvis 0.4.1 -- publication figures for the RHC case study
# Figure 1: adjusted risk difference by confounder strength, crossings marked
# Figure 2: surface over prevalence among RHC patients and odds ratio (fine grid)
# Output: output/figures/ in this folder (PDF and 300-dpi PNG)
# Run from this folder (paper/): source("02_figures.R")
# =============================================================================
dl_dir  <- getwd()
fig_dir <- file.path(getwd(), "output", "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
suppressPackageStartupMessages({ library(confoundvis); library(ggplot2) })

# ---- data: ATbounds covariates + Vanderbilt 30-day death --------------------
data("RHC", package = "ATbounds")
rhc <- as.data.frame(RHC)
v <- if (file.exists(file.path(dl_dir, "rhc.csv"))) read.csv(file.path(dl_dir, "rhc.csv")) else
  read.csv("https://hbiostat.org/data/repo/rhc.csv")
stopifnot(all(rhc$RHC == as.integer(v$swang1 == "RHC")), all(abs(rhc$age - v$age) < 1e-6))
rhc$death30 <- as.integer(v$dth30 == "Yes"); rhc$survival <- NULL
fit <- glm(reformulate(setdiff(names(rhc), "death30"), response = "death30"),
           family = binomial, data = rhc)
rd <- sens_risk_difference(fit, "RHC", p1 = 0.4, p0 = 0.2, threshold = 0.02,
                           lambda = seq(1, 9, length.out = 401))
m <- rd$meta
cr <- data.frame(
  or = c(robustness_points(rd)$lambda_sig, m$lambda_threshold, m$lambda_zero),
  label = c("95% CI includes 0", "Estimate reaches 2 per 100", "Estimate reaches 0"))
save2 <- function(p, name, w, h) for (ext in c("pdf", "png"))
  ggsave(file.path(fig_dir, paste0(name, ".", ext)), p, width = w, height = h, dpi = 300)

# ---- Figure 1: path ------------------------------------------------------------
# crossing points: lower 95% limit reaches 0; estimate reaches 0.02; estimate reaches 0
cr$y <- c(0, 0.02, 0)
cr$txt <- c(sprintf("OR %.2f: 95%% CI reaches 0", cr$or[1]),
            sprintf("OR %.2f: 2 per 100", cr$or[2]),
            sprintf("OR %.2f: zero", cr$or[3]))
cr$vj <- c(1.6, -0.8, -0.8)
cr$hj <- c(-0.05, -0.05, -0.05)
f1 <- plot_robustness_curve(rd, points = FALSE) +
  geom_vline(data = cr, aes(xintercept = or), linetype = "dotted", colour = "grey50") +
  geom_point(data = cr, aes(x = or, y = y), size = 2.2) +
  geom_label(data = cr, aes(x = or, y = y, label = txt, vjust = vj, hjust = hj),
             size = 3, linewidth = 0, fill = "white", alpha = 0.85) +
  annotate("text", x = max(rd$lambda), y = 0.02, label = "Practical threshold (2 per 100)",
           hjust = 1, vjust = -0.6, size = 3) +
  scale_y_continuous(labels = function(x) paste0(100 * x),
                     name = "Adjusted risk difference (deaths per 100 patients)") +
  scale_x_continuous(breaks = 1:9,
                     name = "Odds ratio of the hidden confounder with 30-day death") +
  theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank())
save2(f1, "fig1_rhc_riskdiff_path", 7.5, 4.5)

# ---- Figure 2: surface (fine grid) -------------------------------------------
s <- sens_surface(rd, a = seq(0.2, 1, length.out = 161), b = seq(1, 12, length.out = 161))
f2 <- plot_sensitivity_surface(s) +
  geom_vline(xintercept = 0.4, linetype = "dotted", colour = "grey20") +
  annotate("point", x = 0.4, y = c(m$lambda_threshold, m$lambda_zero), size = 2) +
  annotate("text", x = 0.41, y = c(m$lambda_threshold, m$lambda_zero), hjust = 0,
           size = 3, label = sprintf("OR %.2f", c(m$lambda_threshold, m$lambda_zero))) +
  labs(x = "Prevalence of the confounder among RHC patients (20% among others)",
       y = "Odds ratio of the confounder with 30-day death") +
  theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank())
save2(f2, "fig2_rhc_riskdiff_surface", 8, 5.5)
cat("Saved fig1 and fig2 in", fig_dir, "\n")
