# =============================================================================
# confoundvis 0.4.1 -- validation checks reported in the manuscript
#
#   1. Closed-form stratum risks agree with numerical root finding.
#   2. Worked example: risks 16% vs 12%, confounder in 40% vs 20%.
#   3. Additive version equals RD - delta (p1 - p0).
#   4. Without confounding, the estimate and SE equal marginaleffects.
#   5. Simulation with a known true risk difference.
#   6. Double machine learning: exports the cross-fitted predictions of the
#      R DoubleML fits so that 04_doubleml_check.py can compare the bounds
#      with the Python package doubleml (pip install doubleml==0.11.4).
#
# Run from this folder (paper/): source("03_validation.R"). About 2 minutes.
# Each check prints PASS or FAIL; the script stops if any check fails.
# =============================================================================
suppressPackageStartupMessages(library(confoundvis))
if (utils::packageVersion("confoundvis") < "0.4.1")
  stop("remotes::install_github('subirhait/confoundvis@v0.4.1')")
out_dir <- file.path(getwd(), "output", "validation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
ok <- logical(0)
check <- function(name, cond, detail = "") {
  ok[name] <<- isTRUE(cond)
  cat(sprintf("  %-62s %s %s\n", name, if (isTRUE(cond)) "PASS" else "FAIL", detail))
}

cat("\n1. Closed-form root of the mixture equation\n")
set.seed(1)
r <- runif(200, 0.01, 0.99); p <- runif(200); G <- exp(runif(200, -3, 3))
t_closed <- confoundvis:::.rd_t(r, p, G)
t_root <- mapply(function(ri, pi, gi) exp(uniroot(function(h)
  pi * plogis(h + log(gi)) + (1 - pi) * plogis(h) - ri, c(-50, 50), tol = 1e-14)$root),
  r, p, G)
check("closed form vs uniroot (200 random cases)",
      max(abs(t_closed - t_root) / t_root) < 1e-8,
      sprintf("max rel. diff %.1e", max(abs(t_closed - t_root) / t_root)))

cat("\n2. Worked example (risks 16% vs 12%; U in 40% vs 20%)\n")
toy <- data.frame(a = rep(0:1, each = 10000),
                  y = c(rep(1:0, c(1200, 8800)), rep(1:0, c(1600, 8400))))
ft <- glm(y ~ a, family = binomial, data = toy)
rt <- sens_risk_difference(ft, "a", p1 = 0.4, p0 = 0.2, threshold = 0.02)
rd_hand <- function(g) {                 # direct calculation, no package code
  h <- function(r, p) uniroot(function(z) p * plogis(z + log(g)) + (1 - p) * plogis(z) - r,
                              c(-30, 30), tol = 1e-14)$root
  h1 <- h(0.16, 0.4); h0 <- h(0.12, 0.2); q <- 0.3
  (q * plogis(h1 + log(g)) + (1 - q) * plogis(h1)) -
    (q * plogis(h0 + log(g)) + (1 - q) * plogis(h0))
}
check("odds ratio to reach 0.02 (package vs direct)",
      abs(rd_hand(rt$meta$lambda_threshold) - 0.02) < 1e-7,
      sprintf("OR %.3f", rt$meta$lambda_threshold))
check("odds ratio to reach zero (package vs direct)",
      abs(rd_hand(rt$meta$lambda_zero)) < 1e-7, sprintf("OR %.3f", rt$meta$lambda_zero))

cat("\n3. Additive version\n")
ra <- sens_risk_difference(ft, "a", p1 = 0.4, p0 = 0.2,
                           confounder_effect = "risk_difference")
check("RD - delta (p1 - p0)", max(abs(ra$theta - (0.04 - 0.2 * ra$lambda))) < 1e-10)

cat("\n4. No confounding vs marginaleffects\n")
set.seed(3); n <- 3000
d <- data.frame(x1 = rnorm(n), x2 = factor(sample(letters[1:3], n, TRUE)))
d$a <- rbinom(n, 1, plogis(0.4 * d$x1))
d$y <- rbinom(n, 1, plogis(-1 + 0.5 * d$a + 0.7 * d$x1 + 0.2 * d$a * d$x1))
f <- glm(y ~ a * x1 + x2, family = binomial, data = d)
for (tg in c("ATE", "ATT")) {
  r0 <- sens_risk_difference(f, "a", p1 = 0.5, p0 = 0.3, target = tg)
  me <- marginaleffects::avg_comparisons(f, variables = "a",
          newdata = if (tg == "ATT") subset(d, a == 1) else d)
  check(sprintf("%s estimate and SE vs marginaleffects", tg),
        abs(r0$theta[1] - me$estimate) < 1e-8 && abs(r0$se[1] - me$std.error) < 1e-5)
}

cat("\n5. Simulation with a known true risk difference (n = 400,000)\n")
set.seed(11); n <- 400000; g <- 3; p1 <- 0.45; p0 <- 0.15
x <- rnorm(n); a <- rbinom(n, 1, plogis(0.6 * x))
u <- rbinom(n, 1, ifelse(a == 1, p1, p0))          # prevalence depends on A only
lp <- function(aa) -1.2 + 0.25 * aa + 0.8 * x + log(g) * u
y <- rbinom(n, 1, plogis(lp(a)))
truth <- mean(plogis(lp(1)) - plogis(lp(0)))
fs <- glm(y ~ a + x, family = binomial, data = data.frame(y, a, x))
s <- sens_risk_difference(fs, "a", p1 = p1, p0 = p0, lambda = c(1, g))
cat(sprintf("  model without U %.4f; adjusted at the true OR %.4f; truth %.4f\n",
            s$theta[1], s$theta[2], truth))
check("adjusted estimate within 0.002 of the truth", abs(s$theta[2] - truth) < 0.002)
saveRDS(list(naive = s$theta[1], adjusted = s$theta[2], truth = truth),
        file.path(out_dir, "simulation.rds"))

cat("\n6. Double machine learning: export predictions for 04_doubleml_check.py\n")
k401 <- tryCatch(DoubleML::fetch_401k(return_type = "data.table"), error = function(e) NULL)
if (is.null(k401) && file.exists("sipp1991.dta"))
  k401 <- data.table::as.data.table(readstata13::read.dta13("sipp1991.dta"))
if (!is.null(k401)) {
  lgr::get_logger("mlr3")$set_threshold("warn")
  xc <- c("age", "inc", "educ", "fsize", "marr", "twoearn", "db", "pira", "hown")
  dd <- DoubleML::DoubleMLData$new(k401, y_col = "net_tfa", d_cols = "e401", x_cols = xc)
  res <- list()
  set.seed(2026)
  plr <- DoubleML::DoubleMLPLR$new(dd, ml_l = mlr3::lrn("regr.ranger", num.trees = 200),
           ml_m = mlr3::lrn("regr.ranger", num.trees = 200), n_folds = 5, n_rep = 3)
  plr$fit(store_predictions = TRUE)
  set.seed(2026)
  irm <- DoubleML::DoubleMLIRM$new(dd, ml_g = mlr3::lrn("regr.ranger", num.trees = 200),
           ml_m = mlr3::lrn("classif.ranger", num.trees = 200), n_folds = 5,
           trimming_threshold = 0.01)
  irm$fit(store_predictions = TRUE)
  for (nm in c("plr", "irm")) {
    m <- get(nm); p <- sens_dml(m)
    b <- confoundvis:::.dml_bounds(p$meta$dml_elements, 0.03, 0.03, 1, 0.95)
    res[[nm]] <- c(coef = unname(m$coef), rv = p$meta$rv, rva = p$meta$rva,
                   theta_lower = b$lower$theta, se_lower = b$lower$se,
                   ci_lower = b$lower$ci, ci_upper = b$upper$ci)
    pr <- m$predictions; nr <- m$n_rep
    mats <- lapply(names(pr), function(k) {
      z <- matrix(pr[[k]][, , 1], ncol = nr)
      if (nm == "irm" && k == "ml_m") z <- pmin(pmax(z, 0.01), 0.99)  # as in the R score
      colnames(z) <- paste0(k, "_", seq_len(nr)); z })
    utils::write.csv(do.call(cbind, mats), file.path(out_dir, paste0(nm, "_predictions.csv")),
                     row.names = FALSE)
  }
  utils::write.csv(as.data.frame(k401)[, c("net_tfa", "e401", xc)],
                   file.path(out_dir, "k401_data.csv"), row.names = FALSE)
  utils::write.csv(data.frame(model = names(res), do.call(rbind, res)),
                   file.path(out_dir, "dml_r_results.csv"), row.names = FALSE)
  cat("  Written; now run:  python 04_doubleml_check.py\n")
} else cat("  401(k) data not available: skipped\n")

cat(sprintf("\n%d of %d checks passed\n", sum(ok), length(ok)))
if (!all(ok)) stop("Validation failed: ", paste(names(ok)[!ok], collapse = "; "))
