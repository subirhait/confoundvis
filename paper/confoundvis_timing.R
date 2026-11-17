# Timing of confoundvis 0.3.1 on simulated data (reported in the article).
# Each step is timed with system.time(); the median of `reps` runs is kept.
library(confoundvis); library(ggplot2)
set.seed(20261117)
p <- 20; grid <- seq(0, 0.12, by = 1e-4)          # 1,201 grid points
timeit <- function(expr, reps) {
  e <- substitute(expr); f <- parent.frame()
  median(replicate(reps, system.time(eval(e, f))[["elapsed"]]))
}
out <- NULL
for (n in c(1e4, 1e5, 1e6)) {
  reps <- if (n < 1e6) 3 else 1
  X <- matrix(rnorm(n * p), n, p, dimnames = list(NULL, paste0("x", 1:p)))
  d <- rbinom(n, 1, plogis(X[, 1:5] %*% rep(0.3, 5)))
  y <- 0.1 * d + drop(X[, 1:10] %*% rep(0.2, 10)) + rnorm(n)
  dat <- data.frame(y, d, X)
  t_lm <- timeit(fit <- lm(y ~ ., data = dat), reps)
  t_r2 <- timeit(pr <- sens_path_lm(fit, "d", lambda = grid), reps)
  t_it <- timeit(it <- itcv_lm(fit, "d", lambda = grid), reps)
  t_rr <- timeit(rr <- rir_lm(fit, "d"), reps)
  t_ci <- timeit(ci <- covariate_impacts(fit, "d"), reps)
  t_sf <- timeit(sf <- sens_surface(pr, n = 201), reps)
  t_pf <- timeit(prof <- sens_profile(fit, "d"), reps)
  tf <- tempfile(fileext = ".pdf")
  t_p1 <- timeit(ggsave(tf, plot_robustness_curve(pr, points = FALSE), width = 5, height = 4), reps)
  t_p2 <- timeit(ggsave(tf, plot_sensitivity_surface(sf), width = 5, height = 4), reps)
  t_p3 <- timeit(ggsave(tf, plot_sensitivity_love(ci), width = 5, height = 4), reps)
  out <- rbind(out, data.frame(n = n, covariates = p, lm_fit = t_lm,
    sens_path_lm = t_r2, itcv_lm = t_it, rir_lm = t_rr,
    covariate_impacts = t_ci, sens_surface = t_sf, sens_profile = t_pf,
    plot_curve = t_p1, plot_surface = t_p2, plot_love = t_p3))
  print(out[nrow(out), ])
}
write.csv(out, "timing.csv", row.names = FALSE)
cat(R.version.string, "\n"); cat(Sys.info()[["sysname"]], "\n")
