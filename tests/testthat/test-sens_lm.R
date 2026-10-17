fit <- lm(mpg ~ am + wt + hp + qsec, data = mtcars)

test_that("sens_path_lm returns a confoundsens path starting at the estimate", {
  p <- sens_path_lm(fit, "am")
  expect_s3_class(p, "confoundsens")
  expect_equal(p$theta[1], unname(coef(fit)["am"]))
  expect_identical(p$meta$framework, "partial_r2")
  expect_true(all(diff(sign(p$theta[1]) * p$theta) <= 1e-12))  # moves toward zero
  expect_length(p$se, length(p$lambda))
})

test_that("equal-strength path reaches zero at the robustness value", {
  p <- sens_path_lm(fit, "am", lambda = seq(0, 0.9, by = 1e-4))
  rp <- robustness_points(p)
  expect_equal(rp$lambda_null, p$meta$rv, tolerance = 1e-3)
})

test_that("sens_path_lm validates inputs", {
  expect_error(sens_path_lm(fit, "nope"), "not a column")
  expect_error(sens_path_lm(fit, "am", lambda = c(0, 1.2)), "\\[0, 1\\)")
  expect_error(sens_path_lm(fit, "am", ratio = -1), "positive")
  g <- glm(am ~ wt, data = mtcars, family = binomial)
  expect_error(sens_path_lm(g, "wt"), "stats::lm")
})

test_that("sens_path_lm matches sensemakr adjusted estimates and RV", {
  skip_if_not_installed("sensemakr")
  lam <- c(0.01, 0.05, 0.1, 0.2)
  p <- sens_path_lm(fit, "am", lambda = lam)
  ref_est <- sensemakr::adjusted_estimate(fit, treatment = "am",
                                          r2dz.x = lam, r2yz.dx = lam)
  ref_se  <- sensemakr::adjusted_se(fit, treatment = "am",
                                    r2dz.x = lam, r2yz.dx = lam)
  expect_equal(p$theta, unname(ref_est), tolerance = 1e-8)
  expect_equal(p$se, unname(ref_se), tolerance = 1e-8)
  rv <- sensemakr::robustness_value(fit, covariates = "am", q = 1, alpha = 1)
  expect_equal(p$meta$rv, unname(as.numeric(rv)), tolerance = 1e-8)
})

test_that("ratio != 1 matches sensemakr", {
  skip_if_not_installed("sensemakr")
  lam <- c(0.02, 0.08)
  p <- sens_path_lm(fit, "am", lambda = lam, ratio = 2)
  ref <- sensemakr::adjusted_estimate(fit, treatment = "am",
                                      r2dz.x = lam, r2yz.dx = 2 * lam)
  expect_equal(p$theta, unname(ref), tolerance = 1e-8)
})

test_that("itcv_lm follows Frank (2000) and its path hits r_crit at ITCV", {
  it <- itcv_lm(fit, "am", lambda = seq(0, 0.9, by = 1e-4))
  s <- summary(fit)$coefficients["am", ]
  dof <- fit$df.residual
  r <- s[["t value"]] / sqrt(s[["t value"]]^2 + dof)
  tc <- qt(0.975, dof); rc <- tc / sqrt(tc^2 + dof)
  expect_equal(it$r, r)
  expect_equal(it$itcv, (abs(r) - rc) / (1 - rc))
  rp <- robustness_points(it$path)
  expect_equal(rp$lambda_sig, it$itcv, tolerance = 1e-3)
  expect_identical(rp$sig_rule, "theta threshold")
})

test_that("itcv_lm agrees with konfound", {
  skip_if_not_installed("konfound")
  k <- suppressWarnings(suppressMessages(
    konfound::konfound(fit, am, to_return = "raw_output")))
  it <- itcv_lm(fit, "am")
  expect_equal(it$r, k$obs_r, tolerance = 1e-8)
  expect_equal(it$r_crit, k$critical_r, tolerance = 1e-8)
  expect_equal(it$itcv, k$itcvGz, tolerance = 1e-8)             # conditional
  expect_equal(it$itcv_unconditional, k$itcv, tolerance = 1e-8) # unconditional
})

test_that("covariate_impacts reproduces partial correlations and partial R2", {
  imp <- covariate_impacts(fit, "am")
  expect_s3_class(imp, "data.frame")
  expect_setequal(imp$covariate, c("wt", "hp", "qsec"))
  # partial correlation of wt with am given hp, qsec
  e1 <- resid(lm(am ~ hp + qsec, mtcars)); e2 <- resid(lm(wt ~ hp + qsec, mtcars))
  expect_equal(abs(imp$r_du[imp$covariate == "wt"]), abs(cor(e1, e2)),
               tolerance = 1e-10)
  # partial R2 of wt in the outcome model given am, hp, qsec
  f0 <- lm(mpg ~ am + hp + qsec, mtcars)
  r2 <- 1 - sum(resid(fit)^2) / sum(resid(f0)^2)
  expect_equal(imp$r2yz.dx[imp$covariate == "wt"], r2, tolerance = 1e-10)
  expect_equal(attr(imp, "threshold"), itcv_lm(fit, "am")$itcv)
  imp2 <- covariate_impacts(fit, "am", metric = "partial_r2")
  expect_equal(imp2$impact, imp2$bias_index)
  expect_equal(attr(imp2, "threshold"),
               abs(summary(fit)$coefficients["am", 3]) / sqrt(fit$df.residual))
})

test_that("covariate_impacts partial R2 matches sensemakr partial_r2", {
  skip_if_not_installed("sensemakr")
  imp <- covariate_impacts(fit, "am", metric = "partial_r2")
  ref <- sensemakr::partial_r2(fit, covariates = "wt")
  expect_equal(imp$r2yz.dx[imp$covariate == "wt"], unname(ref),
               tolerance = 1e-10)
})

test_that("covariate_impacts handles factor terms", {
  d <- mtcars; d$cyl <- factor(d$cyl)
  f <- lm(mpg ~ am + wt + cyl, data = d)
  expect_message(imp <- covariate_impacts(f, "am"), "Multi-column")
  expect_true(is.na(imp$impact[imp$covariate == "cyl"]))
  expect_equal(imp$term_df[imp$covariate == "cyl"], 2L)
  imp2 <- covariate_impacts(f, "am", metric = "partial_r2")
  expect_false(anyNA(imp2$impact))
})

test_that("covariate_impacts feeds love and contour plots", {
  imp <- covariate_impacts(fit, "am")
  expect_s3_class(plot_sensitivity_love(imp), "ggplot")
  expect_s3_class(
    plot_sensitivity_contour(attr(imp, "threshold"), benchmarks = imp), "ggplot")
})

test_that("weighted lm is supported", {
  fw <- lm(mpg ~ am + wt, data = mtcars, weights = hp)
  p <- sens_path_lm(fw, "am")
  expect_equal(p$theta[1], unname(coef(fw)["am"]))
  skip_if_not_installed("sensemakr")
  ref <- sensemakr::adjusted_estimate(fw, treatment = "am",
                                      r2dz.x = 0.1, r2yz.dx = 0.1)
  expect_equal(sens_path_lm(fw, "am", lambda = c(0, 0.1))$theta[2],
               unname(ref), tolerance = 1e-8)
})

test_that("unconditional ITCV uses weighted R2 for weighted fits", {
  set.seed(2)
  n <- 200
  x1 <- rnorm(n); x2 <- rnorm(n); D <- 0.5 * x1 + rnorm(n)
  y <- D + x1 + x2 + rnorm(n); w <- runif(n, 0.2, 3)
  it <- itcv_lm(lm(y ~ D + x1 + x2, weights = w), "D")
  expect_equal(it$r2yz, summary(lm(y ~ x1 + x2, weights = w))$r.squared,
               tolerance = 1e-10)
  expect_equal(it$r2xz, summary(lm(D ~ x1 + x2, weights = w))$r.squared,
               tolerance = 1e-10)
})

test_that("ITCV benchmarks follow the practical guide (eqs. 10 and 12)", {
  fit <- lm(mpg ~ am + wt + hp + qsec, data = mtcars)
  it <- itcv_lm(fit, "wt")
  expect_equal(it$benchmark_corr_product, sqrt(it$r2xz * it$r2yz))
  expect_equal(it$itcv_ratio_to_benchmark,
               it$itcv_unconditional / sqrt(it$r2xz * it$r2yz))
  ci <- covariate_impacts(fit, "wt")
  hp <- ci[ci$covariate == "hp", ]
  expect_equal(hp$r_du_zero, cor(mtcars$hp, mtcars$wt), tolerance = 1e-10)
  expect_equal(hp$r_yu_zero, cor(mtcars$hp, mtcars$mpg), tolerance = 1e-10)
  expect_equal(hp$impact_zero, sign(coef(fit)[["wt"]]) * hp$r_du_zero * hp$r_yu_zero)
  ok <- !is.na(ci$impact_zero) & ci$impact_zero > 0
  expect_equal(ci$itcv_ratio[ok], it$itcv_unconditional / ci$impact_zero[ok])
  expect_true(all(is.na(ci$itcv_ratio[!ok])))
})
