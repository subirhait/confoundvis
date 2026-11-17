# Paths with a stored threshold fix alpha when they are built (reviewer report:
# robustness_points() previously ignored a different alpha silently).

test_that("threshold paths use their stored alpha and refuse a different one", {
  fit <- lm(mpg ~ am + wt + hp, data = mtcars)
  g <- seq(0, 0.9, by = 0.001)
  it05 <- itcv_lm(fit, "wt", alpha = 0.05, lambda = g)
  it01 <- itcv_lm(fit, "wt", alpha = 0.01, lambda = g)
  expect_equal(robustness_points(it05$path)$alpha, 0.05)
  expect_equal(robustness_points(it05$path, alpha = 0.05)$alpha, 0.05)
  expect_error(robustness_points(it05$path, alpha = 0.01), "Rebuild it at alpha = 0.01")
  expect_error(sens_report(it05$path, alpha = 0.01), "built with alpha = 0.05")
  p01 <- robustness_points(it01$path)
  expect_equal(p01$alpha, 0.01)
  expect_equal(p01$lambda_sig, it01$itcv, tolerance = 1e-3)
  expect_false(isTRUE(all.equal(p01$lambda_sig,
                                robustness_points(it05$path)$lambda_sig)))
  expect_match(paste(sens_report(it01$path)$text, collapse = " "), "alpha = 0.01")

  rr <- rir_lm(fit, "wt", alpha = 0.01)
  expect_equal(robustness_points(rr$path)$alpha, 0.01)
  expect_error(robustness_points(rr$path, alpha = 0.05), "rir_lm")

  d <- c(1.2, 0.8, 2.1, 0.4, 1.7, -0.3, 0.9, 1.5, 2.2, 0.6, 1.1, 0.7)
  pp <- sens_path_pairs(d, gamma = seq(1, 4, by = 0.01), alpha = 0.025)
  expect_equal(robustness_points(pp)$alpha, 0.025)
  expect_error(robustness_points(pp, alpha = 0.05), "sens_path_pairs")
})

test_that("paths without a stored threshold honour alpha", {
  fit <- lm(mpg ~ am + wt + hp, data = mtcars)
  p <- sens_path_lm(fit, "wt", lambda = seq(0, 0.9, by = 0.001))
  a <- robustness_points(p)
  b <- robustness_points(p, alpha = 0.01)
  expect_equal(a$alpha, 0.05)
  expect_equal(b$alpha, 0.01)
  expect_lt(b$lambda_sig, a$lambda_sig)
  expect_error(robustness_points(p, alpha = 2), "single number")
})

test_that("substantive RIR threshold carries no alpha", {
  fit <- lm(mpg ~ am + wt + hp, data = mtcars)
  rr <- rir_lm(fit, "wt", threshold = 2)
  expect_true(is.na(robustness_points(rr$path)$alpha))
  expect_warning(robustness_points(rr$path, alpha = 0.01), "ignored")
})

test_that("sens_table reports the level used in each row", {
  fit <- lm(mpg ~ am + wt + hp, data = mtcars)
  g <- seq(0, 0.9, by = 0.001)
  tb <- sens_table(a = sens_path_lm(fit, "wt", lambda = g),
                   b = itcv_lm(fit, "wt", alpha = 0.01, lambda = g)$path)
  expect_equal(tb$alpha, c(0.05, 0.01))
})

test_that("glm fits point to konfound and from_konfound()", {
  set.seed(2)
  n <- 400
  x1 <- rnorm(n)
  d <- rbinom(n, 1, plogis(0.3 * x1))
  y <- rbinom(n, 1, plogis(-0.5 + 0.9 * d + 0.4 * x1))
  gfit <- glm(y ~ d + x1, family = binomial)
  expect_error(rir_lm(gfit, "d"), "from_konfound")
  expect_error(itcv_lm(gfit, "d"), "konfound")
  skip_if_not_installed("konfound")
  k <- suppressWarnings(suppressMessages(
    konfound::konfound(gfit, d, index = "RIR", to_return = "raw_output")))
  p <- from_konfound(k, index = "RIR", treatment = "d")
  expect_s3_class(p, "confoundsens")
  expect_equal(p$meta$rir, k$RIR_primary)
  expect_equal(robustness_points(p)$lambda_sig, k$perc_bias_to_change / 100,
               tolerance = 1e-3)
})

test_that("covariate_impacts labels and factor handling", {
  m <- mtcars
  m$cylf <- factor(m$cyl)
  fit <- lm(mpg ~ am + cylf + hp, data = m)
  ci <- covariate_impacts(fit, "am", metric = "partial_r2",
                          labels = c(hp = "Gross horsepower"))
  expect_equal(ci$label, c("cylf", "Gross horsepower"))
  expect_equal(ci$term_df, c(2L, 1L))
  expect_true(is.na(ci$r_du[1]))
  expect_false(is.na(ci$bias_index[1]))
  p <- plot_sensitivity_love(ci)
  expect_true("Gross horsepower" %in% levels(p$data$covariate))
  m$cylo <- ordered(m$cyl)
  co <- covariate_impacts(lm(mpg ~ am + cylo + hp, data = m), "am",
                          metric = "partial_r2")
  expect_equal(co$bias_index, ci$bias_index)
  expect_error(covariate_impacts(fit, "am", labels = c("x")), "named")
})
