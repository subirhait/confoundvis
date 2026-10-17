test_that("robustness_points handles levels, se, and unreached crossings", {
  x <- new_confoundsens(
    lambda = rep(seq(0, 1, by = 0.1), 2),
    theta  = c(1 - seq(0, 1, by = 0.1) * 2, 1 - seq(0, 1, by = 0.1) * 0.5),
    se     = rep(0.1, 22),
    level  = rep(c("a", "b"), each = 11)
  )
  rp <- robustness_points(x)
  expect_equal(rp$lambda_null, c(0.5, NA))
  expect_equal(rp$lambda_sig[1], (1 - qnorm(0.975) * 0.1) / 2, tolerance = 1e-10)
  expect_identical(rp$sig_rule, c("normal CI", "normal CI"))
})

test_that("robustness_points returns first lambda when already crossed", {
  x <- new_confoundsens(lambda = c(0, 0.1, 0.2), theta = c(0.1, 0.05, 0),
                        se = rep(1, 3))
  expect_equal(robustness_points(x)$lambda_sig, 0)
})

test_that("robustness_points with no uncertainty gives NA rule", {
  x <- new_confoundsens(lambda = c(0, 1), theta = c(1, -1))
  rp <- robustness_points(x)
  expect_equal(rp$lambda_null, 0.5)
  expect_true(is.na(rp$sig_rule))
})

test_that("sens_report prints a caveat and framework summary", {
  fit <- lm(mpg ~ am + wt + hp, data = mtcars)
  r <- sens_report(sens_path_lm(fit, "am"))
  expect_s3_class(r, "confoundvis_report")
  expect_true(any(grepl("robustness value", r$text)))
  expect_true(any(grepl("do not show whether such a confounder exists", r$text)))
  expect_output(print(r), "Sensitivity report")
  r2 <- sens_report(itcv_lm(fit, "am")$path)
  expect_true(any(grepl("ITCV", r2$text)))
  expect_error(robustness_points(list()), "confoundsens")
})

test_that("sens_report distinguishes marginal significance from none", {
  fit <- lm(mpg ~ am + wt + hp + qsec, data = mtcars)  # am: t = 2.09, p = 0.046
  r <- sens_report(sens_path_lm(fit, "am"))
  expect_true(any(grepl("unadjusted estimate is significant", r$text)))
  expect_false(any(grepl("not significant \\(alpha = 0.05\\) even without", r$text)))
  fit2 <- lm(mpg ~ am + wt + hp + qsec + disp + drat, data = mtcars)
  if (summary(fit2)$coefficients["am", 4] > 0.05) {
    r2 <- sens_report(sens_path_lm(fit2, "am"))
    expect_true(any(grepl("not significant", r2$text)))
  }
})

test_that("love plot caption rounds the threshold", {
  fit <- lm(mpg ~ am + wt + hp + qsec, data = mtcars)
  p <- plot_sensitivity_love(covariate_impacts(fit, "am"))
  expect_match(p$labels$caption, "= 0.0103$")
})
