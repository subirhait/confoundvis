test_that("rir_lm matches konfound's RIR and percent bias", {
  skip_if_not_installed("konfound")
  fit <- lm(mpg ~ am + wt + hp, data = mtcars)
  rr <- rir_lm(fit, "wt")
  k <- suppressWarnings(suppressMessages(
    konfound::konfound(fit, wt, index = "RIR", to_return = "raw_output")))
  expect_equal(rr$beta_threshold, k$beta_threshold, tolerance = 1e-10)
  expect_equal(rr$perc_bias, k$perc_bias_to_change, tolerance = 1e-8)
  expect_equal(rr$rir, k$RIR_primary)
  p <- from_konfound(k, index = "RIR")
  expect_equal(p$meta$rir, k$RIR_primary)
  expect_equal(p$meta$estimate, unname(coef(fit)["wt"]), tolerance = 1e-8)
  expect_equal(robustness_points(p)$lambda_sig, rr$perc_bias / 100, tolerance = 1e-3)
})

test_that("rir path crosses the threshold at the percent bias", {
  fit <- lm(mpg ~ am + wt + hp, data = mtcars)
  rr <- rir_lm(fit, "hp", lambda = seq(0, 1, by = 1e-4))
  est <- unname(coef(fit)["hp"]); se <- coef(summary(fit))["hp", 2]
  bt <- sign(est) * qt(0.975, df.residual(fit)) * se
  expect_equal(rr$perc_bias, 100 * (1 - bt / est))
  expect_equal(rr$rir, round(32 * (1 - bt / est)))
  rp <- robustness_points(rr$path)
  expect_equal(rp$lambda_sig, 1 - bt / est, tolerance = 1e-6)
  expect_equal(rp$lambda_null, 1)
  r <- sens_report(rr$path)
  expect_true(any(grepl("RIR =", r$text)))
  expect_false(any(grepl("covariate_impacts", r$text)))
  expect_equal(sens_table(rir = rr$path)$summary_value, rr$rir)
})

test_that("rir_lm handles a non-significant estimate and bad input", {
  fit <- lm(mpg ~ am + wt + hp + qsec + drat, data = mtcars)
  rr <- rir_lm(fit, "drat")
  expect_true(is.na(rr$rir) && is.na(rr$perc_bias))
  expect_equal(robustness_points(rr$path)$lambda_sig, 0)
  expect_true(any(grepl("does not exceed", sens_report(rr$path)$text)))
  expect_error(rir_lm(fit, "wt", lambda = c(0, 1.5)), "\\[0, 1\\]")
  expect_error(from_konfound(list(a = 1), index = "RIR"), "beta_threshold")
})
