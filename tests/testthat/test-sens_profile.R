test_that("sens_profile collects the three regression analyses", {
  fit <- lm(mpg ~ am + wt + hp, data = mtcars)
  pr <- sens_profile(fit, "wt")
  expect_s3_class(pr, "confoundvis_profile")
  expect_named(pr$paths, c("partial_r2", "itcv", "rir"))
  expect_equal(pr$table$summary_value[1], pr$paths$partial_r2$meta$rv)
  expect_equal(pr$table$lambda_sig[2], pr$itcv$itcv, tolerance = 1e-4)
  expect_equal(pr$table$lambda_sig[3], pr$rir$perc_bias / 100, tolerance = 1e-4)
  expect_output(print(pr), "confoundvis_profile")
  expect_s3_class(plot(pr), "ggplot")
})

test_that("sens_profile accepts extra paths and checks inputs", {
  fit <- lm(mpg ~ am + wt + hp, data = mtcars)
  ex <- new_confoundsens(lambda = seq(0, 1, 0.1), theta = 1 - seq(0, 1, 0.1),
                         meta = list(null = 0))
  pr <- sens_profile(fit, "wt", extra = list(user = ex))
  expect_equal(nrow(pr$table), 4)
  expect_error(sens_profile(fit, "wt", extra = list(ex)), "named")
  expect_error(sens_profile(fit, "wt", alpha = 3), "alpha")
})
