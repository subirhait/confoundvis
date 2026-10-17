test_that("sens_table collects paths with their own summaries", {
  fit <- lm(mpg ~ am + wt + hp, data = mtcars)
  g <- seq(0, 0.9, by = 0.001)
  p1 <- sens_path_lm(fit, "wt", lambda = g)
  it <- itcv_lm(fit, "wt", lambda = g)
  tb <- sens_table(pr2 = p1, itcv = it$path)
  expect_s3_class(tb, "confoundvis_table")
  expect_equal(nrow(tb), 2L)
  expect_identical(tb$path, c("pr2", "itcv"))
  expect_equal(tb$summary_value[1], p1$meta$rv)
  expect_equal(tb$summary_value[2], it$itcv)
  expect_equal(tb$lambda_sig[2], it$itcv, tolerance = 1e-3)
  expect_equal(tb$lambda_null, c(robustness_points(p1)$lambda_null,
                                 robustness_points(it$path)$lambda_null))
  expect_output(print(tb), "comparable only within a row")
})

test_that("sens_table accepts a list, expands levels, and names defaults", {
  x <- new_confoundsens(lambda = rep(c(0, 0.5, 1), 2), theta = c(1, 0, -1, 2, 1, 0),
                        level = rep(c("a", "b"), each = 3))
  tb <- sens_table(list(x, x))
  expect_equal(nrow(tb), 4L)
  expect_identical(unique(tb$path), c("path1", "path2"))
  expect_true(all(is.na(tb$summary_value)))
  expect_error(sens_table(), "at least one")
  expect_error(sens_table(list(1)), "confoundsens")
})
