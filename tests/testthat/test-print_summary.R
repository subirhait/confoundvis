test_that("print.confoundsens runs without error", {
  x <- new_confoundsens(lambda = 1:5, theta = 5:1)
  expect_output(print(x), "<confoundsens>")
  expect_output(print(x), "n      : 5")
})

test_that("summary.confoundsens produces correct structure", {
  x   <- new_confoundsens(lambda = seq(0, 0.4, length.out = 10),
                          theta  = seq(1, 0.6, length.out = 10),
                          level  = rep(c("female","male"), 5))
  s   <- summary(x)
  expect_s3_class(s, "summary.confoundsens")
  expect_s3_class(s$table, "data.frame")
  expect_equal(nrow(s$table), 2L)   # two levels
  expect_true(all(c("level","n","theta_0") %in% names(s$table)))
})

test_that("summary.confoundsens pooled when no level", {
  x <- new_confoundsens(lambda = 0:4, theta = 4:0)
  s <- summary(x)
  expect_equal(s$table$level, "pooled")
})

test_that("print.summary.confoundsens produces output", {
  x <- new_confoundsens(lambda = 1:3, theta = 3:1)
  expect_output(print(summary(x)), "pooled")
})
