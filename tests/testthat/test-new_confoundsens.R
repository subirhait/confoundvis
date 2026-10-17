test_that("new_confoundsens creates valid object", {
  x <- new_confoundsens(lambda = seq(0, 0.2, length.out = 5),
                        theta  = seq(1, 0.8, length.out = 5))
  expect_s3_class(x, "confoundsens")
  expect_length(x$lambda, 5)
  expect_null(x$level)
  expect_null(x$se)
  expect_null(x$t)
})

test_that("new_confoundsens stores optional fields", {
  x <- new_confoundsens(
    lambda = 1:4,
    theta  = c(0.5, 0.4, 0.3, 0.2),
    se     = rep(0.05, 4),
    t      = c(3.1, 2.8, 2.4, 1.9),
    level  = c("female", "female", "male", "male")
  )
  expect_equal(x$se, rep(0.05, 4))
  expect_equal(x$t,  c(3.1, 2.8, 2.4, 1.9))
  expect_equal(x$level, c("female", "female", "male", "male"))
})

test_that("new_confoundsens coerces level to character", {
  x <- new_confoundsens(lambda = 1:3, theta = 3:1,
                        level = factor(c("a", "b", "a")))
  expect_type(x$level, "character")
})

test_that("new_confoundsens errors on non-numeric lambda", {
  expect_error(new_confoundsens(lambda = "a", theta = 1),
               "`lambda` must be a numeric vector")
})

test_that("new_confoundsens errors on length mismatch", {
  expect_error(new_confoundsens(lambda = 1:5, theta = 1:4),
               "same length")
})

test_that("new_confoundsens warns on non-monotone lambda", {
  expect_warning(
    new_confoundsens(lambda = c(0.1, 0.0, 0.2), theta = 1:3),
    "nondecreasing"
  )
})
