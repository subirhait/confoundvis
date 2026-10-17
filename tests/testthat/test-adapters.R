test_that("from_sensemakr reproduces RV and RV_alpha", {
  skip_if_not_installed("sensemakr")
  data("darfur", package = "sensemakr")
  fit <- lm(peacefactor ~ directlyharmed + age + farmer_dar + herder_dar +
              pastvoted + hhsize_darfur + female + village, data = darfur)
  s <- sensemakr::sensemakr(fit, treatment = "directlyharmed",
                            benchmark_covariates = "female", kd = 1:3)
  p <- from_sensemakr(s, lambda = seq(0, 0.5, by = 1e-4))
  expect_s3_class(p, "confoundsens")
  expect_identical(p$meta$source, "sensemakr")
  expect_s3_class(p$meta$bounds, "data.frame")
  rp <- robustness_points(p)
  expect_equal(rp$lambda_null, as.numeric(s$sensitivity_stats$rv_q), tolerance = 1e-3)
  expect_equal(rp$lambda_sig,  as.numeric(s$sensitivity_stats$rv_qa), tolerance = 1e-3)
  expect_equal(as_confoundsens(s)$theta, from_sensemakr(s)$theta)
})

test_that("from_sensemakr rejects other objects", {
  expect_error(from_sensemakr(list()), "sensemakr")
})

test_that("from_konfound reproduces ITCV from konfound output", {
  skip_if_not_installed("konfound")
  k <- suppressMessages(konfound::pkonfound(
    est_eff = 2, std_err = 0.4, n_obs = 100, n_covariates = 3,
    index = "IT", to_return = "raw_output"))
  p <- from_konfound(k, lambda = seq(0, 0.9, by = 1e-4))
  rc <- abs(k$critical_r); r <- abs(k$obs_r)
  expect_equal(p$meta$itcv, (r - rc) / (1 - rc))
  expect_equal(robustness_points(p)$lambda_sig, p$meta$itcv, tolerance = 1e-3)
  expect_null(p$t)
  expect_length(from_konfound(k, dof = 95)$t, 51L)
})

test_that("from_konfound validates input", {
  expect_error(from_konfound(list(a = 1)), "obs_r")
})

test_that("from_evalue path crosses the null at the E-value", {
  skip_if_not_installed("EValue")
  e <- EValue::evalues.RR(est = 1.8, lo = 1.3, hi = 2.5)
  p <- from_evalue(e, lambda = seq(1, 5, by = 1e-4))
  expect_equal(robustness_points(p)$lambda_null,
               unname(e["E-values", "point"]), tolerance = 1e-3)
  # protective effect
  e2 <- EValue::evalues.RR(est = 0.6, lo = 0.4, hi = 0.9)
  p2 <- from_evalue(e2, lambda = seq(1, 5, by = 1e-4))
  expect_equal(robustness_points(p2)$lambda_null,
               unname(e2["E-values", "point"]), tolerance = 1e-3)
})

test_that("from_evalue validates input", {
  expect_error(from_evalue(matrix(1)), "evalue")
})

test_that("as_confoundsens dispatches on lm and rejects unknown classes", {
  fit <- lm(mpg ~ am + wt, data = mtcars)
  expect_equal(as_confoundsens(fit, treatment = "am")$theta,
               sens_path_lm(fit, "am")$theta)
  expect_error(as_confoundsens(1:3), "must be a data.frame")
})
