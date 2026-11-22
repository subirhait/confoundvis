iv_fixture <- function(alpha = 0.05) {
  data("card", package = "iv.sensemakr", envir = environment())
  xm <- stats::model.matrix(~ exper + expersq + black + south + smsa + reg661 +
    reg662 + reg663 + reg664 + reg665 + reg666 + reg667 + reg668 + smsa66,
    data = card)
  fit <- iv.sensemakr::iv_fit(card$lwage, card$educ, card$nearc4, xm)
  iv.sensemakr::sensemakr(fit, benchmark_covariates = c("black", "smsa"),
                          alpha = alpha)
}

test_that("from_iv_sensemakr() agrees with iv.sensemakr", {
  skip_if_not_installed("iv.sensemakr")
  s <- iv_fixture()
  p <- from_iv_sensemakr(s)
  expect_equal(p$meta$framework, "iv")
  expect_equal(p$theta[1], s$unadjusted$estimates$iv$estimate)
  # robustness values reported by iv.sensemakr
  expect_equal(p$meta$rv, s$sensitivity_stats$iv$rv_qa)
  expect_equal(p$meta$xrv, s$sensitivity_stats$iv$xrv_qa)
  expect_equal(p$meta$rv_point,
               as.numeric(iv.sensemakr::rv(s$unadjusted, alpha = 1)), tolerance = 1e-6)
  rp <- robustness_points(p)
  expect_equal(rp$lambda_sig, p$meta$rv, tolerance = 1e-5)
  expect_equal(rp$lambda_null, p$meta$rv_point, tolerance = 1e-4)
  expect_equal(rp$sig_rule, "adjusted confidence limit")
  # at lambda = 0 the adjusted interval is the Anderson-Rubin interval, up to
  # the df - 1 convention of sensemakr::adjusted_critical_value()
  expect_equal(c(p$meta$ci_lower[1], p$meta$ci_upper[1]),
               s$unadjusted$estimates$iv$conf.int, tolerance = 1e-3)
  # benchmark bounds
  L <- confoundvis:::.iv_limits(p$meta$iv_parts, s$bounds$iv$r2zw.x,
                                s$bounds$iv$r2y0w.zx, 0.05)
  expect_equal(L$lwr, s$bounds$iv$lwr, tolerance = 1e-10)
  expect_equal(L$upr, s$bounds$iv$upr, tolerance = 1e-10)
})

test_that("from_iv_sensemakr() keeps the level of the iv.sensemakr object", {
  skip_if_not_installed("iv.sensemakr")
  s <- iv_fixture(alpha = 0.10)
  p <- from_iv_sensemakr(s)
  expect_equal(robustness_points(p)$alpha, 0.10)
  expect_equal(robustness_points(p)$lambda_sig, s$sensitivity_stats$iv$rv_qa,
               tolerance = 1e-5)
  expect_error(robustness_points(p, alpha = 0.05), "iv.sensemakr::sensemakr")
  expect_error(from_iv_sensemakr(lm(mpg ~ wt, mtcars)), "iv.sensemakr")
})

test_that("IV surface, plots, report, and confoundvis()", {
  skip_if_not_installed("iv.sensemakr")
  s <- iv_fixture()
  p <- from_iv_sensemakr(s)
  sf <- sens_surface(p, n = 41)
  expect_s3_class(sf, "confoundsurf")
  # every point of the significance boundary puts the adjusted limit at h0
  sb <- sf$boundary[grepl("^Significance", sf$boundary$kind), ]
  v <- confoundvis:::.iv_near(p$meta$iv_parts, sb$a, sb$b, 0.05, 1)
  expect_lt(max(abs(v)), 1e-6)
  expect_s3_class(plot_sensitivity_surface(sf), "ggplot")
  expect_s3_class(plot_robustness_curve(p), "ggplot")
  expect_match(paste(sens_report(p)$text, collapse = " "), "XRV")
  r <- confoundvis(s)
  expect_equal(r$type, "iv")
  expect_output(print(r), "Anderson-Rubin")
  r2 <- confoundvis(s, alpha = 0.10)
  expect_equal(r2$table$alpha, 0.10)
  expect_equal(r2$table$lambda_sig, iv_fixture(0.10)$sensitivity_stats$iv$rv_qa,
               tolerance = 1e-5)
})
