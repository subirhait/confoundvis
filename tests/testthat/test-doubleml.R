# Reference values from Python doubleml 0.11.4 (DoubleMLPLR with the same
# cross-fitted predictions passed as external_predictions, then
# sensitivity_analysis(cf_y = 0.04, cf_d = 0.03)).

dml_fixture <- function(store = TRUE) {
  set.seed(2026)
  n <- 300
  X <- matrix(stats::rnorm(n * 3), n, 3); colnames(X) <- paste0("x", 1:3)
  d <- 0.5 * X[, 1] - 0.3 * X[, 2] + stats::rnorm(n)
  y <- 0.6 * d + X[, 1] + 0.5 * X[, 3] + stats::rnorm(n)
  dt <- data.table::data.table(y = y, d = d, X)
  obj <- DoubleML::DoubleMLData$new(dt, y_col = "y", d_cols = "d")
  folds <- rep(1:3, length.out = n)
  smpls <- list(list(train_ids = lapply(1:3, function(k) which(folds != k)),
                     test_ids  = lapply(1:3, function(k) which(folds == k))))
  plr <- DoubleML::DoubleMLPLR$new(obj, ml_l = mlr3::lrn("regr.lm"),
                                   ml_m = mlr3::lrn("regr.lm"), n_folds = 3,
                                   draw_sample_splitting = FALSE)
  plr$set_sample_splitting(smpls)
  plr$fit(store_predictions = store)
  plr
}

test_that("sens_dml() reproduces Python doubleml", {
  skip_if_not_installed("DoubleML")
  skip_if_not_installed("mlr3")
  skip_if_not_installed("data.table")
  skip_if_not_installed("mlr3learners")
  if (requireNamespace("lgr", quietly = TRUE)) lgr::get_logger("mlr3")$set_threshold("warn")
  plr <- dml_fixture()
  p <- sens_dml(plr)
  expect_s3_class(p, "confoundsens")
  expect_equal(p$meta$framework, "dml")
  expect_equal(p$theta[1], 0.586868079205, tolerance = 1e-10)
  b <- confoundvis:::.dml_bounds(p$meta$dml_elements, 0.04, 0.03, 1, 0.95)
  expect_equal(b$lower$theta, 0.5506453709,   tolerance = 1e-9)
  expect_equal(b$lower$se,    0.062225679285, tolerance = 1e-9)
  expect_equal(b$lower$ci,    0.448293236638, tolerance = 1e-9)
  expect_equal(b$upper$ci,    0.725689051287, tolerance = 1e-9)
  # doubleml finds RV and RVa with minimize_scalar (tolerance about 1e-5)
  expect_equal(p$meta$rv,  0.430167714158, tolerance = 1e-5)
  expect_equal(p$meta$rva, 0.36847958698,  tolerance = 1e-5)
})

test_that("sens_dml() path crossings agree with RV and RVa", {
  skip_if_not_installed("DoubleML")
  skip_if_not_installed("mlr3")
  skip_if_not_installed("data.table")
  skip_if_not_installed("mlr3learners")
  if (requireNamespace("lgr", quietly = TRUE)) lgr::get_logger("mlr3")$set_threshold("warn")
  plr <- dml_fixture()
  p <- sens_dml(plr)
  rp <- robustness_points(p)
  expect_equal(rp$alpha, 0.10)
  expect_equal(rp$lambda_null, p$meta$rv, tolerance = 1e-4)
  expect_equal(rp$lambda_sig,  p$meta$rva, tolerance = 1e-4)
  # two-sided 0.05 equals one-sided level 0.975
  p975 <- sens_dml(plr, level = 0.975)
  expect_equal(robustness_points(p, alpha = 0.05)$lambda_sig, p975$meta$rva,
               tolerance = 1e-4)
  # rho scales the bound; a null above the estimate uses the upper bound
  pu <- sens_dml(plr, null = 1)
  expect_equal(pu$meta$side, "upper")
  expect_true(pu$theta[2] > pu$theta[1])
  expect_true(sens_dml(plr, rho = 0.5)$meta$rv > p$meta$rv)
})

test_that("sens_dml() surface, report, and confoundvis()", {
  skip_if_not_installed("DoubleML")
  skip_if_not_installed("mlr3")
  skip_if_not_installed("data.table")
  skip_if_not_installed("mlr3learners")
  if (requireNamespace("lgr", quietly = TRUE)) lgr::get_logger("mlr3")$set_threshold("warn")
  plr <- dml_fixture()
  p <- sens_dml(plr)
  s <- sens_surface(p, n = 41)
  expect_s3_class(s, "confoundsurf")
  # every point of the null boundary puts the bound at the null
  nb <- s$boundary[grepl("^Null", s$boundary$kind), ]
  v <- confoundvis:::.dml_bounds(p$meta$dml_elements, nb$b, nb$a, 1, 0.95)$lower$theta
  expect_lt(max(abs(v)), 1e-8)
  sb <- s$boundary[grepl("^Significance", s$boundary$kind), ]
  v <- confoundvis:::.dml_bounds(p$meta$dml_elements, sb$b, sb$a, 1, 0.95)$lower$ci
  expect_lt(max(abs(v)), 1e-6)
  expect_s3_class(plot_sensitivity_surface(s), "ggplot")
  expect_s3_class(plot_robustness_curve(p), "ggplot")
  expect_match(paste(sens_report(p)$text, collapse = " "), "RVa")
  r <- confoundvis(plr)
  expect_equal(r$type, "dml")
  expect_equal(r$table$alpha, 0.05)
  expect_output(print(r), "DoubleML")
})

test_that("sens_dml() input checks", {
  expect_error(sens_dml(lm(mpg ~ wt, mtcars)), "DoubleML")
  skip_if_not_installed("DoubleML")
  skip_if_not_installed("mlr3")
  skip_if_not_installed("data.table")
  skip_if_not_installed("mlr3learners")
  if (requireNamespace("lgr", quietly = TRUE)) lgr::get_logger("mlr3")$set_threshold("warn")
  plr <- dml_fixture()
  expect_error(sens_dml(plr, rho = 2), "rho")
  expect_error(sens_dml(plr, level = 1), "level")
  expect_error(sens_dml(dml_fixture(store = FALSE)), "store_predictions")
})
