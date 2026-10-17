test_that("plot_reversal_cone returns ggplot and validates", {
  p <- plot_reversal_cone(theta0 = 0.4, delta = 0.2, grid_n = 40)
  expect_s3_class(p, "ggplot")
  expect_s3_class(plot_reversal_cone(-0.3, 0.5, grid_n = 40,
                                     show_boundary = FALSE,
                                     show_volume = FALSE), "ggplot")
  expect_error(plot_reversal_cone(theta0 = 0), "nonzero")
  expect_error(plot_reversal_cone(delta = -1), "positive")
  expect_error(plot_reversal_cone(q_range = c(1, 0)), "q_range")
  expect_error(plot_reversal_cone(p_range = 1), "p_range")
  expect_error(plot_reversal_cone(grid_n = 5), "grid_n")
})

test_that("plot_cone_comparison draws and returns both panels", {
  pdf(NULL); on.exit(dev.off())
  out <- plot_cone_comparison(grid_n = 30)
  expect_named(out, c("within", "between"))
  expect_s3_class(out$within, "ggplot")
})

test_that("plot_taylor_panels returns three panels; old name is deprecated", {
  pdf(NULL); on.exit(dev.off())
  out <- plot_taylor_panels(delta_max = 0.5, step = 0.05)
  expect_named(out, c("A", "B", "C"))
  expect_warning(plot_figure2_taylor_panels(delta_max = 0.5, step = 0.05),
                 "deprecated")
  expect_error(plot_taylor_panels(local_max_delta = 5), "<=")
})

test_that("plot_robustness_curve draws reference lines from meta", {
  fit <- lm(mpg ~ am + wt + hp, data = mtcars)
  it <- itcv_lm(fit, "am")$path
  p <- plot_robustness_curve(it)
  geoms <- vapply(p$layers, function(l) class(l$geom)[1], "")
  expect_equal(sum(geoms == "GeomHline"), 2L)
  p0 <- plot_robustness_curve(it, reference = FALSE)
  expect_false(any(vapply(p0$layers, function(l) class(l$geom)[1], "") == "GeomHline"))
  x <- new_confoundsens(lambda = 1:3, theta = 3:1)
  expect_false(any(vapply(plot_robustness_curve(x)$layers,
                          function(l) class(l$geom)[1], "") == "GeomHline"))
})

test_that("new_confoundsens meta is validated and back-compatible", {
  expect_error(new_confoundsens(1:2, 2:1, meta = "x"), "list")
  old <- structure(list(lambda = 1:2, theta = 2:1), class = "confoundsens")
  expect_output(print(old), "confoundsens")
})
