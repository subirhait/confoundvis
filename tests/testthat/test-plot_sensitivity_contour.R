test_that("plot_sensitivity_contour returns a ggplot", {
  p <- plot_sensitivity_contour(threshold = 0.02)
  expect_s3_class(p, "ggplot")
})

test_that("plot_sensitivity_contour with benchmarks returns ggplot", {
  b <- data.frame(r_yu = c(0.1, 0.2), r_du = c(0.15, 0.08),
                  label = c("SES", "BMI"))
  p <- plot_sensitivity_contour(threshold = 0.02, benchmarks = b)
  expect_s3_class(p, "ggplot")
})

test_that("plot_sensitivity_contour errors on non-positive threshold", {
  expect_error(plot_sensitivity_contour(threshold = -0.01), "positive")
  expect_error(plot_sensitivity_contour(threshold = 0),     "positive")
})

test_that("plot_sensitivity_contour errors on bad benchmarks", {
  expect_error(
    plot_sensitivity_contour(threshold = 0.02,
                             benchmarks = data.frame(x = 1)),
    "columns `r_yu` and `r_du`"
  )
})

test_that("plot_sensitivity_contour boundary has four groups", {
  p    <- plot_sensitivity_contour(threshold = 0.05)
  bdat <- p$data
  expect_true("group" %in% names(bdat))
  expect_equal(length(unique(bdat$group)), 4L)
})

test_that("exact boundary passes through the equal-magnitude ITCV point", {
  fit <- lm(mpg ~ am + wt + hp, data = mtcars)
  it  <- itcv_lm(fit, "wt")
  p   <- plot_sensitivity_contour(itcv = it)
  expect_s3_class(p, "ggplot")
  r <- it$r; a <- sign(r) * sqrt(it$itcv)   # adverse equal-magnitude pair
  radj <- (r - a * (-a)) / (1 - a^2)        # opposite signs reduce a negative r
  expect_equal(abs(radj), it$r_crit, tolerance = 1e-10)
  layers <- vapply(p$layers, function(l) class(l$geom)[1], character(1))
  expect_true(all(c("GeomRaster", "GeomPath") %in% layers))
  expect_error(plot_sensitivity_contour(itcv = list(r = 2, r_crit = 0.1)), "in \\(-1, 1\\)")
  expect_error(plot_sensitivity_contour(0.02, limits = c(1, -1)), "limits")
})

test_that("closed-form exact boundaries have |r_adj| equal to the critical r", {
  for (rc in list(c(0.4, 0.2), c(-0.5, 0.3))) {
    ex <- confoundvis:::.exact_boundary(rc[1], rc[2], 2000)
    radj <- (rc[1] - ex$r_yu * ex$r_du) /
      sqrt((1 - ex$r_yu^2) * (1 - ex$r_du^2))
    o <- ex$kind == "Original direction"
    expect_true(any(o) && any(!o))
    expect_equal(sign(rc[1]) * radj[o], rep(rc[2], sum(o)), tolerance = 1e-8)
    expect_equal(sign(rc[1]) * radj[!o], rep(-rc[2], sum(!o)), tolerance = 1e-8)
  }
})
