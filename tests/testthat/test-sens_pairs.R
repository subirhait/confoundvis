test_that("sens_path_pairs (wilcoxon) matches DOS2::senWilcox", {
  skip_if_not_installed("DOS2")
  set.seed(11)
  d <- round(rnorm(70, 0.35), 2)
  d[1:4] <- 0
  g <- c(1, 1.25, 1.5, 2, 3)
  p <- sens_path_pairs(d, gamma = g)
  ours <- stats::pnorm(p$theta, lower.tail = FALSE)
  ref  <- vapply(g, function(G) DOS2::senWilcox(d, gamma = G)$pval, numeric(1))
  expect_equal(ours, ref, tolerance = 1e-12)
  expect_equal(p$meta$pvalue, ref, tolerance = 1e-12)
})

test_that("sens_path_pairs (binary) matches sensitivity2x2xk::mh", {
  skip_if_not_installed("sensitivity2x2xk")
  d <- c(rep(1, 30), rep(-1, 14), rep(0, 40))
  g <- c(1, 1.3, 1.8, 2.5)
  p <- sens_path_pairs(d, gamma = g, method = "binary")
  tab <- array(0, c(2, 2, length(d)))
  for (k in seq_along(d)) {
    yt <- as.integer(d[k] == 1); yc <- as.integer(d[k] == -1)
    tab[1, 2 - yt, k] <- 1; tab[2, 2 - yc, k] <- 1
  }
  ref <- vapply(g, function(G) sensitivity2x2xk::mh(tab, Gamma = G)$pval,
                numeric(1))
  expect_equal(p$meta$pvalue, ref, tolerance = 1e-10)
})

test_that("binary bound equals the closed-form binomial tail", {
  d <- c(rep(1, 22), rep(-1, 8), rep(0, 5))
  G <- 1.7
  p <- sens_path_pairs(d, gamma = c(1, G), method = "binary")
  expect_equal(p$meta$pvalue[2],
               stats::pbinom(21, 30, G / (1 + G), lower.tail = FALSE))
})

test_that("sensitivity value is where the bound equals alpha", {
  set.seed(3)
  d <- rnorm(50, 0.6)
  p <- sens_path_pairs(d, gamma = seq(1, 6, by = 0.001), alpha = 0.05)
  sv <- p$meta$sens_value
  expect_true(is.finite(sv) && sv > 1)
  expect_equal(confoundvis:::.wilcox_dev(d, sv), stats::qnorm(0.95),
               tolerance = 1e-8)
  rp <- robustness_points(p)
  expect_equal(rp$lambda_sig, sv, tolerance = 1e-3)
  expect_identical(rp$sig_rule, "theta threshold")
})

test_that("sens_path_pairs handles non-significant data and bad input", {
  d <- c(0.1, -0.2, 0.05, -0.1, 0.02, -0.03)
  p <- sens_path_pairs(d)
  expect_true(is.na(p$meta$sens_value))
  r <- sens_report(p)
  expect_true(any(grepl("not significant", r$text)))
  expect_error(sens_path_pairs(c(1, NA)), "numeric vector")
  expect_error(sens_path_pairs(rnorm(5), gamma = c(0.5, 1)), ">= 1")
  expect_error(sens_path_pairs(c(1, 2, 0), method = "binary"), "-1, 0, and 1")
  expect_error(sens_path_pairs(c(0, 0, 0)), "zero")
})

test_that("sens_report describes a Rosenbaum path", {
  set.seed(5)
  r <- sens_report(sens_path_pairs(rnorm(40, 0.7), gamma = seq(1, 5, by = 0.01)))
  expect_true(any(grepl("sensitivity value Gamma", r$text)))
  expect_true(any(grepl("upper-bound P-value reaches 0.05", r$text)))
  expect_false(any(grepl("covariate_impacts", r$text)))
})

test_that("a negative-direction effect is reported as not significant", {
  set.seed(4)
  p <- sens_path_pairs(rnorm(40, -0.5), gamma = seq(1, 3, by = 0.01))
  rp <- robustness_points(p)
  expect_equal(rp$lambda_sig, 1)
  expect_true(is.na(p$meta$sens_value))
  r <- sens_report(p)
  expect_true(any(grepl("not significant", r$text)))
  expect_false(any(grepl("stays below 0.05", r$text)))
})
