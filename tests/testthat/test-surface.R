test_that("partial R2 surface boundaries match the RV and sensemakr", {
  fit <- lm(mpg ~ am + wt + hp, data = mtcars)
  p <- sens_path_lm(fit, "wt")
  s <- sens_surface(p, n = 101)
  expect_s3_class(s, "confoundsurf")
  f2 <- (p$meta$estimate / p$meta$se0)^2 / p$meta$dof
  nb <- s$boundary[grepl("^Null", s$boundary$kind), ]
  expect_equal(nb$b, f2 * (1 - nb$a) / nb$a)
  sb <- s$boundary[grepl("^Significance \\(alpha", s$boundary$kind), ]
  i <- which.min(abs(sb$a - sb$b))
  expect_equal(sb$a[i], robustness_points(sens_path_lm(fit, "wt",
                 lambda = seq(0, 0.9, by = 1e-4)))$lambda_sig, tolerance = 2e-3)
  skip_if_not_installed("sensemakr")
  sm <- sensemakr::sensemakr(fit, "wt")
  expect_equal(sb$a[i], as.numeric(sm$sensitivity_stats$rv_qa), tolerance = 2e-3)
})

test_that("E-value surface passes through (E, E) and (E_CI, E_CI)", {
  e <- structure(rbind(RR = c(point = 1.8, lower = 1.3, upper = 2.5),
                       `E-values` = c(point = 1.8 + sqrt(1.8 * 0.8),
                                      lower = 1.3 + sqrt(1.3 * 0.3), upper = NA)),
                 class = c("evalue", "matrix"))
  s <- sens_surface(from_evalue(e, lambda = seq(1, 4, by = 0.01)))
  bf <- function(a, b) a * b / (a + b - 1)
  nb <- s$boundary[grepl("^Null", s$boundary$kind), ]
  expect_equal(bf(nb$a, nb$b), rep(1.8, nrow(nb)), tolerance = 1e-10)
  cb <- s$boundary[s$boundary$kind == "Confidence limit", ]
  expect_equal(bf(cb$a, cb$b), rep(1.3, nrow(cb)), tolerance = 1e-10)
  E <- 1.8 + sqrt(1.8 * 0.8)
  expect_equal(bf(E, E), 1.8)
})

test_that("Rosenbaum surface reproduces sensitivitymv::amplify", {
  set.seed(9)
  p <- sens_path_pairs(rnorm(60, 0.6), gamma = seq(1, 6, by = 0.01))
  s <- sens_surface(p, n = 101)
  sb <- s$boundary[grepl("^Signif", s$boundary$kind), ]
  sv <- p$meta$sens_value
  lam <- sv + c(1, 2)
  ours <- (sv * lam - 1) / (lam - sv)
  expect_equal(stats::approx(sb$a, sb$b, xout = lam)$y, ours, tolerance = 1e-3)
  skip_if_not_installed("sensitivitymv")
  expect_equal(unname(sensitivitymv::amplify(sv, lam)), ours)
})

test_that("ITCV surface uses the exact boundary; RIR and bad input error", {
  fit <- lm(mpg ~ am + wt + hp, data = mtcars)
  it <- itcv_lm(fit, "wt")
  s <- sens_surface(it$path, n = 81)
  ob <- s$boundary[s$boundary$kind == "Significance (original direction)", ]
  radj <- (it$r - ob$a * ob$b) / sqrt((1 - ob$a^2) * (1 - ob$b^2))
  expect_equal(abs(radj), rep(it$r_crit, nrow(ob)), tolerance = 1e-8)
  expect_s3_class(plot_sensitivity_surface(s, covariate_impacts(fit, "wt")), "ggplot")
  expect_error(sens_surface(rir_lm(fit, "wt")$path), "single sensitivity parameter")
  expect_error(sens_surface(new_confoundsens(c(0, 1), c(1, 0))), "needs a path")
  expect_error(plot_sensitivity_surface(1), "confoundsurf")
  expect_output(print(s), "confoundsurf")
})

test_that("small user grids that miss a boundary do not fail", {
  fit <- lm(mpg ~ am + wt + hp, data = mtcars)
  expect_s3_class(sens_surface(sens_path_lm(fit, "wt"), a = c(0, 0.02), b = c(0, 0.02)),
                  "confoundsurf")
  e <- structure(rbind(RR = c(point = 1.8, lower = 1.3, upper = 2.5),
                       `E-values` = c(point = 3, lower = 1.92, upper = NA)),
                 class = c("evalue", "matrix"))
  expect_s3_class(sens_surface(from_evalue(e), a = c(1, 1.2), b = c(1, 1.2)), "confoundsurf")
  set.seed(9)
  p <- sens_path_pairs(rnorm(60, 0.6))
  expect_s3_class(sens_surface(p, a = c(1, 1.5), b = c(1, 1.5)), "confoundsurf")
})

test_that("Rosenbaum zones: P = 0.5 boundary by root finding, NA sens value", {
  set.seed(9)
  p <- sens_path_pairs(rnorm(60, 0.6))            # default grid 1..3
  s <- sens_surface(p)
  gh <- s$meta$gamma_half
  expect_true(is.finite(gh) && gh > p$meta$sens_value)
  expect_equal(confoundvis:::.wilcox_dev(p$meta$d, gh), 0, tolerance = 1e-8)
  expect_true(any(s$boundary$kind == "Upper-bound P-value = 0.5"))
  expect_true("Upper-bound P-value at least 0.5" %in% s$grid$zone)
  p0 <- sens_path_pairs(c(0.1, -0.2, 0.05, -0.1, 0.02, -0.03))  # not significant
  s0 <- sens_surface(p0)
  expect_true(is.na(p0$meta$sens_value))
  expect_true(all(!is.na(s0$grid$zone)))
  expect_warning(sens_surface(p, alpha = 0.1), "ignored")
})

test_that("E-value: CI including 1 and protective risk ratios", {
  e1 <- structure(rbind(RR = c(point = 1.3, lower = 0.9, upper = 1.9),
                        `E-values` = c(point = 1.92, lower = 1, upper = NA)),
                  class = c("evalue", "matrix"))
  s1 <- sens_surface(from_evalue(e1))
  expect_true(all(!is.na(s1$grid$zone)))
  e2 <- structure(rbind(RR = c(point = 0.6, lower = 0.4, upper = 0.8),
                        `E-values` = c(point = NA, lower = NA, upper = NA)),
                  class = c("evalue", "matrix"))
  s2 <- sens_surface(from_evalue(e2))
  bf <- function(a, b) a * b / (a + b - 1)
  cb <- s2$boundary[s2$boundary$kind == "Confidence limit", ]
  expect_equal(bf(cb$a, cb$b), rep(1 / 0.8, nrow(cb)), tolerance = 1e-10)
})

test_that("non-significant starts show where significance is gained", {
  fit <- lm(mpg ~ am + wt + hp, data = mtcars)    # 'am' is not significant
  s <- sens_surface(itcv_lm(fit, "am")$path, n = 81)
  expect_true("Significant, original direction" %in% s$grid$zone)
  s2 <- sens_surface(sens_path_lm(fit, "am"), n = 81)
  expect_true("Significant, original direction" %in% s2$grid$zone)
  s3 <- sens_surface(sens_path_lm(fit, "wt"), n = 81)
  expect_true(any(grepl("reversed", s3$boundary$kind)))
  expect_true(all(c("solid", "dashed") %in% s3$boundary$style))
})
