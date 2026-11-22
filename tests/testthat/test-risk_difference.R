toy_fit <- function() {
  # risks 16% (treated) and 12% (untreated), equal groups, no covariates
  toy <- data.frame(a = rep(0:1, each = 500),
                    y = c(rep(1:0, c(60, 440)), rep(1:0, c(80, 420))))
  glm(y ~ a, family = binomial, data = toy)
}

test_that("odds-ratio model reproduces the worked example", {
  r <- sens_risk_difference(toy_fit(), "a", p1 = 0.4, p0 = 0.2,
                            threshold = 0.02, lambda = c(1, 2.16, 3, 4.47))
  expect_equal(r$theta[1], 0.04, tolerance = 1e-7)
  expect_equal(round(100 * r$theta, 1), c(4.0, 2.0, 1.1, 0.0))
  expect_equal(r$meta$lambda_threshold, 2.1575, tolerance = 1e-4)
  expect_equal(r$meta$lambda_zero, 4.4722, tolerance = 1e-4)
})

test_that("closed-form solution agrees with root finding", {
  r <- c(0.01, 0.2, 0.5, 0.93); p <- 0.35; G <- 3.7
  t <- confoundvis:::.rd_t(r, p, G)
  expect_equal(p * G * t / (1 + G * t) + (1 - p) * t / (1 + t), r, tolerance = 1e-12)
  t2 <- confoundvis:::.rd_t(r, p, 1 / G)
  expect_equal(p * t2 / G / (1 + t2 / G) + (1 - p) * t2 / (1 + t2), r, tolerance = 1e-12)
})

test_that("risk-difference model is the closed-form correction", {
  r <- sens_risk_difference(toy_fit(), "a", p1 = 0.4, p0 = 0.2,
                            confounder_effect = "risk_difference",
                            lambda = c(0, 0.05, 0.1))
  expect_equal(r$theta, 0.04 - c(0, 0.05, 0.1) * 0.2, tolerance = 1e-7)
  expect_equal(r$meta$lambda_zero, 0.2, tolerance = 1e-8)
})

test_that("no confounding gives the standardized risk difference and SE", {
  set.seed(3); n <- 1500
  d <- data.frame(x1 = rnorm(n), x2 = factor(sample(letters[1:3], n, TRUE)))
  d$a <- rbinom(n, 1, plogis(0.4 * d$x1))
  d$y <- rbinom(n, 1, plogis(-1 + 0.5 * d$a + 0.7 * d$x1 + 0.2 * d$a * d$x1))
  f <- glm(y ~ a * x1 + x2, family = binomial, data = d)
  r <- sens_risk_difference(f, "a", p1 = 0.5, p0 = 0.3)
  p1 <- predict(f, transform(d, a = 1), type = "response")
  p0 <- predict(f, transform(d, a = 0), type = "response")
  expect_equal(r$theta[1], mean(p1 - p0), tolerance = 1e-7)
  rt <- sens_risk_difference(f, "a", p1 = 0.5, p0 = 0.3, target = "ATT")
  expect_equal(rt$theta[1], mean((p1 - p0)[d$a == 1]), tolerance = 1e-7)
  skip_if_not_installed("marginaleffects")
  me <- marginaleffects::avg_comparisons(f, variables = "a")
  expect_equal(r$theta[1], me$estimate, tolerance = 1e-8)
  expect_equal(r$se[1], me$std.error, tolerance = 1e-5)
  me2 <- marginaleffects::avg_comparisons(f, variables = "a",
                                          newdata = subset(d, a == 1))
  expect_equal(rt$se[1], me2$std.error, tolerance = 1e-5)
})

test_that("direction, crossings, factor treatment, and errors", {
  set.seed(5); n <- 3000
  x <- rnorm(n); a <- rbinom(n, 1, plogis(0.5 * x))
  y <- rbinom(n, 1, plogis(-0.5 - 0.6 * a + 0.5 * x))
  d <- data.frame(y, a, x, af = factor(ifelse(a == 1, "yes", "no")))
  f <- glm(y ~ a + x, family = binomial, data = d)
  r <- sens_risk_difference(f, "a", p1 = 0.8, p0 = 0.1, threshold = -0.03)
  expect_lt(r$theta[1], 0)
  expect_true(all(diff(r$theta[r$theta < 0]) > 0))     # moves toward zero
  expect_equal(r$meta$direction, -1)
  rp <- robustness_points(r)
  expect_equal(rp$lambda_null, r$meta$lambda_zero, tolerance = 1e-2)
  expect_lt(r$meta$lambda_threshold, r$meta$lambda_zero)
  # a weaker prevalence gap cannot remove this estimate: reported as NA
  expect_true(is.na(sens_risk_difference(f, "a", p1 = 0.3, p0 = 0.1)$meta$lambda_zero))
  ff <- glm(y ~ af + x, family = binomial, data = d)
  rf <- sens_risk_difference(ff, "af", p1 = 0.3, p0 = 0.1)
  expect_equal(rf$theta, sens_risk_difference(f, "a", p1 = 0.3, p0 = 0.1)$theta)
  expect_error(sens_risk_difference(lm(y ~ a, d), "a", 0.3, 0.1), "binomial")
  expect_error(sens_risk_difference(f, "a", 0.3, 0.3), "equals")
  expect_error(sens_risk_difference(f, "z", 0.3, 0.1), "treatment")
  expect_error(sens_risk_difference(f, "a", 1.3, 0.1), "p1")
})

test_that("surface, report, table, plot, and confoundvis()", {
  r <- sens_risk_difference(toy_fit(), "a", p1 = 0.4, p0 = 0.2, threshold = 0.02)
  s <- sens_surface(r, n = 31)
  expect_s3_class(s, "confoundsurf")
  expect_true(any(grepl("Practical threshold", s$boundary$kind)))
  expect_s3_class(plot_sensitivity_surface(s), "ggplot")
  expect_s3_class(plot_robustness_curve(r), "ggplot")
  txt <- paste(sens_report(r)$text, collapse = " ")
  expect_match(txt, "practical threshold")
  expect_match(txt, "not\\s+sampling error")
  tb <- sens_table(rd = r)
  expect_equal(tb$summary_value, r$meta$lambda_threshold)
  skip_if_not_installed("konfound")
  set.seed(2); n <- 800
  d <- data.frame(x = rnorm(n)); d$a <- rbinom(n, 1, 0.5)
  d$y <- rbinom(n, 1, plogis(-1 + 0.8 * d$a + 0.5 * d$x))
  f <- glm(y ~ a + x, family = binomial, data = d)
  res <- suppressWarnings(confoundvis(f, "a", risk = list(p1 = 0.4, p0 = 0.2,
                                                           threshold = 0.05)))
  expect_named(res$paths, c("rir", "risk_difference"))
  expect_output(print(res), "Risk difference")
})
