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

test_that("infeasible additive scenarios are rejected", {
  f <- toy_fit()
  r <- sens_risk_difference(f, "a", p1 = 0.4, p0 = 0.2,
                            confounder_effect = "risk_difference")
  expect_equal(r$meta$max_feasible, 0.4, tolerance = 1e-8)   # 0.16 / 0.4
  expect_lte(max(r$lambda), r$meta$max_feasible)
  expect_warning(r2 <- sens_risk_difference(f, "a", p1 = 0.4, p0 = 0.2,
                                            confounder_effect = "risk_difference",
                                            lambda = c(0, 0.3, 0.5)),
                 "outside")
  expect_true(is.na(r2$theta[3]))
  expect_false(anyNA(r2$theta[1:2]))
})

test_that("thresholds beyond zero, offsets, and the confidence level are honoured", {
  f <- toy_fit()
  r <- sens_risk_difference(f, "a", p1 = 0.4, p0 = 0.2, threshold = -0.02)
  expect_gt(r$meta$lambda_threshold, r$meta$lambda_zero)
  expect_equal(sens_risk_difference(f, "a", 0.4, 0.2, threshold = -0.02,
                                    lambda = c(1, r$meta$lambda_threshold))$theta[2],
               -0.02, tolerance = 1e-6)
  set.seed(1); n <- 2000
  d <- data.frame(x = rnorm(n), a = rbinom(n, 1, 0.5), o = rnorm(n))
  d$y <- rbinom(n, 1, plogis(-1 + 0.8 * d$a + d$o + 0.3 * d$x))
  fo <- glm(y ~ a + x + offset(o), family = binomial, data = d)
  ro <- sens_risk_difference(fo, "a", 0.4, 0.2)
  direct <- mean(predict(fo, transform(d, a = 1), type = "response") -
                 predict(fo, transform(d, a = 0), type = "response"))
  expect_equal(ro$theta[1], direct, tolerance = 1e-10)
  r1 <- sens_risk_difference(f, "a", 0.4, 0.2, threshold = 0.02, alpha = 0.01)
  g <- plot_robustness_curve(r1, points = FALSE)
  k <- which(vapply(g$layers, function(l) inherits(l$geom, "GeomRibbon"), logical(1)))
  b <- ggplot2::layer_data(g, k)
  expect_equal(b$ymax[1] - r1$theta[1], stats::qnorm(0.995) * r1$se[1],
               tolerance = 1e-8)
  expect_match(paste(sens_report(r1)$text, collapse = " "), "alpha = 0.01")
  expect_false(any(grepl("alpha = 0.05", sens_report(r1)$text)))
})

test_that("prevalences at 0 and 1 give the closed-form answers", {
  f <- toy_fit()
  # additive, p1 = 1, p0 = 0: RD - delta; feasible while delta <= r1 = 0.16
  r <- sens_risk_difference(f, "a", 1, 0, confounder_effect = "risk_difference",
                            threshold = 0.02)
  expect_equal(r$meta$lambda_zero, 0.04, tolerance = 1e-7)
  expect_equal(r$meta$lambda_threshold, 0.02, tolerance = 1e-7)
  expect_equal(r$meta$max_feasible, 0.16, tolerance = 1e-7)
  # odds ratio, p1 = 1, p0 = 0 (equal groups, q = 0.5):
  # RD* = 0.02 + 0.5 [expit(logit .16 - g) - expit(logit .12 + g)]
  ro <- sens_risk_difference(f, "a", 1, 0)
  g <- log(ro$meta$lambda_zero)
  expect_equal(0.02 + 0.5 * (plogis(qlogis(0.16) - g) - plogis(qlogis(0.12) + g)),
               0, tolerance = 1e-7)
  # swapping the coding of U gives the same answer
  expect_equal(sens_risk_difference(f, "a", 0, 1)$meta$lambda_zero,
               ro$meta$lambda_zero, tolerance = 1e-7)
})

test_that("report wording for thresholds already met or never reached", {
  f <- toy_fit()
  t1 <- paste(sens_report(sens_risk_difference(f, "a", 0.4, 0.2, threshold = 0.05))$text,
              collapse = " ")
  expect_match(t1, "already no\\s+farther from zero than the practical threshold")
  expect_false(grepl("odds ratio of\\s+1 brings", t1))
  t2 <- paste(sens_report(sens_risk_difference(f, "a", 0.25, 0.2, threshold = 0.02))$text,
              collapse = " ")
  expect_match(t2, "No odds ratio brings the risk\\s+difference to zero")
  expect_false(grepl("not reached brings", t2))
})

test_that("additive feasibility is the average of person-level limits", {
  # one person with a fitted risk near 0 must not rule out every scenario
  set.seed(9); n <- 1500
  x <- rnorm(n); x[1] <- -12
  a <- rbinom(n, 1, 0.5)
  y <- rbinom(n, 1, plogis(-0.5 + 0.6 * a + 0.8 * x))
  f <- glm(y ~ a + x, family = binomial, data = data.frame(y, a, x))
  r <- sens_risk_difference(f, "a", 0.4, 0.2, confounder_effect = "risk_difference")
  R <- r$meta$rd_setup
  U <- pmin(R$r1 / 0.4, (1 - R$r1) / 0.6, R$r0 / 0.2, (1 - R$r0) / 0.8)
  expect_equal(r$meta$max_feasible, mean(U), tolerance = 1e-10)
  expect_gt(r$meta$max_feasible, 0.1)
  expect_lt(min(U), 1e-3)
  # endpoint prevalences: a zero denominator imposes no limit
  r1 <- sens_risk_difference(f, "a", 1, 0, confounder_effect = "risk_difference")
  expect_equal(r1$meta$max_feasible, mean(pmin(R$r1, 1 - R$r0)), tolerance = 1e-10)
  r0 <- sens_risk_difference(f, "a", 0, 1, confounder_effect = "risk_difference")
  expect_equal(r0$meta$max_feasible, mean(pmin(R$r1, 1 - R$r0)), tolerance = 1e-10)
  expect_false(anyNA(r1$theta))
})

test_that("baseline interval past the threshold is told apart from overlap", {
  f <- toy_fit()
  past <- paste(sens_report(sens_risk_difference(f, "a", 0.4, 0.2, threshold = 0.2))$text,
                collapse = " ")
  expect_match(past, "entirely past the\\s+threshold")
  over <- paste(sens_report(sens_risk_difference(f, "a", 0.4, 0.2, threshold = 0.04))$text,
                collapse = " ")
  expect_match(over, "includes the\\s+threshold")
})

test_that("surface draws a threshold on the far side of zero", {
  r <- sens_risk_difference(toy_fit(), "a", 0.4, 0.2, threshold = -0.02)
  s <- sens_surface(r, a = seq(0.3, 0.6, by = 0.025), b = seq(1, 20, by = 0.5))
  expect_true(any(grepl("^Practical threshold", s$boundary$kind)))
  expect_true(any(s$grid$zone == "Past practical threshold", na.rm = TRUE))
})
