# ---------------------------------------------------------------------------
# Rosenbaum-type sensitivity analysis for matched pairs
#
# Under the sensitivity model of Rosenbaum (2002, Ch. 4), two units in a
# matched pair may differ in their odds of treatment by at most a factor
# Gamma >= 1. For each Gamma, the one-sided P-value of a randomization test
# has a sharp upper bound. The path stores, for each Gamma, the standard
# normal deviate equivalent to that upper bound, so the generic functions
# (robustness_points(), plot_robustness_curve(), sens_report(), sens_table())
# apply unchanged.
# ---------------------------------------------------------------------------

# Internal: upper-bound deviate for Wilcoxon's signed rank statistic.
# Large-sample approximation of Rosenbaum (2002, Sec. 4.3.3), computed as in
# DOS2::senWilcox(): all |d| are ranked (ties receive average ranks), then
# zero differences receive rank zero. Zeros are therefore not dropped before
# ranking.
.wilcox_dev <- function(d, gamma) {
  ad <- abs(d)
  rk <- rank(ad) * (ad > 0)
  ts <- sum((d > 0) * rk)
  vapply(gamma, function(g) {
    pr <- g / (1 + g)
    (ts - pr * sum(rk)) / sqrt(sum(rk^2) * pr * (1 - pr))
  }, numeric(1))
}

# Internal: exact upper-bound P-value for matched pairs with a binary
# outcome (McNemar's statistic; Rosenbaum 2002, Sec. 4.3.2). Among the
# discordant pairs, the number in which the treated unit had the event is
# stochastically no larger than Binomial(D, Gamma / (1 + Gamma)).
.binary_p <- function(d, gamma) {
  n_disc <- sum(d != 0)
  ts <- sum(d > 0)
  vapply(gamma, function(g)
    stats::pbinom(ts - 1, n_disc, g / (1 + g), lower.tail = FALSE),
    numeric(1))
}

# Internal: P-value to an upper-tail standard normal deviate, kept finite.
.p_to_dev <- function(p) {
  stats::qnorm(pmin(pmax(p, 1e-300), 1 - 1e-15), lower.tail = FALSE)
}


#' Rosenbaum sensitivity path for matched pairs
#'
#' Computes a sensitivity path for a matched-pair observational study under
#' the sensitivity model of Rosenbaum (2002): within a pair, the odds of
#' treatment for the two units may differ by at most a factor
#' \eqn{\Gamma \ge 1}{Gamma >= 1}. For each value of \eqn{\Gamma}{Gamma}, the
#' function computes the sharp upper bound on the one-sided P-value of a
#' randomization test of no effect and stores it as an equivalent standard
#' normal deviate, so that the path can be drawn, located, tabulated, and
#' reported with the same functions as the other frameworks.
#'
#' Two statistics are available. `method = "wilcoxon"` uses Wilcoxon's signed
#' rank statistic with the large-sample approximation of Rosenbaum (2002,
#' Sec. 4.3.3), computed as in `DOS2::senWilcox()`: all absolute differences
#' are ranked (ties receive average ranks) and zero differences then receive
#' rank zero. When some differences are zero, the statistic and P-values
#' therefore differ from those of `wilcox.test()` and `rbounds::psens()`,
#' which drop zeros before ranking. `method = "binary"` is for a binary outcome
#' (differences in \{-1, 0, 1\}) and uses the exact binomial bound for
#' McNemar's statistic (Rosenbaum 2002, Sec. 4.3.2), which agrees with
#' `sensitivity2x2xk::mh()` applied to the pairs as strata.
#'
#' The test is one-sided in the direction of `d > 0` (treated minus control).
#' Supply `-d` to test the opposite direction.
#'
#' @param d Numeric vector of treated-minus-control differences, one per
#'   matched pair. For `method = "binary"`, values must be in \{-1, 0, 1\}.
#' @param gamma Numeric vector of sensitivity parameters, all `>= 1`.
#' @param method `"wilcoxon"` (default) or `"binary"`.
#' @param alpha One-sided significance level used for the stored threshold
#'   and the sensitivity value.
#'
#' @return A `confoundsens` object with `lambda` = \eqn{\Gamma}{Gamma} and
#'   `theta` = the standard normal deviate equivalent to the upper-bound
#'   P-value. `meta$threshold` is `qnorm(1 - alpha)`, so
#'   [robustness_points()] reports as `lambda_sig` the value of
#'   \eqn{\Gamma}{Gamma} at which the upper bound reaches `alpha`, and as
#'   `lambda_null` the value at which it reaches 0.5. `meta$sens_value` holds
#'   the sensitivity value computed by root finding rather than on the grid,
#'   and `meta$pvalue` the upper-bound P-values along the grid.
#'
#' @references Rosenbaum, P. R. (2002). *Observational Studies* (2nd ed.).
#'   Springer. \doi{10.1007/978-1-4757-3692-2}
#' @seealso [robustness_points()], [sens_table()]
#' @export
#'
#' @examples
#' set.seed(1)
#' d <- rnorm(60, mean = 0.5)
#' p <- sens_path_pairs(d, gamma = seq(1, 4, by = 0.01))
#' robustness_points(p)
#' sens_report(p)
#'
#' # Binary outcome: 30 discordant pairs, 22 with the event in the treated unit
#' db <- c(rep(1, 22), rep(-1, 8), rep(0, 50))
#' robustness_points(sens_path_pairs(db, seq(1, 3, by = 0.01), method = "binary"))
sens_path_pairs <- function(d, gamma = seq(1, 3, by = 0.01),
                            method = c("wilcoxon", "binary"),
                            alpha = 0.05) {
  method <- match.arg(method)
  if (!is.numeric(d) || length(d) < 2L || anyNA(d)) {
    stop("`d` must be a numeric vector of pair differences (length >= 2, no NA).",
         call. = FALSE)
  }
  if (!is.numeric(gamma) || length(gamma) < 2L || anyNA(gamma) || any(gamma < 1)) {
    stop("`gamma` must be a numeric vector (length >= 2) with values >= 1.",
         call. = FALSE)
  }
  if (!is.numeric(alpha) || length(alpha) != 1L || alpha <= 0 || alpha >= 1) {
    stop("`alpha` must be a single number in (0, 1).", call. = FALSE)
  }
  if (method == "binary" && !all(d %in% c(-1, 0, 1))) {
    stop("For method = \"binary\", `d` must contain only -1, 0, and 1.",
         call. = FALSE)
  }
  if (all(d == 0)) {
    stop("All pair differences are zero; the test is undefined.", call. = FALSE)
  }

  gamma <- sort(unique(gamma))
  if (method == "wilcoxon") {
    dev_fun <- function(g) .wilcox_dev(d, g)
    stat <- sum((d > 0) * rank(abs(d)) * (abs(d) > 0))
    stat_label <- "Wilcoxon signed rank statistic"
  } else {
    dev_fun <- function(g) .p_to_dev(.binary_p(d, g))
    stat <- sum(d > 0)
    stat_label <- "Discordant pairs with the event in the treated unit"
  }
  dev  <- dev_fun(gamma)
  pval <- stats::pnorm(dev, lower.tail = FALSE)
  if (method == "binary") pval <- .binary_p(d, gamma)

  crit <- stats::qnorm(1 - alpha)
  sens_value <- if (dev_fun(1) <= crit) NA_real_ else {
    hi <- 2
    while (dev_fun(hi) > crit && hi < 1e6) hi <- hi * 2
    if (dev_fun(hi) > crit) Inf else
      stats::uniroot(function(g) dev_fun(g) - crit, c(1, hi), tol = 1e-10)$root
  }

  new_confoundsens(
    lambda = gamma, theta = dev,
    meta = list(
      framework    = "rosenbaum",
      source       = "pairs",
      method       = method,
      d            = d,
      n_pairs      = length(d),
      n_nonzero    = sum(d != 0),
      statistic    = stat,
      statistic_label = stat_label,
      estimate     = dev[1L],
      alpha        = alpha,
      sens_value   = sens_value,
      pvalue       = pval,
      null         = 0,
      threshold    = crit,
      lambda_label = "Gamma (odds of differential treatment assignment)",
      theta_label  = "Deviate equivalent to upper-bound P-value"
    )
  )
}
