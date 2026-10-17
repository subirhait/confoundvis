# Internal: first lambda at which g(lambda) <= 0, by linear interpolation.
.first_crossing <- function(lambda, g) {
  ok <- !is.na(g)
  lambda <- lambda[ok]; g <- g[ok]
  if (length(g) == 0L) return(NA_real_)
  if (g[1L] <= 0) return(lambda[1L])
  i <- which(g <= 0)
  if (length(i) == 0L) return(NA_real_)
  i <- i[1L]
  lambda[i - 1L] + (0 - g[i - 1L]) * (lambda[i] - lambda[i - 1L]) /
    (g[i] - g[i - 1L])
}


#' Locate robustness thresholds along a sensitivity path
#'
#' Finds, for each level of a `confoundsens` object, the confounder strength
#' \eqn{\lambda} at which (i) the point estimate reaches the null value and
#' (ii) statistical significance is lost. Crossings are located by linear
#' interpolation between grid points, so their precision depends on the
#' spacing of `lambda`.
#'
#' The significance criterion is chosen from what the object carries, in this
#' order: a critical value on the `theta` scale in `meta$threshold` (ITCV
#' paths); adjusted t statistics with residual degrees of freedom in
#' `meta$dof` (partial \eqn{R^2} paths; critical value from a t distribution
#' with `dof - 1` degrees of freedom, as in Cinelli and Hazlett, 2020);
#' standard errors (normal critical value); or t statistics without degrees
#' of freedom (normal critical value).
#'
#' @param x A `confoundsens` object.
#' @param alpha Two-sided significance level.
#'
#' @return A data frame with columns `level`, `theta_start`, `lambda_null`,
#'   `lambda_sig`, and `sig_rule`. `NA` means the crossing is not reached
#'   within the supplied `lambda` range.
#' @export
#'
#' @examples
#' fit <- lm(mpg ~ am + wt + hp, data = mtcars)
#' robustness_points(sens_path_lm(fit, "am", lambda = seq(0, 0.9, by = 0.001)))
robustness_points <- function(x, alpha = 0.05) {
  if (!inherits(x, "confoundsens")) {
    stop("`x` must be an object of class 'confoundsens'.", call. = FALSE)
  }
  x <- validate_confoundsens(x)
  meta <- x$meta
  null <- if (is.null(meta$null)) 0 else meta$null
  lev  <- if (is.null(x$level)) rep("pooled", length(x$lambda)) else x$level

  res <- lapply(unique(lev), function(L) {
    i  <- which(lev == L)
    o  <- order(x$lambda[i])
    i  <- i[o]
    la <- x$lambda[i]
    th <- x$theta[i]
    # Direction of the original conclusion: from the stored threshold when
    # there is one (ITCV, Rosenbaum), otherwise from the starting value.
    s0 <- if (!is.null(meta$threshold) && meta$threshold != null)
      sign(meta$threshold - null) else sign(th[1L] - null)
    if (s0 == 0) s0 <- 1
    lam_null <- .first_crossing(la, s0 * (th - null))

    rule <- NA_character_
    lam_sig <- NA_real_
    if (!is.null(meta$threshold)) {
      lam_sig <- .first_crossing(la, s0 * (th - meta$threshold))
      rule <- "theta threshold"
    } else if (!is.null(x$t) && !is.null(meta$dof) && is.finite(meta$dof)) {
      crit <- stats::qt(1 - alpha / 2, meta$dof - 1)
      lam_sig <- .first_crossing(la, s0 * x$t[i] - crit)
      rule <- sprintf("t, df = %s", meta$dof - 1)
    } else if (!is.null(x$se)) {
      z <- stats::qnorm(1 - alpha / 2)
      lam_sig <- .first_crossing(la, s0 * (th - null) - z * x$se[i])
      rule <- "normal CI"
    } else if (!is.null(x$t)) {
      lam_sig <- .first_crossing(la, s0 * x$t[i] - stats::qnorm(1 - alpha / 2))
      rule <- "normal t"
    }
    data.frame(level = L, theta_start = th[1L], lambda_null = lam_null,
               lambda_sig = lam_sig, sig_rule = rule,
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, res)
  rownames(out) <- NULL
  out
}


#' Plain-language sensitivity report
#'
#' Summarises a `confoundsens` path in a few sentences: the starting
#' estimate, the framework-specific summary (robustness value, ITCV, or
#' E-value) when available, and the confounder strengths at which the point
#' estimate reaches the null and significance is lost. The report closes
#' with a statement of what the analysis cannot establish.
#'
#' @inheritParams robustness_points
#' @param digits Number of significant digits.
#'
#' @return An object of class `confoundvis_report`: a list with `text` (a
#'   character vector) and `points` (the output of [robustness_points()]).
#'   It prints as text.
#' @export
#'
#' @examples
#' fit <- lm(mpg ~ am + wt + hp, data = mtcars)
#' sens_report(sens_path_lm(fit, "am", lambda = seq(0, 0.9, by = 0.001)))
sens_report <- function(x, alpha = 0.05, digits = 3) {
  pts  <- robustness_points(x, alpha = alpha)
  meta <- x$meta
  fmt  <- function(v) formatC(v, digits = digits, format = "g")
  lam_lab <- if (is.null(meta$lambda_label)) "lambda" else meta$lambda_label
  max_l <- max(x$lambda)

  txt <- sprintf("Sensitivity report (framework: %s%s)",
                 if (is.null(meta$framework)) "unspecified" else meta$framework,
                 if (is.null(meta$source)) "" else paste0("; source: ", meta$source))
  if (!is.null(meta$treatment) && !is.na(meta$treatment)) {
    txt <- c(txt, sprintf("Treatment: %s", meta$treatment))
  }
  fw <- if (is.null(meta$framework)) "" else meta$framework
  if (fw == "partial_r2" && !is.null(meta$rv)) {
    txt <- c(txt, sprintf(
      "Estimate %s (SE %s, df %s); robustness value RV = %s.",
      fmt(meta$estimate), fmt(meta$se0), meta$dof, fmt(meta$rv)))
  } else if (fw == "itcv" && !is.null(meta$itcv)) {
    txt <- c(txt, sprintf(
      "Partial correlation r = %s; critical r = %s; ITCV = %s.",
      fmt(meta$estimate), fmt(meta$r_crit), fmt(meta$itcv)))
  } else if (fw == "evalue" && !is.null(meta$evalue)) {
    txt <- c(txt, sprintf("Risk ratio %s; E-value (point) = %s.",
                          fmt(meta$estimate), fmt(meta$evalue)))
  } else if (fw == "rosenbaum" && !is.null(meta$statistic)) {
    txt <- c(txt, sprintf(
      "Matched pairs: %s (%s with nonzero difference); %s = %s; sensitivity value %s (one-sided alpha = %s).",
      meta$n_pairs, meta$n_nonzero, meta$statistic_label,
      format(signif(meta$statistic, 6)),
      if (is.na(meta$sens_value)) "not reached (not significant at Gamma = 1)"
      else if (!is.finite(meta$sens_value)) "Gamma > 1e6"
      else paste("Gamma =", fmt(meta$sens_value)), meta$alpha))
  }
  txt <- c(txt, sprintf("Confounder strength (lambda): %s.", lam_lab))

  for (k in seq_len(nrow(pts))) {
    p  <- pts[k, ]
    pre <- if (nrow(pts) > 1L) sprintf("[%s] ", p$level) else ""
    if (fw == "rosenbaum") {
      a1 <- if (is.null(meta$alpha)) alpha else meta$alpha
      r1 <- if (is.na(p$lambda_sig))
        sprintf("The upper-bound P-value stays below %s for Gamma <= %s", a1, fmt(max_l))
      else if (p$lambda_sig <= min(x$lambda))
        sprintf("The test is not significant (one-sided alpha = %s) even at Gamma = %s",
                a1, fmt(min(x$lambda)))
      else sprintf("The upper-bound P-value reaches %s at Gamma = %s", a1, fmt(p$lambda_sig))
      r2 <- if (is.na(p$lambda_null))
        sprintf("it stays below 0.5 for Gamma <= %s", fmt(max_l))
      else if (p$lambda_null <= min(x$lambda))
        sprintf("it is already at least 0.5 at Gamma = %s", fmt(min(x$lambda)))
      else sprintf("it reaches 0.5 at Gamma = %s", fmt(p$lambda_null))
      txt <- c(txt, sprintf("%s%s; %s.", pre, r1, r2))
      next
    }
    s1 <- if (is.na(p$lambda_null))
      sprintf("the point estimate does not reach the null for lambda <= %s", fmt(max_l))
    else sprintf("the point estimate reaches the null at lambda = %s", fmt(p$lambda_null))
    s2 <- if (is.na(p$sig_rule)) "no uncertainty information is available"
    else if (is.na(p$lambda_sig))
      sprintf("significance (alpha = %s) holds for lambda <= %s", alpha, fmt(max_l))
    else if (p$lambda_sig <= min(x$lambda)) {
      p0 <- if (fw == "partial_r2" && !is.null(meta$se0) && !is.null(meta$dof))
        2 * stats::pt(-abs(meta$estimate / meta$se0), meta$dof) else NA_real_
      if (!is.na(p0) && p0 < alpha)
        sprintf(paste0(
          "the unadjusted estimate is significant (t = %s, p = %s), but only ",
          "marginally: with the df - 1 adjustment used for the robustness value ",
          "(as in sensemakr), it is below the critical value at lambda = 0, so ",
          "any confounding removes significance (RV at alpha = %s is 0)"),
          fmt(meta$estimate / meta$se0), fmt(p0), alpha)
      else
        sprintf("the estimate is not significant (alpha = %s) even without confounding", alpha)
    }
    else sprintf("significance (alpha = %s) is lost at lambda = %s", alpha, fmt(p$lambda_sig))
    txt <- c(txt, sprintf("%sStarting at %s, %s; %s.", pre,
                          fmt(p$theta_start), s1, s2))
  }
  bench <- if (fw %in% c("partial_r2", "itcv"))
    "substantive knowledge and observed benchmarks (see covariate_impacts())."
  else "substantive knowledge and, where available, observed benchmarks."
  txt <- c(txt, paste(
    "These values describe how strong an omitted confounder would need to be.",
    "They do not show whether such a confounder exists, and they rest on the",
    "framework's own modelling assumptions. Judge plausibility against",
    bench))

  structure(list(text = txt, points = pts), class = "confoundvis_report")
}

#' @export
print.confoundvis_report <- function(x, ...) {
  cat(strwrap(x$text, exdent = 2), sep = "\n")
  invisible(x)
}
