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

# Internal: TRUE when the significance rule of a path is fixed at the level
# it was built with (a stored threshold, or stored confidence limits).
.fixed_level <- function(meta) {
  !is.null(meta$threshold) || !is.null(meta$ci_near)
}

# Internal: the significance level a path is evaluated at. Paths with a
# stored threshold carry the level it was computed from; a different
# requested level is refused rather than silently ignored.
.resolve_alpha <- function(meta, alpha) {
  if (!is.null(alpha) &&
      (!is.numeric(alpha) || length(alpha) != 1L || is.na(alpha) ||
       alpha <= 0 || alpha >= 1)) {
    stop("`alpha` must be NULL or a single number in (0, 1).", call. = FALSE)
  }
  if (!.fixed_level(meta)) {
    if (!is.null(alpha)) return(alpha)
    return(if (is.null(meta$default_alpha)) 0.05 else meta$default_alpha)
  }
  stored <- if (is.null(meta$alpha)) NA_real_ else as.numeric(meta$alpha)
  if (is.null(alpha)) return(stored)
  if (is.na(stored)) {
    warning("This path has a fixed threshold that was not built from a ",
            "significance level (substantive threshold or imported without a ",
            "level); `alpha` is ignored.", call. = FALSE)
    return(NA_real_)
  }
  if (abs(alpha - stored) > 1e-12) {
    fw <- if (is.null(meta$framework)) "" else meta$framework
    how <- switch(fw,
      itcv      = "itcv_lm(..., alpha = %s) or from_konfound(..., alpha = %s)",
      rir       = "rir_lm(..., alpha = %s) or from_konfound(..., alpha = %s)",
      rosenbaum = "sens_path_pairs(..., alpha = %s)%s",
      iv        = "iv.sensemakr::sensemakr(..., alpha = %s)%s, then from_iv_sensemakr()",
      "the function that created it with alpha = %s%s")
    how <- if (!fw %in% c("itcv", "rir"))
      sprintf(how, alpha, "") else sprintf(how, alpha, alpha)
    stop(sprintf(paste0(
      "This %s path was built with alpha = %s, and its threshold is fixed at ",
      "that level. Rebuild it at alpha = %s with %s."),
      if (fw == "") "sensitivity" else fw, stored, alpha, how), call. = FALSE)
  }
  stored
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
#' order: bias-adjusted confidence limits in `meta$ci_near` (instrumental
#' variable paths from [from_iv_sensemakr()]); a critical value on the
#' `theta` scale in `meta$threshold` (ITCV paths); adjusted t statistics with residual degrees of freedom in
#' `meta$dof` (partial \eqn{R^2} paths; critical value from a t distribution
#' with `dof - 1` degrees of freedom, as in Cinelli and Hazlett, 2020);
#' standard errors (normal critical value); or t statistics without degrees
#' of freedom (normal critical value).
#'
#' @param x A `confoundsens` object.
#' Paths with a stored threshold (ITCV, RIR, and Rosenbaum paths) fix the
#' significance level when they are built, because the threshold is computed
#' from it. For these paths `alpha` defaults to the stored level. Asking for a
#' different level is an error: rebuild the path at the new level instead
#' (for example `itcv_lm(model, treatment, alpha = 0.01)`). Instrumental
#' variable paths are fixed at the level of the `iv.sensemakr` object in the
#' same way. For other paths the default is the level stored in
#' `meta$default_alpha` ([sens_dml()] stores the two-sided level that matches
#' the one-sided `doubleml` level), or 0.05. The level actually used is returned in column `alpha`.
#'
#' @param alpha Significance level, or `NULL` (default) to use the level
#'   stored in the path, if any, and 0.05 otherwise. Two-sided, except for
#'   Rosenbaum paths, which use the one-sided level of Rosenbaum's test.
#'
#' @return A data frame with columns `level`, `theta_start`, `lambda_null`,
#'   `lambda_sig`, `sig_rule`, and `alpha`. `NA` means the crossing is not
#'   reached within the supplied `lambda` range.
#' @export
#'
#' @examples
#' fit <- lm(mpg ~ am + wt + hp, data = mtcars)
#' robustness_points(sens_path_lm(fit, "am", lambda = seq(0, 0.9, by = 0.001)))
robustness_points <- function(x, alpha = NULL) {
  if (!inherits(x, "confoundsens")) {
    stop("`x` must be an object of class 'confoundsens'.", call. = FALSE)
  }
  x <- validate_confoundsens(x)
  meta <- x$meta
  alpha <- .resolve_alpha(meta, alpha)
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
    if (!is.null(meta$ci_near)) {
      lam_sig <- .first_crossing(la, s0 * (meta$ci_near[i] - null))
      rule <- "adjusted confidence limit"
    } else if (!is.null(meta$threshold)) {
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
               lambda_sig = lam_sig, sig_rule = rule, alpha = alpha,
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
sens_report <- function(x, alpha = NULL, digits = 3) {
  pts  <- robustness_points(x, alpha = alpha)
  alpha <- pts$alpha[1L]
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
  } else if (fw == "rir" && !is.null(meta$beta_threshold)) {
    txt <- c(txt, if (!is.na(meta$perc_bias)) sprintf(
      "Estimate %s; threshold %s. To nullify the inference, %s%% of the estimate would have to be due to bias%s.",
      fmt(meta$estimate), fmt(meta$beta_threshold),
      formatC(meta$perc_bias, digits = 4, format = "g"),
      if (!is.na(meta$rir)) {
        if (!is.null(meta$n_obs) && !is.na(meta$n_obs))
          sprintf(": replacing %s of %s cases with cases for which the effect is zero (RIR = %s)",
                  meta$rir, meta$n_obs, meta$rir)
        else sprintf(" (RIR = %s cases replaced with null-effect cases)", meta$rir)
      } else "")
    else sprintf(
      "Estimate %s does not exceed the threshold %s; see konfound for the replacement needed to sustain an inference.",
      fmt(meta$estimate), fmt(meta$beta_threshold)))
  } else if (fw == "iv" && !is.null(meta$rv)) {
    txt <- c(txt, sprintf(paste0(
      "IV estimate %s (Anderson-Rubin %s%% CI [%s, %s]); robustness value ",
      "RV = %s; extreme robustness value XRV = %s (alpha = %s)."),
      fmt(meta$estimate), 100 * (1 - meta$alpha), fmt(meta$ci_unadj[1]),
      fmt(meta$ci_unadj[2]), fmt(meta$rv), fmt(meta$xrv), meta$alpha))
  } else if (fw == "dml" && !is.null(meta$rv)) {
    txt <- c(txt, sprintf(paste0(
      "DoubleML %s (%s score) estimate %s; bias scale S = %s; rho = %s. ",
      "RV = %s; RVa = %s at one-sided level %s (two-sided alpha = %s)."),
      toupper(meta$model), meta$score, fmt(meta$estimate), fmt(meta$max_bias),
      meta$rho, if (is.na(meta$rv)) "not reached" else fmt(meta$rv),
      if (is.na(meta$rva)) "not reached" else fmt(meta$rva), meta$level,
      meta$default_alpha))
  } else if (fw == "risk_difference" && !is.null(meta$estimate)) {
    pp <- function(v) paste0(formatC(100 * v, format = "f", digits = 1), "%")
    st <- function(v) if (is.na(v)) "not reached" else fmt(v)
    eff <- if (meta$confounder_effect == "odds_ratio") "odds ratio" else "risk difference"
    txt <- c(txt,
      sprintf(paste0(
        "Standardized risks (%s): %s treated vs %s untreated; risk difference ",
        "%s (SE %s)."), meta$target, pp(meta$risk1_obs), pp(meta$risk0_obs),
        fmt(meta$estimate), fmt(meta$se0)),
      sprintf(paste0(
        "Hypothetical binary confounder: prevalence %s among the treated and ",
        "%s among the untreated. Under this scenario, a confounder-outcome %s ",
        "of %s brings the risk difference to the practical threshold %s, and ",
        "%s brings it to zero. %s"),
        pp(meta$p1), pp(meta$p0), eff, st(meta$lambda_threshold),
        fmt(meta$rd_threshold), st(meta$lambda_zero), {
          lc <- .rd_ci_cross_meta(meta, alpha)
          lo <- if (meta$confounder_effect == "odds_ratio") 1 else 0
          if (is.na(lc)) sprintf(paste0(
            "Allowing for sampling error (alpha = %s), the confidence limit ",
            "nearest zero does not reach the threshold."), alpha)
          else if (lc <= lo) sprintf(paste0(
            "Allowing for sampling error (alpha = %s), the confidence interval ",
            "already includes the threshold without confounding."), alpha)
          else sprintf(paste0(
            "Allowing for sampling error (alpha = %s), the confidence limit ",
            "nearest zero reaches the threshold at %s."), alpha, fmt(lc))
        }),
      if (meta$confounder_effect == "risk_difference" && is.finite(meta$max_feasible))
        sprintf(paste0("Scenarios with a confounder-outcome risk difference ",
                       "above %s imply risks outside [0, 1] and are not used."),
                fmt(meta$max_feasible)),
      paste("The path shows the range across confounding scenarios, not",
            "sampling error; standard errors are given for each scenario."))
  } else if (fw == "rosenbaum" && !is.null(meta$statistic)) {
    txt <- c(txt, sprintf(
      "Matched pairs: %s (%s with nonzero difference); %s = %s; sensitivity value %s (one-sided alpha = %s).",
      meta$n_pairs, meta$n_nonzero, meta$statistic_label,
      format(signif(meta$statistic, 6)),
      if (is.na(meta$sens_value)) "not reached (not significant at Gamma = 1)"
      else if (!is.finite(meta$sens_value)) "Gamma > 1e6"
      else paste("Gamma =", fmt(meta$sens_value)), meta$alpha))
  }
  txt <- c(txt, sprintf("%s (lambda): %s.",
                        if (fw == "rir") "Sensitivity parameter" else "Confounder strength",
                        lam_lab))

  for (k in seq_len(nrow(pts))) {
    p  <- pts[k, ]
    pre <- if (nrow(pts) > 1L) sprintf("[%s] ", p$level) else ""
    if (fw == "rir") {
      q1 <- if (is.na(p$lambda_sig))
        sprintf("the estimate stays above the threshold for lambda <= %s", fmt(max_l))
      else if (p$lambda_sig <= min(x$lambda))
        "the estimate does not exceed the threshold even before replacement"
      else sprintf("the estimate falls to the threshold when a proportion %s of the cases is replaced",
                   fmt(p$lambda_sig))
      txt <- c(txt, sprintf("%sStarting at %s, %s.", pre, fmt(p$theta_start), q1))
      next
    }
    if (fw == "rosenbaum") {
      a1 <- alpha
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
    what <- if (fw %in% c("iv", "dml")) "the bias-adjusted bound" else
      if (fw == "risk_difference") "the adjusted risk difference" else "the point estimate"
    s1 <- if (is.na(p$lambda_null))
      sprintf("%s does not reach the null for lambda <= %s", what, fmt(max_l))
    else sprintf("%s reaches the null at lambda = %s", what, fmt(p$lambda_null))
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
  if (fw == "rir") {
    txt <- c(txt, paste(
      "RIR describes how many cases would have to be replaced with cases in",
      "which the treatment has no effect. It does not show that such cases",
      "exist; judge plausibility against knowledge of the sample and design."))
    return(structure(list(text = txt, points = pts), class = "confoundvis_report"))
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
