# ---------------------------------------------------------------------------
# Instrumental variables: sensitivity of Anderson-Rubin inference
# (Cinelli and Hazlett, 2025), imported from 'iv.sensemakr'.
#
# The quantities below are those of iv.sensemakr. The adjusted critical value
# comes from sensemakr::adjusted_critical_value(); the confidence set is the
# solution of the quadratic inequality in iv.sensemakr's abc_builder(). They
# are restated here because iv.sensemakr does not export them.
# ---------------------------------------------------------------------------

# Coefficient, standard error, and residual df of `term` in an lm fit.
.lm_term <- function(fit, term) {
  cf <- summary(fit)$coefficients
  if (!term %in% rownames(cf)) {
    stop(sprintf("Term '%s' not found in the iv.sensemakr model.", term),
         call. = FALSE)
  }
  list(est = unname(cf[term, 1L]), se = unname(cf[term, 2L]),
       dof = fit$df.residual)
}

# First-stage and reduced-form quantities used by the Anderson-Rubin set.
.iv_parts <- function(fit) {
  fs <- fit$models$fs; rf <- fit$models$rf
  if (is.null(fs) || is.null(rf)) {
    stop("The iv.sensemakr object has no first-stage and reduced-form models.",
         call. = FALSE)
  }
  a <- .lm_term(fs, "z"); b <- .lm_term(rf, "z")
  if (a$dof != b$dof) {
    stop("First-stage and reduced-form degrees of freedom differ.", call. = FALSE)
  }
  list(fs = a$est, fs_se = a$se, rf = b$est, rf_se = b$se,
       rho = stats::cor(stats::resid(fs), stats::resid(rf)), dof = a$dof)
}

# Limits of the bias-adjusted Anderson-Rubin confidence set at partial R2
# values (r2zw, r2y0w) and level alpha. alpha = 1 gives the set of point
# values compatible with the confounding (no sampling uncertainty).
.iv_limits <- function(p, r2zw, r2y0w, alpha) {
  td <- sensemakr::adjusted_critical_value(r2dz.x = r2zw, r2yz.dx = r2y0w,
                                           dof = p$dof, alpha = alpha)
  qa <- p$fs^2 - p$fs_se^2 * td^2
  qb <- 2 * p$rho * p$rf_se * p$fs_se * td^2 - 2 * p$rf * p$fs
  qc <- p$rf^2 - p$rf_se^2 * td^2
  dl <- qb^2 - 4 * qa * qc
  lwr <- upr <- rep(NA_real_, length(td))
  lwr[qa <= 0] <- -Inf
  upr[qa < 0] <- Inf
  ok <- qa > 0 & dl >= 0
  lwr[ok] <- (-qb[ok] - sqrt(dl[ok])) / (2 * qa[ok])
  upr[ok] <- (-qb[ok] + sqrt(dl[ok])) / (2 * qa[ok])
  list(lwr = lwr, upr = upr, t_dagger = td)
}

# The limit of the set nearest to h0, on the side of the estimate. When the
# set becomes unbounded or contains h0 the returned value is on the other
# side of h0 (or infinite), so a sign change marks the crossing.
.iv_near <- function(p, r2zw, r2y0w, alpha, s0) {
  L <- .iv_limits(p, r2zw, r2y0w, alpha)
  if (s0 > 0) L$lwr else L$upr
}

#' Convert an iv.sensemakr result to a confoundsens path
#'
#' Imports an instrumental-variable sensitivity analysis from
#' `iv.sensemakr` (Cinelli and Hazlett, 2025). The confounder in this
#' framework is an omitted variable \eqn{W} that makes the instrument
#' invalid, described by two partial \eqn{R^2} values: \eqn{R^2_{Z \sim W
#' \mid X}} (with the instrument) and \eqn{R^2_{Y_0 \sim W \mid Z, X}} (with
#' the potential outcome under no treatment). Inference uses the
#' Anderson--Rubin confidence set, which stays valid when the instrument is
#' weak.
#'
#' The path follows the equal-strength line \eqn{\lambda = R^2_{Z \sim W \mid
#' X} = R^2_{Y_0 \sim W \mid Z, X}}, the line on which `iv.sensemakr`
#' defines its robustness value. `theta` is the limit, nearest the null
#' \eqn{h_0}, of the set of IV estimates compatible with the confounding
#' (`iv.sensemakr` at `alpha = 1`); it starts at the IV point estimate and
#' reaches \eqn{h_0} at the robustness value for the point estimate.
#' `meta$ci_near` holds the bias-adjusted Anderson--Rubin confidence limit
#' nearest \eqn{h_0} at the level of the `iv.sensemakr` object; it reaches
#' \eqn{h_0} at the robustness value \eqn{RV_{q,\alpha}} that
#' `iv.sensemakr` reports. [robustness_points()] uses it for the
#' significance crossing. An infinite limit means that the adjusted
#' confidence set is unbounded.
#'
#' The robustness values and benchmark bounds of the `iv.sensemakr` object
#' are kept in `meta` (`rv`, `xrv`, `rv_point`, `bounds`). The significance
#' level is fixed when the `iv.sensemakr` object is created; to use another
#' level, call `iv.sensemakr::sensemakr(..., alpha = )` again.
#'
#' @param x An object of class `"iv.sensemakr"`, as returned by
#'   `iv.sensemakr::sensemakr()` applied to an `iv_fit()` result.
#' @param lambda Optional grid of equal partial \eqn{R^2} values in
#'   \eqn{[0, 1)}. By default 2001 points from 0 to three times the larger
#'   robustness value (at least 0.02).
#'
#' @return A `confoundsens` object with `meta$framework = "iv"`.
#' @references Cinelli, C., & Hazlett, C. (2025). An omitted variable bias
#'   framework for sensitivity analysis of instrumental variables.
#'   *Biometrika*, 112(2), asaf004. \doi{10.1093/biomet/asaf004}
#' @seealso [sens_surface()] for the two-parameter plane.
#' @export
#'
#' @examplesIf requireNamespace("iv.sensemakr", quietly = TRUE)
#' data("card", package = "iv.sensemakr")
#' xm <- model.matrix(~ exper + expersq + black + south + smsa + reg661 +
#'   reg662 + reg663 + reg664 + reg665 + reg666 + reg667 + reg668 + smsa66,
#'   data = card)
#' fit <- iv.sensemakr::iv_fit(card$lwage, card$educ, card$nearc4, xm)
#' s <- iv.sensemakr::sensemakr(fit, benchmark_covariates = c("black", "smsa"))
#' path <- from_iv_sensemakr(s)
#' robustness_points(path)
#' sens_report(path)
from_iv_sensemakr <- function(x, lambda = NULL) {
  if (!inherits(x, "iv.sensemakr")) {
    stop("`x` must be an 'iv.sensemakr' object from iv.sensemakr::sensemakr().",
         call. = FALSE)
  }
  if (!requireNamespace("sensemakr", quietly = TRUE)) {
    stop("Package 'sensemakr' is needed for from_iv_sensemakr().", call. = FALSE)
  }
  fit <- x$unadjusted
  est <- fit$estimates$iv
  alpha <- if (!is.null(x$pars$alpha)) x$pars$alpha else est$alpha
  h0 <- if (!is.null(x$pars$h0)) x$pars$h0 else est$h0
  st <- x$sensitivity_stats$iv
  p <- .iv_parts(fit)
  theta0 <- est$estimate
  s0 <- if (theta0 >= h0) 1 else -1

  # Exact equal-strength crossings, by root finding.
  rv_solve <- function(a) {
    g <- function(l) s0 * (.iv_near(p, l, l, a, s0) - h0)
    g0 <- g(0)
    if (!is.finite(g0) || g0 <= 0) return(0)
    hi <- 1 - 1e-9
    gh <- g(hi)
    if (is.finite(gh) && gh > 0) return(1)
    # bracket: g is decreasing and becomes -Inf when the set is unbounded
    lo <- 0; up <- hi
    for (k in 1:200) {
      mid <- (lo + up) / 2
      gm <- g(mid)
      if (is.finite(gm) && gm > 0) lo <- mid else up <- mid
      if (up - lo < 1e-13) break
    }
    if (is.finite(g(up))) {
      stats::uniroot(g, c(lo, up), tol = 1e-13)$root
    } else up
  }
  rv <- if (!is.null(st$rv_qa)) st$rv_qa else rv_solve(alpha)
  rv_point <- rv_solve(1)

  if (is.null(lambda)) {
    top <- min(0.99, max(0.02, 3 * max(rv, rv_point)))
    lambda <- seq(0, top, length.out = 2001L)
  } else {
    .check_lambda_unit(lambda)
    lambda <- sort(unique(lambda))
  }
  theta <- .iv_near(p, lambda, lambda, 1, s0)
  ci <- .iv_limits(p, lambda, lambda, alpha)
  near <- if (s0 > 0) ci$lwr else ci$upr

  bounds <- NULL
  if (!is.null(x$bounds$iv) && nrow(x$bounds$iv)) {
    bi <- x$bounds$iv
    bounds <- data.frame(label = bi$bound_label, a = bi$r2zw.x,
                         b = bi$r2y0w.zx, lwr = bi$lwr, upr = bi$upr,
                         stringsAsFactors = FALSE)
  }
  new_confoundsens(
    lambda = lambda, theta = theta,
    meta = list(
      framework    = "iv",
      source       = "iv.sensemakr",
      treatment    = if (!is.null(fit$pars$d_name)) fit$pars$d_name else NA_character_,
      instrument   = if (!is.null(fit$pars$z_name)) fit$pars$z_name else NA_character_,
      estimate     = theta0,
      ci_unadj     = est$conf.int,
      alpha        = alpha,
      null         = h0,
      dof          = p$dof,
      rv           = rv,
      rv_point     = rv_point,
      xrv          = if (!is.null(st$xrv_qa)) st$xrv_qa else NA_real_,
      iv_parts     = p,
      ci_near      = near,
      ci_lower     = ci$lwr,
      ci_upper     = ci$upr,
      bounds       = bounds,
      lambda_label = "Partial R2 of W with instrument and with Y(0) (equal)",
      theta_label  = "IV estimate bound nearest the null"
    )
  )
}
