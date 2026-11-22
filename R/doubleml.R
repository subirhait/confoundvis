# ---------------------------------------------------------------------------
# DoubleML: omitted-variable-bias sensitivity for debiased machine learning
# (Chernozhukov, Cinelli, Newey, Sharma, and Syrgkanis, 2024).
#
# The R package 'DoubleML' (1.0.x) has no sensitivity analysis; the Python
# package 'doubleml' does (DoubleML.sensitivity_analysis()). The functions
# below compute the same quantities from an R fit with stored predictions,
# following doubleml 0.11 line by line: the sensitivity elements of
# DoubleMLPLR._sensitivity_element_est() and DoubleMLIRM.
# _sensitivity_element_est(), the bias bound of
# utils._sensitivity._compute_sensitivity_bias(), and the aggregation over
# repetitions of DoubleMLFramework._calc_sensitivity_analysis().
# ---------------------------------------------------------------------------

# Per-repetition sensitivity elements: theta, max_bias, and the moments of
# the scaled score and the influence function of max_bias that give the
# standard error of the bound for any confounding strength.
.dml_elements <- function(x, i_treat) {
  cls <- class(x)
  model <- if ("DoubleMLPLR" %in% cls) "plr" else if ("DoubleMLIRM" %in% cls) "irm" else
    stop("`x` must be a fitted DoubleMLPLR or DoubleMLIRM object.", call. = FALSE)
  if (!is.character(x$score)) {
    stop("Sensitivity analysis needs one of the built-in scores, not a ",
         "user-supplied score function.", call. = FALSE)
  }
  if (model == "irm" && x$score != "ATE") {
    stop("For DoubleMLIRM only score = \"ATE\" is supported. The R ATTE score ",
         "uses fold-wise treatment shares, for which no published sensitivity ",
         "formula is implemented.", call. = FALSE)
  }
  if (inherits(x$data, "DoubleMLClusterData")) {
    stop("Clustered DoubleML data are not supported.", call. = FALSE)
  }
  pr <- x$predictions
  if (is.null(pr) || !length(pr)) {
    stop("`x` has no stored predictions. Refit with ",
         "`x$fit(store_predictions = TRUE)`.", call. = FALSE)
  }
  if (is.null(x$all_coef) || all(is.na(x$all_coef))) {
    stop("`x` has not been fitted.", call. = FALSE)
  }
  dat <- x$data$data_model
  y <- dat[[x$data$y_col]]
  dcol <- x$data$d_cols[i_treat]
  d <- dat[[dcol]]
  n <- length(y)
  n_rep <- x$n_rep
  get_pred <- function(nm, r) {
    p <- pr[[nm]]
    if (is.null(p)) stop(sprintf("Stored prediction '%s' is missing.", nm), call. = FALSE)
    if (length(dim(p)) == 3L) p[, r, i_treat] else if (is.matrix(p)) p[, r] else p
  }
  out <- lapply(seq_len(n_rep), function(r) {
    theta <- x$all_coef[i_treat, r]
    if (model == "plr") {
      m_hat <- get_pred("ml_m", r)
      e2 <- if (x$score == "partialling out") {
        (y - get_pred("ml_l", r) - theta * (d - m_hat))^2
      } else {
        (y - get_pred("ml_g", r) - theta * d)^2
      }
      sigma2 <- mean(e2)
      psi_sigma2 <- e2 - sigma2
      res <- d - m_hat
      nu2 <- 1 / mean(res^2)
      psi_nu2 <- nu2 - res^2 * nu2^2
    } else {
      m_hat <- get_pred("ml_m", r)
      thr <- x$trimming_threshold
      if (identical(x$trimming_rule, "truncate") && !is.null(thr) && thr > 0) {
        m_hat <- pmin(pmax(m_hat, thr), 1 - thr)   # as in R score_elements()
      }
      g0 <- get_pred("ml_g0", r); g1 <- get_pred("ml_g1", r)
      e2 <- (y - d * g1 - (1 - d) * g0)^2
      sigma2 <- mean(e2)
      psi_sigma2 <- e2 - sigma2
      m_alpha <- 1 / m_hat + 1 / (1 - m_hat)       # ATE: weights = 1
      rr <- d / m_hat - (1 - d) / (1 - m_hat)
      nu_el <- 2 * m_alpha - rr^2
      nu2 <- mean(nu_el)
      psi_nu2 <- nu_el - nu2
    }
    if (!is.finite(nu2) || nu2 <= 0 || !is.finite(sigma2)) {
      stop("The estimated Riesz representer variance is not positive; the ",
           "bias bound is undefined for this fit.", call. = FALSE)
    }
    mb <- sqrt(sigma2 * nu2)
    psi_mb <- (sigma2 * psi_nu2 + nu2 * psi_sigma2) / (2 * mb)
    psi <- if (length(dim(x$psi)) == 3L) x$psi[, r, i_treat] else x$psi[, r]
    psi_a <- if (length(dim(x$psi_a)) == 3L) x$psi_a[, r, i_treat] else x$psi_a[, r]
    sp <- psi / mean(psi_a)
    list(theta = theta, max_bias = mb, sigma2 = sigma2, nu2 = nu2,
         A = mean(sp^2), B = mean(sp * psi_mb), C = mean(psi_mb^2), n = n)
  })
  list(model = model, reps = out, n = n, treatment = dcol, score = x$score)
}

# Row medians without a loop over rows.
.row_median <- function(m) {
  k <- ncol(m)
  if (k == 1L) return(m[, 1L])
  s <- matrix(m[order(row(m), m)], ncol = k, byrow = TRUE)
  if (k %% 2L) s[, (k + 1L) / 2L] else (s[, k / 2L] + s[, k / 2L + 1L]) / 2
}

# Bounds for vectors cf_y, cf_d (same length): aggregated over repetitions as
# in doubleml (median theta; SE from the median upper 1.96-limit; median CI).
.dml_bounds <- function(el, cf_y, cf_d, rho, level) {
  cs <- abs(rho) * sqrt(cf_y * cf_d / (1 - cf_d))
  q <- stats::qnorm(level)
  one <- function(side) {
    sgn <- if (side == "lower") -1 else 1
    th <- sapply(el$reps, function(e) e$theta + sgn * cs * e$max_bias)
    sg <- sapply(el$reps, function(e)
      sqrt(pmax(e$A + 2 * sgn * cs * e$B + cs^2 * e$C, 0) / e$n))
    th <- matrix(th, nrow = length(cs)); sg <- matrix(sg, nrow = length(cs))
    theta <- .row_median(th)
    se <- (.row_median(th + 1.96 * sg) - theta) / 1.96
    ci <- .row_median(th + sgn * q * sg)
    list(theta = theta, se = se, ci = ci)
  }
  list(lower = one("lower"), upper = one("upper"))
}

# Robustness value: equal cf_y = cf_d at which the bound (param = "theta")
# or the confidence limit (param = "ci") reaches the null.
.dml_rv <- function(el, rho, level, null, side, param) {
  sgn <- if (side == "lower") 1 else -1
  g <- function(v) sgn * (.dml_bounds(el, v, v, rho, level)[[side]][[param]] - null)
  if (g(0) <= 0) return(0)
  hi <- 0.9999
  if (g(hi) > 0) return(NA_real_)
  stats::uniroot(g, c(0, hi), tol = 1e-12)$root
}

#' Sensitivity analysis for a DoubleML fit
#'
#' Computes the omitted-variable-bias sensitivity analysis of Chernozhukov
#' et al. (2024) for a partially linear (`DoubleMLPLR`) or interactive
#' (`DoubleMLIRM`, average treatment effect) model fitted with the R
#' package `DoubleML`, and returns it as a `confoundsens` path. The R
#' package does not provide this analysis; the Python package `doubleml`
#' does, and `sens_dml()` reproduces its `sensitivity_analysis()` from the
#' stored cross-fitted predictions, following its formulas and its
#' aggregation over repeated sample splits.
#'
#' The confounder is described by two parameters: `cf_y`, the share of the
#' residual variation of the outcome explained by the omitted confounder,
#' and `cf_d`, the relative gain in the variation of the Riesz
#' representer, together with the correlation `rho` between the two bias
#' components (\eqn{|\rho| = 1} is the worst case). The bias bound is
#' \eqn{|\rho| \sqrt{c_y c_d / (1 - c_d)} \, S}, where \eqn{S^2 =
#' \sigma^2 \nu^2} is estimated from the cross-fitted residuals.
#'
#' The path follows the equal-strength line \eqn{\lambda = c_y = c_d}, the
#' line on which `doubleml` defines its robustness values. `theta` is the
#' bound nearest the null (the lower bound when the estimate is above
#' `null`); it reaches `null` at the robustness value RV. `se` is the
#' standard error of that bound, defined so that `theta - qnorm(level) *
#' se` equals the `doubleml` confidence limit. `doubleml` uses a one-sided
#' limit at confidence `level`; this corresponds to a two-sided level
#' \eqn{\alpha = 2(1 - \text{level})}, which is stored as the default
#' significance level of the path (0.10 for `level = 0.95`), so that
#' [robustness_points()] reproduces the robustness value RVa. Supply
#' `alpha` to [robustness_points()] to use another level; with several
#' repeated sample splits the result is then an approximation, because
#' `doubleml` takes medians of the confidence limits.
#'
#' The fit must have been made with `fit(store_predictions = TRUE)`.
#' Clustered data, user-supplied score functions, and the ATTE score are not
#' supported. For `DoubleMLIRM`, propensity scores are truncated with the
#' fit's own `trimming_threshold`, as in its score. The analysis depends on
#' the learners and the sample split, like the estimate itself.
#'
#' @param x A fitted `DoubleMLPLR` or `DoubleMLIRM` object.
#' @param lambda Optional grid of equal `cf_y = cf_d` values in
#'   \eqn{[0, 1)}. By default 2001 points from 0 to 1.5 times the larger
#'   robustness value (at least 0.05, at most 0.99).
#' @param rho Correlation between the bias components, in \eqn{[-1, 1]}.
#' @param level One-sided confidence level, as in `doubleml`.
#' @param null Null value of the effect.
#' @param treatment Treatment column, by name or position, when the model
#'   has several treatments.
#'
#' @return A `confoundsens` object with `meta$framework = "dml"`, carrying
#'   the robustness values `rv` and `rva`, the bias scale `max_bias`
#'   (\eqn{S}), the confidence limits (`ci_lower`, `ci_upper`), and the
#'   elements needed by [sens_surface()].
#' @references Chernozhukov, V., Cinelli, C., Newey, W., Sharma, A., &
#'   Syrgkanis, V. (2024). Long story short: Omitted variable bias in causal
#'   machine learning. NBER Working Paper 30302. \doi{10.3386/w30302}
#'
#'   Bach, P., Chernozhukov, V., Kurz, M. S., & Spindler, M. (2022).
#'   DoubleML: An object-oriented implementation of double machine learning
#'   in Python. *Journal of Machine Learning Research*, 23(53), 1--6.
#' @seealso [sens_surface()] for the plane of `cf_d` and `cf_y`.
#' @export
#'
#' @examplesIf requireNamespace("DoubleML", quietly = TRUE) && requireNamespace("mlr3learners", quietly = TRUE)
#' library(DoubleML)
#' set.seed(1)
#' dat <- make_plr_CCDDHNR2018(n_obs = 500, alpha = 0.5, return_type = "data.table")
#' obj <- DoubleMLData$new(dat, y_col = "y", d_cols = "d")
#' plr <- DoubleMLPLR$new(obj, ml_l = mlr3::lrn("regr.lm"),
#'                        ml_m = mlr3::lrn("regr.lm"), n_folds = 5)
#' plr$fit(store_predictions = TRUE)
#' path <- sens_dml(plr)
#' robustness_points(path)
#' sens_report(path)
sens_dml <- function(x, lambda = NULL, rho = 1, level = 0.95, null = 0,
                     treatment = NULL) {
  if (!inherits(x, "DoubleML")) {
    stop("`x` must be a fitted DoubleML object (DoubleMLPLR or DoubleMLIRM).",
         call. = FALSE)
  }
  if (!is.numeric(rho) || length(rho) != 1L || is.na(rho) || abs(rho) > 1) {
    stop("`rho` must be a single number in [-1, 1].", call. = FALSE)
  }
  if (!is.numeric(level) || length(level) != 1L || is.na(level) ||
      level <= 0 || level >= 1) {
    stop("`level` must be a single number in (0, 1).", call. = FALSE)
  }
  if (!is.numeric(null) || length(null) != 1L || !is.finite(null)) {
    stop("`null` must be a single finite number.", call. = FALSE)
  }
  dcols <- x$data$d_cols
  i_treat <- if (is.null(treatment)) {
    if (length(dcols) > 1L) {
      stop("The model has several treatments; choose one with `treatment`.",
           call. = FALSE)
    }
    1L
  } else if (is.character(treatment)) {
    k <- match(treatment, dcols)
    if (is.na(k)) stop("`treatment` is not a treatment of the model.", call. = FALSE)
    k
  } else as.integer(treatment)

  el <- .dml_elements(x, i_treat)
  b0 <- .dml_bounds(el, 0, 0, rho, level)
  theta0 <- b0$lower$theta
  side <- if (null > theta0) "upper" else "lower"
  rv  <- .dml_rv(el, rho, level, null, side, "theta")
  rva <- .dml_rv(el, rho, level, null, side, "ci")

  if (is.null(lambda)) {
    top <- suppressWarnings(min(0.99, max(0.05, 1.5 * max(c(rv, rva), na.rm = TRUE))))
    if (!is.finite(top)) top <- 0.99
    lambda <- seq(0, top, length.out = 2001L)
  } else {
    .check_lambda_unit(lambda)
    lambda <- sort(unique(lambda))
  }
  bd <- .dml_bounds(el, lambda, lambda, rho, level)
  near <- bd[[side]]
  q <- stats::qnorm(level)
  se_eff <- abs(near$theta - near$ci) / q
  new_confoundsens(
    lambda = lambda, theta = near$theta, se = se_eff,
    meta = list(
      framework     = "dml",
      source        = "DoubleML (formulas of Python doubleml)",
      model         = el$model,
      score         = el$score,
      treatment     = el$treatment,
      estimate      = theta0,
      se0           = b0$lower$se,
      null          = null,
      rho           = rho,
      level         = level,
      default_alpha = round(2 * (1 - level), 12),
      side          = side,
      rv            = rv,
      rva           = rva,
      max_bias      = stats::median(vapply(el$reps, `[[`, numeric(1), "max_bias")),
      n_rep         = length(el$reps),
      dml_elements  = el,
      ci_lower      = bd$lower$ci,
      ci_upper      = bd$upper$ci,
      lambda_label  = "cf_y = cf_d (equal confounding strength)",
      theta_label   = sprintf("Bias-adjusted %s bound", side)
    )
  )
}

# ---- surface: cf_d (a) by cf_y (b) -------------------------------------------
.surf_dml <- function(m, a, b, n, alpha) {
  el <- m$dml_elements; rho <- m$rho; null <- m$null; side <- m$side
  if (is.null(alpha)) alpha <- m$default_alpha
  lev <- 1 - alpha / 2
  sgn <- if (side == "lower") 1 else -1
  top <- suppressWarnings(min(0.99, max(0.05, 2 * max(c(m$rv, m$rva), na.rm = TRUE))))
  if (!is.finite(top)) top <- 0.99
  if (is.null(a)) a <- seq(0, top, length.out = n) else a <- .check_grid(a, 0, 0.999, "a")
  if (is.null(b)) b <- seq(0, top, length.out = n) else b <- .check_grid(b, 0, 0.999, "b")
  ev <- function(aa, bb) .dml_bounds(el, cf_y = bb, cf_d = aa, rho, lev)[[side]]
  gr <- expand.grid(a = a, b = b)
  v <- ev(gr$a, gr$b)
  gr$value <- v$theta
  gr$zone <- ifelse(sgn * (v$theta - null) <= 0, "Estimate at or past zero",
             ifelse(sgn * (v$ci - null) <= 0, "Not significant", NA_character_))
  # Null boundary in closed form for one repetition: cs = |theta - null| / S.
  K <- 2L * n
  if (length(el$reps) == 1L) {
    e <- el$reps[[1L]]
    k2 <- ((e$theta - null) / (abs(rho) * e$max_bias))^2
    af <- seq(1e-6, min(max(a), 0.999), length.out = K)
    bn <- k2 * (1 - af) / af
    nul <- .bnd_piece(af, bn, bn <= max(b), sprintf("Null (bound = %s)", null),
                      "n", "dashed")
  } else {
    gn <- function(aa, bb) sgn * (ev(aa, bb)$theta - null)
    nul <- .trace_zero(gn, seq(0, min(max(a), 0.999), length.out = K),
                       seq(0, min(max(b), 0.999), length.out = K),
                       sprintf("Null (bound = %s)", null), "n", "dashed")
  }
  gs <- function(aa, bb) sgn * (ev(aa, bb)$ci - null)
  sig <- .trace_zero(gs, seq(0, min(max(a), 0.999), length.out = K),
                     seq(0, min(max(b), 0.999), length.out = K),
                     sprintf("Significance (alpha = %s, one-sided level %s)", alpha, lev),
                     "s", "solid")
  list(grid = gr, boundary = rbind(nul, sig),
       meta = list(a_label = "cf_d (gain in Riesz representer variation)",
                   b_label = "cf_y (share of outcome residual variation)",
                   value_label = sprintf("Bias-adjusted %s bound", side),
                   estimate = m$estimate, null = null, rv = m$rv, rva = m$rva,
                   rho = rho, alpha = alpha))
}
