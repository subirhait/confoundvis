# ---------------------------------------------------------------------------
# Adapters: results from other sensitivity packages -> confoundsens
#
# These functions read only documented, stable fields of the source objects
# and never call the source packages, so those packages stay in Suggests.
# ---------------------------------------------------------------------------

#' Convert a sensemakr result to a confoundsens path
#'
#' Builds the partial \eqn{R^2} sensitivity path from the summary statistics
#' stored in a `sensemakr` object (estimate, standard error, and residual
#' degrees of freedom), using the same bias formulas as [sens_path_lm()].
#' Benchmark bounds stored in the object are attached as `meta$bounds`.
#'
#' @param x An object of class `"sensemakr"` (from `sensemakr::sensemakr()`).
#' @inheritParams sens_path_lm
#'
#' @return A `confoundsens` object.
#' @seealso [sens_path_lm()]
#' @export
#'
#' @examplesIf requireNamespace("sensemakr", quietly = TRUE)
#' data("darfur", package = "sensemakr")
#' fit <- lm(peacefactor ~ directlyharmed + age + farmer_dar + herder_dar +
#'             pastvoted + hhsize_darfur + female + village, data = darfur)
#' s <- sensemakr::sensemakr(fit, treatment = "directlyharmed",
#'                           benchmark_covariates = "female", kd = 1:3)
#' path <- from_sensemakr(s)
#' sens_report(path)
from_sensemakr <- function(x, lambda = seq(0, 0.5, by = 0.01), ratio = 1) {
  if (!inherits(x, "sensemakr")) {
    stop("`x` must be an object of class 'sensemakr'.", call. = FALSE)
  }
  st <- x$sensitivity_stats
  need <- c("treatment", "estimate", "se", "dof")
  if (!is.data.frame(st) || !all(need %in% names(st))) {
    stop("`x$sensitivity_stats` does not have the expected fields: ",
         paste(need, collapse = ", "), ".", call. = FALSE)
  }
  if (nrow(st) != 1L) {
    stop("Only single-treatment sensemakr objects are supported.", call. = FALSE)
  }
  .check_lambda_unit(lambda)
  out <- .path_from_stats(st$estimate, st$se, st$dof, lambda, ratio,
                          treatment = st$treatment, source = "sensemakr")
  if (!is.null(x$bounds)) out$meta$bounds <- x$bounds
  out
}


#' Convert konfound output to a confoundsens path
#'
#' Builds a path from the output of `konfound::konfound()` or
#' `konfound::pkonfound()` called with `to_return = "raw_output"`.
#'
#' With `index = "IT"` (default), the impact threshold (ITCV) path is built
#' from the observed and critical correlations (`obs_r`, `critical_r`). With
#' `index = "RIR"`, the robustness-of-inference-to-replacement path is built
#' from `beta_threshold` and `perc_bias_to_change`, and konfound's
#' `RIR_primary` is carried in `meta$rir`; see [rir_lm()] for the path.
#'
#' @param x A list returned by `konfound()`/`pkonfound()` with
#'   `to_return = "raw_output"`.
#' @param dof Optional residual degrees of freedom (ITCV only). If supplied,
#'   the path also carries adjusted t statistics.
#' @param lambda Numeric vector of impacts in \[0, 1) for `index = "IT"`, or
#'   of replacement proportions in \[0, 1\] for `index = "RIR"` (default
#'   `seq(0, 1, by = 0.001)` when not supplied).
#' @param treatment Optional label for the treatment.
#' @param index `"IT"` (default) or `"RIR"`.
#'
#' @return A `confoundsens` object. For `index = "IT"` it is on the
#'   correlation scale and `meta$itcv` holds the conditional ITCV (konfound's
#'   `itcvGz`); for `index = "RIR"` it is on the scale of the estimate and
#'   `meta$perc_bias` and `meta$rir` hold konfound's percent bias and RIR.
#' @seealso [itcv_lm()], [rir_lm()]
#' @export
#'
#' @examplesIf requireNamespace("konfound", quietly = TRUE)
#' k <- konfound::pkonfound(est_eff = 2, std_err = 0.4, n_obs = 100,
#'                          n_covariates = 3, index = "IT",
#'                          to_return = "raw_output")
#' path <- from_konfound(k)
#' path$meta$itcv
#' rir <- from_konfound(k, index = "RIR")
#' rir$meta$rir
from_konfound <- function(x, dof = NULL, lambda = seq(0, 0.5, by = 0.01),
                          treatment = NA_character_, index = c("IT", "RIR")) {
  index <- match.arg(index)
  if (index == "RIR") {
    need <- c("beta_threshold", "perc_bias_to_change")
    if (!is.list(x) || !all(need %in% names(x)) ||
        is.null(x$beta_threshold) || is.null(x$perc_bias_to_change)) {
      stop("`x` must be konfound raw output containing `beta_threshold` and ",
           "`perc_bias_to_change` (call konfound()/pkonfound() with ",
           "to_return = \"raw_output\").", call. = FALSE)
    }
    if (missing(lambda)) lambda <- seq(0, 1, by = 0.001)
    .check_lambda_prop(lambda)
    bt <- as.numeric(x$beta_threshold)
    pb <- as.numeric(x$perc_bias_to_change)
    if (!is.finite(pb) || pb <= 0 || pb >= 100) {
      stop("konfound reports no RIR to invalidate the inference ",
           "(perc_bias_to_change = ", pb, "); the estimate may not be significant.",
           call. = FALSE)
    }
    est <- bt / (1 - pb / 100)
    rir <- if (is.null(x$RIR_primary)) NA_real_ else as.numeric(x$RIR_primary)
    return(.rir_core(est, bt, n_obs = NA_real_, lambda = lambda, alpha = NA_real_,
                     source = "konfound", treatment = treatment, rir = rir))
  }
  if (!is.list(x) || !all(c("obs_r", "critical_r") %in% names(x)) ||
      is.null(x$obs_r) || is.null(x$critical_r)) {
    stop("`x` must be konfound raw output containing `obs_r` and `critical_r` ",
         "(call konfound()/pkonfound() with to_return = \"raw_output\").",
         call. = FALSE)
  }
  .check_lambda_unit(lambda)
  .itcv_core(r = as.numeric(x$obs_r), dof = if (is.null(dof)) NA_real_ else dof,
             alpha = NA_real_, lambda = lambda, treatment = treatment,
             source = "konfound", r_crit = as.numeric(x$critical_r))$path
}


#' Convert an E-value result to a confoundsens path
#'
#' Builds a bias-adjusted risk-ratio path from an `evalue` object (as returned
#' by `EValue::evalues.RR()`, `evalues.OR()`, `evalues.HR()`, `evalues.OLS()`,
#' and related functions). Confounder strength \eqn{\lambda} is the common
#' risk ratio relating the confounder to both treatment and outcome; the
#' corresponding bias factor is \eqn{B = \lambda^2/(2\lambda - 1)}
#' (VanderWeele and Ding, 2017). The path is on the log risk-ratio scale, so
#' it crosses the null (0) exactly at the E-value.
#'
#' @param x An object of class `"evalue"`.
#' @param lambda Numeric vector of confounder risk ratios (>= 1).
#'
#' @return A `confoundsens` object with `theta` = log of the bias-adjusted
#'   point risk ratio and, when the confidence interval is available,
#'   `se` on the log scale.
#' @references VanderWeele, T. J., & Ding, P. (2017). Sensitivity analysis in
#'   observational research: Introducing the E-value. *Annals of Internal
#'   Medicine*, 167(4), 268--274. \doi{10.7326/M16-2607}
#' @export
#'
#' @examplesIf requireNamespace("EValue", quietly = TRUE)
#' e <- EValue::evalues.RR(est = 1.8, lo = 1.3, hi = 2.5)
#' path <- from_evalue(e)
#' robustness_points(path)
from_evalue <- function(x, lambda = seq(1, 5, by = 0.02)) {
  if (!inherits(x, "evalue") || !is.matrix(x) || !"RR" %in% rownames(x)) {
    stop("`x` must be an 'evalue' object with an 'RR' row.", call. = FALSE)
  }
  if (!is.numeric(lambda) || length(lambda) < 2L || anyNA(lambda) ||
      any(lambda < 1)) {
    stop("`lambda` must be a numeric vector (length >= 2) with values >= 1.",
         call. = FALSE)
  }
  rr <- x["RR", ]
  pt <- unname(rr[["point"]])
  lo <- unname(rr[["lower"]]); hi <- unname(rr[["upper"]])
  lambda <- sort(unique(lambda))
  B <- lambda^2 / (2 * lambda - 1)
  log_adj <- if (pt >= 1) log(pt) - log(B) else log(pt) + log(B)
  se <- if (!is.na(lo) && !is.na(hi) && lo > 0 && hi > 0) {
    rep((log(hi) - log(lo)) / (2 * stats::qnorm(0.975)), length(lambda))
  } else NULL
  ev_pt <- unname(x["E-values", "point"])
  new_confoundsens(
    lambda = lambda, theta = log_adj, se = se,
    meta = list(
      framework    = "evalue",
      source       = "EValue",
      estimate     = pt,
      evalue       = ev_pt,
      null         = 0,
      lambda_label = "Confounder risk ratio with treatment and outcome",
      theta_label  = "log bias-adjusted risk ratio"
    )
  )
}
