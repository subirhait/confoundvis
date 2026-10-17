# ---------------------------------------------------------------------------
# Robustness of inference to replacement (RIR; Frank et al., 2013)
#
# If a proportion lambda of the cases were replaced by cases in which the
# treatment has no effect, the estimate would become (1 - lambda) * estimate.
# The inference is invalidated when this reaches the threshold for
# significance, beta# = sign(estimate) * t_crit * SE, holding the standard
# error fixed as in konfound. The crossing is the proportion of the estimate
# that would have to be due to bias, and RIR = n times that proportion.
# ---------------------------------------------------------------------------

.rir_core <- function(est, beta_thr, n_obs, lambda, alpha, source, treatment,
                      rir = NULL) {
  lambda <- sort(unique(lambda))
  sig <- abs(est) > abs(beta_thr) && est * beta_thr > 0
  prop <- if (sig) 1 - beta_thr / est else NA_real_
  if (is.null(rir)) {
    rir <- if (sig && !is.na(n_obs)) round(n_obs * prop) else NA_real_
  }
  new_confoundsens(
    lambda = lambda, theta = (1 - lambda) * est,
    meta = list(
      framework     = "rir",
      source        = source,
      treatment     = treatment,
      estimate      = est,
      beta_threshold = beta_thr,
      perc_bias     = 100 * prop,
      rir           = rir,
      n_obs         = n_obs,
      alpha         = alpha,
      null          = 0,
      threshold     = beta_thr,
      lambda_label  = "Proportion of cases replaced with null-effect cases",
      theta_label   = "Estimate after replacement"
    )
  )
}

.check_lambda_prop <- function(lambda) {
  if (!is.numeric(lambda) || length(lambda) < 2L || anyNA(lambda) ||
      any(lambda < 0) || any(lambda > 1)) {
    stop("`lambda` must be a numeric vector (length >= 2) with values in [0, 1].",
         call. = FALSE)
  }
}


#' Robustness of inference to replacement (RIR) from a linear model
#'
#' Computes the robustness of inference to replacement of Frank et al. (2013)
#' for a coefficient in a fitted linear model, and the corresponding path.
#' Suppose a proportion \eqn{\lambda}{lambda} of the cases were replaced by
#' cases in which the treatment has no effect. The estimate would then be
#' \eqn{(1 - \lambda)\hat\beta}{(1 - lambda) * beta}. The inference is
#' invalidated when this reaches the threshold for statistical significance
#' \eqn{\beta^\# = \mathrm{sign}(\hat\beta)\, t_{\mathrm{crit}}\, SE}{beta# = sign(beta) * t_crit * SE},
#' with the standard error held fixed. The crossing
#' \eqn{1 - \beta^\# / \hat\beta}{1 - beta# / beta} is the proportion of the
#' estimate that would have to be due to bias, and RIR is the number of
#' cases this represents. The calculation follows `konfound::konfound()` with
#' `index = "RIR"`, using the residual degrees of freedom of the model.
#'
#' @inheritParams sens_path_lm
#' @param alpha Two-sided significance level.
#' @param lambda Numeric vector of replacement proportions in \[0, 1\].
#'
#' @return A list with `estimate`, `beta_threshold`, `perc_bias` (percent of
#'   the estimate that would have to be due to bias), `rir` (number of cases),
#'   `n_obs`, and `path`, a `confoundsens` object whose `lambda_sig` in
#'   [robustness_points()] equals `perc_bias / 100`. For an estimate that is
#'   not significant, `perc_bias` and `rir` are `NA`; `konfound` reports the
#'   corresponding replacement needed to sustain an inference.
#' @references Frank, K. A., Maroulis, S. J., Duong, M. Q., & Kelcey, B. M.
#'   (2013). What would it take to change an inference? Using Rubin's causal
#'   model to interpret the robustness of causal inferences. *Educational
#'   Evaluation and Policy Analysis*, 35(4), 437--460.
#'   \doi{10.3102/0162373713493129}
#' @seealso [itcv_lm()], [from_konfound()]
#' @export
#'
#' @examples
#' fit <- lm(mpg ~ am + wt + hp, data = mtcars)
#' rr <- rir_lm(fit, "wt")
#' rr[c("beta_threshold", "perc_bias", "rir")]
#' sens_report(rr$path)
rir_lm <- function(model, treatment, alpha = 0.05,
                   lambda = seq(0, 1, by = 0.001)) {
  pr <- .lm_parts(model, treatment)
  .check_lambda_prop(lambda)
  if (!is.numeric(alpha) || length(alpha) != 1L || alpha <= 0 || alpha >= 1) {
    stop("`alpha` must be a single number in (0, 1).", call. = FALSE)
  }
  n_obs <- nrow(pr$X)
  beta_thr <- sign(pr$est) * stats::qt(1 - alpha / 2, pr$dof) * pr$se
  path <- .rir_core(pr$est, beta_thr, n_obs, lambda, alpha,
                    source = "lm", treatment = treatment)
  m <- path$meta
  list(estimate = pr$est, beta_threshold = beta_thr, perc_bias = m$perc_bias,
       rir = m$rir, n_obs = n_obs, path = path)
}
