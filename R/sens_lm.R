# ---------------------------------------------------------------------------
# From a fitted linear model to sensitivity quantities
#
# These functions compute, directly from an `lm` fit, the regression-based
# quantities that the plotting layer consumes. They use only the fitted
# design matrix and residual degrees of freedom, so no external sensitivity
# package is required. Results are cross-checked against 'sensemakr' and
# 'konfound' in the test suite when those packages are installed.
# ---------------------------------------------------------------------------

# Internal: pull the pieces we need from an lm fit.
.lm_parts <- function(model, treatment) {
  if (!inherits(model, "lm") || inherits(model, "glm") || inherits(model, "mlm")) {
    stop("`model` must be a single-outcome linear model fitted by stats::lm().",
         call. = FALSE)
  }
  if (!is.character(treatment) || length(treatment) != 1L) {
    stop("`treatment` must be a single character string naming a coefficient.",
         call. = FALSE)
  }
  X <- stats::model.matrix(model)
  if (!treatment %in% colnames(X)) {
    stop(sprintf(
      "`treatment` '%s' is not a column of the model matrix. Available: %s.",
      treatment, paste(colnames(X), collapse = ", ")), call. = FALSE)
  }
  cf <- summary(model)$coefficients
  if (!treatment %in% rownames(cf)) {
    stop(sprintf("Coefficient for '%s' is not estimable (aliased?).", treatment),
         call. = FALSE)
  }
  y <- stats::model.response(stats::model.frame(model))
  w <- stats::weights(model)
  if (!is.null(w)) {
    sw <- sqrt(w)
    X  <- X * sw
    y  <- y * sw
  }
  list(
    X        = X,
    y        = as.numeric(y),
    assign   = attr(stats::model.matrix(model), "assign"),
    labels   = attr(stats::terms(model), "term.labels"),
    est      = unname(cf[treatment, 1L]),
    se       = unname(cf[treatment, 2L]),
    t        = unname(cf[treatment, 3L]),
    dof      = stats::df.residual(model),
    treatment = treatment
  )
}

# Internal: residual sum of squares of y on X (X may have zero columns).
.rss <- function(y, X) {
  if (is.null(X) || NCOL(X) == 0L) return(sum(y^2))
  fit <- stats::lm.fit(as.matrix(X), y)
  sum(fit$residuals^2)
}

# Internal: coefficient of column `j` of X in the regression of y on X.
.coef_col <- function(y, X, j) {
  fit <- stats::lm.fit(as.matrix(X), y)
  unname(fit$coefficients[j])
}

# Internal: equal-strength (or fixed-ratio) partial R^2 path.
# Bias and adjusted SE follow Cinelli and Hazlett (2020), Eqs. 4-5.
.pr2_path <- function(est, se, dof, lambda, ratio) {
  r2dz <- lambda
  r2yz <- pmin(ratio * lambda, 1)
  bias <- se * sqrt(dof) * sqrt(r2yz * r2dz / (1 - r2dz))
  adj  <- sign(est) * (abs(est) - bias)            # reduces |estimate|
  adj_se <- se * sqrt((1 - r2yz) / (1 - r2dz)) * sqrt(dof / (dof - 1))
  list(theta = adj, se = adj_se, t = adj / adj_se)
}

# Internal: robustness value (Cinelli and Hazlett 2020, Eq. 18) for q = 1.
.rv <- function(t, dof) {
  f <- abs(t) / sqrt(dof)
  0.5 * (sqrt(f^4 + 4 * f^2) - f^2)
}


#' Sensitivity path from a fitted linear model (partial R-squared framework)
#'
#' Computes the bias-adjusted treatment coefficient as a function of the
#' strength of a hypothetical omitted confounder, using the omitted-variable
#' bias formulas of Cinelli and Hazlett (2020). Confounder strength is indexed
#' by \eqn{\lambda}, the partial \eqn{R^2} of the confounder with the
#' treatment (\eqn{R^2_{D \sim Z | X}}); the partial \eqn{R^2} with the outcome
#' is \eqn{R^2_{Y \sim Z | D, X} = \min(\code{ratio} \times \lambda, 1)}. The
#' default `ratio = 1` gives the equal-strength path along which the
#' robustness value is defined.
#'
#' The confounder is assumed to act in the direction that shrinks the estimate
#' toward zero (the conventional worst case for a reported effect).
#'
#' @param model A fitted `lm` object.
#' @param treatment Character string naming the treatment coefficient (a
#'   column of the model matrix).
#' @param lambda Numeric vector of confounder strengths in \[0, 1). Default
#'   is a grid from 0 to 0.5.
#' @param ratio Positive scalar; ratio of the outcome-side to the
#'   treatment-side partial \eqn{R^2}.
#'
#' @return A `confoundsens` object with `theta` (adjusted estimate), `se`
#'   (adjusted standard error), and `t` (adjusted t statistic), plus `meta`
#'   recording the original estimate, degrees of freedom, and robustness value.
#'
#' @references Cinelli, C., & Hazlett, C. (2020). Making sense of sensitivity:
#'   Extending omitted variable bias. *Journal of the Royal Statistical
#'   Society: Series B*, 82(1), 39--67. \doi{10.1111/rssb.12348}
#'
#' @seealso [from_sensemakr()], [itcv_lm()], [covariate_impacts()]
#' @export
#'
#' @examples
#' fit <- lm(mpg ~ am + wt + hp, data = mtcars)
#' path <- sens_path_lm(fit, treatment = "am")
#' path
#' plot_robustness_curve(path, points = FALSE)
sens_path_lm <- function(model, treatment,
                         lambda = seq(0, 0.5, by = 0.01),
                         ratio  = 1) {
  pr <- .lm_parts(model, treatment)
  .check_lambda_unit(lambda)
  if (!is.numeric(ratio) || length(ratio) != 1L || is.na(ratio) || ratio <= 0) {
    stop("`ratio` must be a single positive number.", call. = FALSE)
  }
  .path_from_stats(pr$est, pr$se, pr$dof, lambda, ratio,
                   treatment = treatment, source = "lm")
}

.check_lambda_unit <- function(lambda) {
  if (!is.numeric(lambda) || length(lambda) < 2L || anyNA(lambda) ||
      any(lambda < 0) || any(lambda >= 1)) {
    stop("`lambda` must be a numeric vector (length >= 2) with values in [0, 1).",
         call. = FALSE)
  }
}

.path_from_stats <- function(est, se, dof, lambda, ratio, treatment, source) {
  lambda <- sort(unique(lambda))
  p <- .pr2_path(est, se, dof, lambda, ratio)
  new_confoundsens(
    lambda = lambda, theta = p$theta, se = p$se, t = p$t,
    meta = list(
      framework    = "partial_r2",
      source       = source,
      treatment    = treatment,
      estimate     = est,
      se0          = se,
      dof          = dof,
      ratio        = ratio,
      rv           = .rv(est / se, dof),
      null         = 0,
      lambda_label = if (ratio == 1)
        "Confounder partial R2 with treatment and outcome"
      else sprintf("Confounder partial R2 with treatment (outcome side = %g x)", ratio),
      theta_label  = "Bias-adjusted estimate"
    )
  )
}


#' Impact threshold for a confounding variable from a fitted linear model
#'
#' Computes Frank's (2000) impact threshold for a confounding variable (ITCV)
#' for one coefficient of a linear model, together with the ITCV sensitivity
#' path: the partial correlation between treatment and outcome after
#' conditioning on a confounder whose impact \eqn{k = r_{DU} r_{YU}} is split
#' equally between the two correlations (\eqn{r_{DU} = r_{YU} = \sqrt{k}}).
#'
#' Along this path the adjusted correlation is \eqn{(r - k)/(1 - k)}, which
#' equals the critical correlation \eqn{r^\#} exactly at
#' \eqn{k = \mathrm{ITCV} = (|r| - |r^\#|)/(1 - |r^\#|)}.
#'
#' @param model A fitted `lm` object.
#' @param treatment Character string naming the treatment coefficient.
#' @param alpha Two-sided significance level used for the critical value.
#' @param lambda Numeric vector of impacts in \[0, 1) at which to evaluate
#'   the path.
#'
#' @return A list with elements `r` (observed partial correlation), `r_crit`
#'   (critical partial correlation), `itcv` (the ITCV for a confounder's
#'   correlations *conditional on* the model covariates), `dof`, `path` (a
#'   `confoundsens` object on the correlation scale), `r2xz` and `r2yz` (the
#'   \eqn{R^2} of the covariates with treatment and outcome), and
#'   `itcv_unconditional` (the ITCV expressed in unconditional correlations,
#'   `itcv * sqrt((1 - r2xz) * (1 - r2yz))`).
#'
#' @details The conditional ITCV equals `itcvGz` and the unconditional ITCV
#'   equals `itcv` in the raw output of `konfound::konfound()` (version 1.0).
#'
#' @references Frank, K. A. (2000). Impact of a confounding variable on a
#'   regression coefficient. *Sociological Methods & Research*, 29(2),
#'   147--194. \doi{10.1177/0049124100029002001}
#'
#'   Frank, K. A., Maroulis, S. J., Duong, M. Q., & Kelcey, B. M. (2013).
#'   What would it take to change an inference? Using Rubin's causal model to
#'   interpret the robustness of causal inferences. *Educational Evaluation
#'   and Policy Analysis*, 35(4), 437--460. \doi{10.3102/0162373713493129}
#' @seealso [from_konfound()], [sens_path_lm()]
#' @export
#'
#' @examples
#' fit <- lm(mpg ~ am + wt + hp, data = mtcars)
#' it <- itcv_lm(fit, "am")
#' it$itcv
#' plot_robustness_curve(it$path, points = FALSE)
itcv_lm <- function(model, treatment, alpha = 0.05,
                    lambda = seq(0, 0.5, by = 0.01)) {
  pr <- .lm_parts(model, treatment)
  .check_lambda_unit(lambda)
  if (!is.numeric(alpha) || length(alpha) != 1L || alpha <= 0 || alpha >= 1) {
    stop("`alpha` must be a single number in (0, 1).", call. = FALSE)
  }
  out <- .itcv_core(r = pr$t / sqrt(pr$t^2 + pr$dof), dof = pr$dof,
                    alpha = alpha, lambda = lambda, treatment = treatment,
                    source = "lm")
  # Unconditional ITCV (Frank et al. 2013): rescale by the covariates' R^2
  # with treatment and outcome.
  jD <- match(treatment, colnames(pr$X))
  D  <- pr$X[, jD]
  W  <- pr$X[, -jD, drop = FALSE]
  has_int <- any(pr$assign == 0L)
  tss <- function(v) if (has_int) sum((v - mean(v))^2) else sum(v^2)
  r2xz <- if (ncol(W) > (if (has_int) 1L else 0L)) 1 - .rss(D, W) / tss(D) else 0
  r2yz <- if (ncol(W) > (if (has_int) 1L else 0L)) 1 - .rss(pr$y, W) / tss(pr$y) else 0
  out$r2xz <- r2xz
  out$r2yz <- r2yz
  out$itcv_unconditional <- out$itcv * sqrt((1 - r2xz) * (1 - r2yz))
  out$path$meta$itcv_unconditional <- out$itcv_unconditional
  out
}

.itcv_core <- function(r, dof, alpha, lambda, treatment, source,
                       r_crit = NULL) {
  if (is.null(r_crit)) {
    tc <- stats::qt(1 - alpha / 2, dof)
    r_crit <- tc / sqrt(tc^2 + dof)
  }
  r_crit <- abs(r_crit)
  itcv <- (abs(r) - r_crit) / (1 - r_crit)
  lambda <- sort(unique(lambda))
  r_adj <- sign(r) * (abs(r) - lambda) / (1 - lambda)
  t_adj <- if (!is.null(dof) && is.finite(dof))
    r_adj * sqrt(dof - 1) / sqrt(pmax(1 - r_adj^2, .Machine$double.eps))
  else NULL
  path <- new_confoundsens(
    lambda = lambda, theta = r_adj, t = t_adj,
    meta = list(
      framework    = "itcv",
      source       = source,
      treatment    = treatment,
      estimate     = r,
      r_crit       = r_crit,
      itcv         = itcv,
      dof          = dof,
      alpha        = alpha,
      null         = 0,
      threshold    = sign(r) * r_crit,
      lambda_label = "Confounder impact k = r(D,U) x r(Y,U)",
      theta_label  = "Adjusted partial correlation"
    )
  )
  list(r = r, r_crit = r_crit, itcv = itcv, dof = dof, path = path)
}


#' Benchmark observed covariates against a sensitivity threshold
#'
#' For each covariate term in a fitted linear model (other than the treatment
#' and the intercept), computes how strongly it is related to the treatment
#' and to the outcome. These are the observed benchmarks that analysts compare
#' with a sensitivity threshold ("would an omitted variable as strong as this
#' covariate overturn the result?").
#'
#' Two metrics are returned:
#' \describe{
#'   \item{`impact` (ITCV metric)}{\eqn{r_{DZ|W} \times r_{YZ|W}}, the product
#'     of partial correlations of the covariate with the treatment and with
#'     the outcome, each conditioning on the remaining covariates \eqn{W}
#'     (as in Frank, 2000). Compared with the ITCV. Only positive values (same
#'     sign as the estimate) push the estimate toward zero; the returned
#'     `impact` is signed so that positive means "reduces the estimate".
#'     Defined only for single-column terms.}
#'   \item{`bias_index` (partial \eqn{R^2} metric)}{
#'     \eqn{\sqrt{R^2_{Y \sim Z|D,W} R^2_{D \sim Z|W} / (1 - R^2_{D \sim Z|W})}}.
#'     An omitted confounder with this pair of partial \eqn{R^2} values would
#'     move the estimate by `bias_index` \eqn{\times\, se \sqrt{df}}, so it
#'     would bring the point estimate to zero when `bias_index` exceeds
#'     \eqn{|t|/\sqrt{df}}. Defined for multi-column terms (factors) as well.}
#' }
#'
#' The returned data frame carries an attribute `threshold` for the selected
#' `metric`, so it can be passed straight to [plot_sensitivity_love()]. Its
#' `r_yu`, `r_du`, and `label` columns can be passed as `benchmarks` to
#' [plot_sensitivity_contour()].
#'
#' Covariate benchmarks describe *observed* variables. Whether an unobserved
#' confounder could be as strong as a given covariate is a substantive
#' judgment that no benchmark can settle.
#'
#' @param model A fitted `lm` object.
#' @param treatment Character string naming the treatment coefficient.
#' @param metric Which metric populates the `impact` column used by
#'   [plot_sensitivity_love()]: `"itcv"` (default) or `"partial_r2"`.
#' @param alpha Two-sided significance level used for the ITCV threshold.
#'
#' @return A data frame with one row per covariate term and columns
#'   `covariate`, `term_df`, `r_du`, `r_yu`, `itcv_impact`, `r2dz.x`,
#'   `r2yz.dx`, `bias_index`, `impact`, and `label`, and attributes
#'   `threshold`, `metric`, `itcv`, and `f_threshold`.
#' @export
#'
#' @examples
#' fit <- lm(mpg ~ am + wt + hp + qsec, data = mtcars)
#' imp <- covariate_impacts(fit, "am")
#' imp
#' plot_sensitivity_love(imp)
covariate_impacts <- function(model, treatment,
                              metric = c("itcv", "partial_r2"),
                              alpha  = 0.05) {
  metric <- match.arg(metric)
  pr <- .lm_parts(model, treatment)
  X  <- pr$X
  y  <- pr$y
  jD <- match(treatment, colnames(X))
  D  <- X[, jD]
  W  <- X[, -jD, drop = FALSE]
  asg_W <- pr$assign[-jD]
  trt_term <- pr$assign[jD]

  term_ids <- setdiff(unique(asg_W), c(0L, trt_term))
  if (length(term_ids) == 0L) {
    stop("The model has no covariate terms besides the treatment.", call. = FALSE)
  }

  rss_d_full <- .rss(D, W)
  rss_y_full <- .rss(y, cbind(D, W))
  rss_yw_full <- .rss(y, W)

  rows <- lapply(term_ids, function(k) {
    cols  <- which(asg_W == k)
    W_mk  <- W[, -cols, drop = FALSE]
    r2dz  <- (.rss(D, W_mk) - rss_d_full) / .rss(D, W_mk)
    r2yz  <- (.rss(y, cbind(D, W_mk)) - rss_y_full) / .rss(y, cbind(D, W_mk))
    r2yzw <- (.rss(y, W_mk) - rss_yw_full) / .rss(y, W_mk)
    if (length(cols) == 1L) {
      r_du <- sign(.coef_col(D, W, cols)) * sqrt(max(r2dz, 0))
      r_yu <- sign(.coef_col(y, W, cols)) * sqrt(max(r2yzw, 0))
    } else {
      r_du <- NA_real_
      r_yu <- NA_real_
    }
    data.frame(
      covariate  = pr$labels[k],
      term_df    = length(cols),
      r_du       = r_du,
      r_yu       = r_yu,
      itcv_impact = sign(pr$est) * r_du * r_yu,
      r2dz.x     = r2dz,
      r2yz.dx    = r2yz,
      bias_index = sqrt(r2yz * r2dz / (1 - r2dz)),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)

  itcv <- .itcv_core(r = pr$t / sqrt(pr$t^2 + pr$dof), dof = pr$dof,
                     alpha = alpha, lambda = c(0, 0.5),
                     treatment = treatment, source = "lm")$itcv
  f_thr <- abs(pr$t) / sqrt(pr$dof)

  out$impact <- if (metric == "itcv") out$itcv_impact else out$bias_index
  out$label  <- out$covariate
  rownames(out) <- NULL
  if (metric == "itcv" && anyNA(out$impact)) {
    message("Multi-column terms have no ITCV impact (NA); ",
            "use metric = \"partial_r2\" to include them.")
  }
  attr(out, "threshold")   <- if (metric == "itcv") itcv else f_thr
  attr(out, "metric")      <- metric
  attr(out, "itcv")        <- itcv
  attr(out, "f_threshold") <- f_thr
  out
}
