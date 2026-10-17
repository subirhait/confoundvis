#' Create a confoundsens object
#'
#' Constructs a `confoundsens` object — a lightweight container for storing
#' sensitivity paths (e.g., \eqn{\theta(\lambda)}) and optional uncertainty or
#' stratification information used by `confoundvis` plotting functions.
#'
#' @param lambda Numeric vector of sensitivity strength values (e.g., ITCV
#'   lambda or delta). Should be nondecreasing; a warning is issued otherwise.
#' @param theta Numeric vector of effect estimates along the sensitivity path.
#'   Must be the same length as `lambda`.
#' @param level Optional character (or coercible) vector of subgroup (stratum) labels
#'   (e.g., `"female"`, `"male"`). Must be the same length as `lambda`.
#' @param se Optional numeric vector of standard errors for `theta`. Must be
#'   the same length as `lambda`.
#' @param t Optional numeric vector of test statistics along the path. Must be
#'   the same length as `lambda`.
#'
#' @param meta Optional named list of metadata describing the path, for
#'   example `framework` (`"partial_r2"`, `"itcv"`, `"evalue"`, ...),
#'   `lambda_label`, `theta_label`, `null` (the value of `theta` that
#'   represents no effect), and `estimate`. Used by printing, plotting, and
#'   [sens_report()]; ignored otherwise.
#'
#' @return A `confoundsens` object (a list with class `"confoundsens"`).
#' @export
#'
#' @examples
#' x <- new_confoundsens(
#'   lambda = seq(0, 0.2, length.out = 10),
#'   theta  = seq(1, 0.6, length.out = 10),
#'   se     = rep(0.1, 10)
#' )
#' x
new_confoundsens <- function(lambda,
                             theta,
                             level = NULL,
                             se    = NULL,
                             t     = NULL,
                             meta  = list()) {
  if (!is.numeric(lambda)) stop("`lambda` must be a numeric vector.", call. = FALSE)
  if (!is.numeric(theta))  stop("`theta` must be a numeric vector.",  call. = FALSE)

  x <- list(
    lambda = lambda,
    theta  = theta,
    level  = level,
    se     = se,
    t      = t,
    meta   = meta
  )
  class(x) <- "confoundsens"

  # Capture coerced result (validate modifies x$level via as.character)
  x <- validate_confoundsens(x)
  x
}
