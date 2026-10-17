#' Convert a data frame to a confoundsens object
#'
#' Converts a data frame containing sensitivity analysis results into a
#' `confoundsens` object suitable for use with `confoundvis` plotting functions.
#'
#' @param data A data.frame containing sensitivity analysis results, a fitted
#'   `lm` model, or a `sensemakr` object.
#' @param lambda Character string; name of the column containing lambda values.
#'   Default `"lambda"`.
#' @param theta Character string; name of the column containing theta values.
#'   Default `"theta"`.
#' @param level Optional character string; name of the column containing subgroup
#'   (stratum) labels (e.g., `"female"` / `"male"`).
#' @param se Optional character string; name of the column containing standard
#'   errors for `theta`.
#' @param t Optional character string; name of the column containing test
#'   statistics.
#'
#' @details `as_confoundsens()` is generic. The data-frame method (the 0.1.0
#' behaviour) maps columns to fields. Methods for fitted `lm` models and
#' `sensemakr` objects compute the partial \eqn{R^2} path via
#' [sens_path_lm()] and [from_sensemakr()]; see also [from_konfound()] and
#' [from_evalue()].
#'
#' @param ... Passed to methods. For `lm` objects, `treatment` (required),
#'   `lambda`, and `ratio`; for `sensemakr` objects, `lambda` and `ratio`.
#' @return A `confoundsens` object.
#'
#' @examples
#' df <- data.frame(
#'   lambda = seq(0, 0.2, length.out = 10),
#'   theta  = seq(1, 0.5, length.out = 10),
#'   se     = rep(0.1, 10),
#'   level  = rep(c("female", "male"), length.out = 10)
#' )
#' x <- as_confoundsens(df)
#' x
#'
#' # Directly from a fitted model
#' fit <- lm(mpg ~ am + wt + hp, data = mtcars)
#' as_confoundsens(fit, treatment = "am")
#' @export
as_confoundsens <- function(data, ...) {
  UseMethod("as_confoundsens")
}

#' @rdname as_confoundsens
#' @export
as_confoundsens.default <- function(data, ...) {
  stop("`data` must be a data.frame, an lm fit, or a sensemakr object; got ",
       paste(class(data), collapse = "/"), ".", call. = FALSE)
}

#' @rdname as_confoundsens
#' @export
as_confoundsens.lm <- function(data, ...) {
  sens_path_lm(data, ...)
}

#' @rdname as_confoundsens
#' @export
as_confoundsens.sensemakr <- function(data, ...) {
  from_sensemakr(data, ...)
}

#' @rdname as_confoundsens
#' @export
as_confoundsens.data.frame <- function(data,
                            lambda = "lambda",
                            theta  = "theta",
                            level  = NULL,
                            se     = NULL,
                            t      = NULL,
                            ...) {

  .check_col <- function(nm, arg) {
    if (!is.character(nm) || length(nm) != 1L) {
      stop(sprintf("`%s` must be a single character string.", arg), call. = FALSE)
    }
    if (!nm %in% names(data)) {
      stop(sprintf("Column '%s' not found in `data`.", nm), call. = FALSE)
    }
  }

  .check_col(lambda, "lambda")
  .check_col(theta,  "theta")

  lam <- data[[lambda]]
  th  <- data[[theta]]
  if (!is.numeric(lam)) stop("`lambda` column must be numeric.", call. = FALSE)
  if (!is.numeric(th))  stop("`theta` column must be numeric.",  call. = FALSE)

  lvl <- NULL
  if (!is.null(level)) {
    .check_col(level, "level")
    lvl <- as.character(data[[level]])
  }

  se_ <- NULL
  if (!is.null(se)) {
    .check_col(se, "se")
    se_ <- as.numeric(data[[se]])
  }

  t_ <- NULL
  if (!is.null(t)) {
    .check_col(t, "t")
    t_ <- as.numeric(data[[t]])
  }

  new_confoundsens(lambda = lam, theta = th, level = lvl, se = se_, t = t_)
}
