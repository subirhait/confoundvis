#' Print method for confoundsens objects
#'
#' Prints a concise summary of the key fields in a `confoundsens` object.
#'
#' @param x A `confoundsens` object.
#' @param ... Unused; included for S3 method consistency.
#'
#' @return `x`, invisibly.
#' @export
#'
#' @examples
#' x <- new_confoundsens(
#'   lambda = seq(0, 0.2, length.out = 5),
#'   theta  = seq(1, 0.8, length.out = 5)
#' )
#' print(x)
print.confoundsens <- function(x, ...) {
  x <- validate_confoundsens(x)

  n         <- length(x$lambda)
  has_level <- !is.null(x$level)
  has_se    <- !is.null(x$se)
  has_t     <- !is.null(x$t)

  cat("<confoundsens>\n")
  cat("  n      :", n, "\n")
  cat("  lambda :", sprintf("[%g, %g]",
      min(x$lambda, na.rm = TRUE), max(x$lambda, na.rm = TRUE)), "\n")
  cat("  theta  :", sprintf("[%g, %g]",
      min(x$theta,  na.rm = TRUE), max(x$theta,  na.rm = TRUE)), "\n")

  if (has_level) {
    lvls <- unique(x$level)
    cat("  level  :", paste0(length(lvls), " level(s): ",
        paste(lvls, collapse = ", ")), "\n")
  } else {
    cat("  level  : <none>\n")
  }

  cat("  se     :", if (has_se) "yes" else "no", "\n")
  cat("  t      :", if (has_t)  "yes" else "no", "\n")
  if (length(x$meta) && !is.null(x$meta$framework)) {
    cat("  frame  :", x$meta$framework,
        if (!is.null(x$meta$source)) paste0("(from ", x$meta$source, ")") else "",
        "\n")
  }

  invisible(x)
}
