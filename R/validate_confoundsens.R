#' Validate a confoundsens object
#'
#' Internal validator used by constructors and methods. Checks that required
#' fields are present and consistently sized, coerces optional fields, and
#' issues a warning when `lambda` is not nondecreasing.
#'
#' @param x A `confoundsens` object.
#'
#' @return The validated (and possibly coerced) object, invisibly.
#'   Errors if the object is structurally invalid.
#' @keywords internal
validate_confoundsens <- function(x) {
  if (!inherits(x, "confoundsens")) {
    stop("`x` must inherit from class 'confoundsens'.", call. = FALSE)
  }

  # Required fields
  if (is.null(x$lambda) || !is.numeric(x$lambda)) {
    stop("`x$lambda` must be a numeric vector.", call. = FALSE)
  }
  if (is.null(x$theta) || !is.numeric(x$theta)) {
    stop("`x$theta` must be a numeric vector.", call. = FALSE)
  }

  n <- length(x$lambda)
  if (n == 0L) {
    stop("`x$lambda` must not be empty.", call. = FALSE)
  }
  if (length(x$theta) != n) {
    stop("`x$theta` must have the same length as `x$lambda`.", call. = FALSE)
  }
  if (anyNA(x$lambda)) {
    stop("`x$lambda` must not contain NA.", call. = FALSE)
  }

  # Optional: level — coerce and store back
  if (!is.null(x$level)) {
    if (length(x$level) != n) {
      stop("`x$level` must have the same length as `x$lambda`.", call. = FALSE)
    }
    x$level <- as.character(x$level)
  }

  # Optional: se
  if (!is.null(x$se)) {
    if (!is.numeric(x$se) || length(x$se) != n) {
      stop("`x$se` must be a numeric vector with the same length as `x$lambda`.",
           call. = FALSE)
    }
  }

  # Optional: t
  if (!is.null(x$t)) {
    if (!is.numeric(x$t) || length(x$t) != n) {
      stop("`x$t` must be a numeric vector with the same length as `x$lambda`.",
           call. = FALSE)
    }
  }

  # Optional: meta (added in 0.2.0; absent in objects built by 0.1.0)
  if (is.null(x$meta)) x$meta <- list()
  if (!is.list(x$meta)) {
    stop("`x$meta` must be a list.", call. = FALSE)
  }

  # Non-fatal warning for non-monotone lambda
  grp <- if (is.null(x$level)) rep("", n) else x$level
  non_mono <- any(vapply(split(x$lambda, grp),
                         function(l) length(l) >= 2L && any(diff(l) < 0),
                         logical(1)))
  if (non_mono) {
    warning("`x$lambda` is not nondecreasing within level; plots may look odd.",
            call. = FALSE)
  }

  invisible(x)
}
