#' Summary method for confoundsens objects
#'
#' Produces a concise numerical summary of the sensitivity path, optionally
#' broken down by level.
#'
#' @param object A `confoundsens` object.
#' @param ... Unused; included for S3 method consistency.
#'
#' @return An object of class `summary.confoundsens` containing a `table`
#'   data.frame with per-level summary statistics.
#' @export
#'
#' @examples
#' x <- new_confoundsens(
#'   lambda = seq(0, 0.2, length.out = 10),
#'   theta  = seq(1, 0.6, length.out = 10),
#'   se     = rep(0.08, 10),
#'   level  = rep(c("female", "male"), length.out = 10)
#' )
#' summary(x)
summary.confoundsens <- function(object, ...) {
  x <- validate_confoundsens(object)
  n <- length(x$lambda)

  by_level <- !is.null(x$level)
  levs     <- if (by_level) unique(x$level) else "pooled"

  make_one <- function(idx, level_name) {
    lam <- x$lambda[idx]
    th  <- x$theta[idx]
    se  <- if (is.null(x$se)) rep(NA_real_, length(idx)) else x$se[idx]
    tv  <- if (is.null(x$t))  rep(NA_real_, length(idx)) else x$t[idx]

    # Use [1L] to guard against ties in lambda
    i0 <- which.min(lam)[1L]

    data.frame(
      level      = level_name,
      n          = length(idx),
      lambda_min = min(lam, na.rm = TRUE),
      lambda_max = max(lam, na.rm = TRUE),
      theta_min  = min(th,  na.rm = TRUE),
      theta_max  = max(th,  na.rm = TRUE),
      theta_0    = th[i0],
      se_0       = se[i0],
      t_0        = tv[i0],
      stringsAsFactors = FALSE
    )
  }

  out <- do.call(rbind, lapply(levs, function(L) {
    idx <- if (by_level) which(x$level == L) else seq_len(n)
    make_one(idx, L)
  }))

  res <- list(table = out)
  class(res) <- "summary.confoundsens"
  res
}

#' @export
print.summary.confoundsens <- function(x, ...) {
  cat("<summary.confoundsens>\n\n")
  print(x$table, row.names = FALSE)
  invisible(x)
}
