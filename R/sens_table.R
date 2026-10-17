# ---------------------------------------------------------------------------
# One table for several sensitivity paths
# ---------------------------------------------------------------------------

# Internal: the framework's own summary statistic, with its name.
.framework_summary <- function(meta) {
  fw <- if (is.null(meta$framework)) "" else meta$framework
  switch(fw,
    partial_r2 = list(name = "RV (point null)",          value = meta$rv),
    itcv       = list(name = "ITCV (conditional)",       value = meta$itcv),
    evalue     = list(name = "E-value (point)",          value = meta$evalue),
    rosenbaum  = list(name = "Sensitivity value (Gamma)", value = meta$sens_value),
    list(name = NA_character_, value = NA_real_))
}

#' Summarize several sensitivity paths in one table
#'
#' Collects the decision information of several `confoundsens` objects into
#' one data frame: for each path (and each level within a path), the
#' framework, the meaning of its sensitivity parameter, the framework's own
#' summary statistic, and the confounder strengths at which the null and the
#' significance boundary are reached, together with the rule used to locate
#' the latter. The table is intended for reporting several analyses of the
#' same study side by side. It does not convert the sensitivity parameters to
#' a common scale: values in the `lambda_null` and `lambda_sig` columns are
#' comparable only within a row.
#'
#' @param ... `confoundsens` objects, or a single list of them. Names (from
#'   `name = object` arguments or list names) are used as row labels.
#' @param alpha Two-sided significance level passed to [robustness_points()]
#'   for paths without a stored threshold. Paths with a stored threshold
#'   (ITCV, Rosenbaum) use the level they were built with.
#'
#' @return A data frame of class `confoundvis_table` with columns `path`,
#'   `level`, `framework`, `source`, `parameter` (meaning of lambda),
#'   `summary_name`, `summary_value`, `theta_start`, `lambda_null`,
#'   `lambda_sig`, and `sig_rule`. It prints with a reminder that rows use
#'   different scales.
#' @seealso [robustness_points()], [sens_report()]
#' @export
#'
#' @examples
#' fit <- lm(mpg ~ am + wt + hp, data = mtcars)
#' g <- seq(0, 0.9, by = 0.001)
#' sens_table(partial_r2 = sens_path_lm(fit, "wt", lambda = g),
#'            itcv       = itcv_lm(fit, "wt", lambda = g)$path)
sens_table <- function(..., alpha = 0.05) {
  paths <- list(...)
  if (length(paths) == 1L && is.list(paths[[1L]]) &&
      !inherits(paths[[1L]], "confoundsens")) {
    paths <- paths[[1L]]
  }
  if (length(paths) == 0L) {
    stop("Supply at least one `confoundsens` object.", call. = FALSE)
  }
  ok <- vapply(paths, inherits, logical(1), what = "confoundsens")
  if (!all(ok)) {
    stop("All inputs must be `confoundsens` objects.", call. = FALSE)
  }
  nm <- names(paths)
  if (is.null(nm)) nm <- rep("", length(paths))
  nm[nm == ""] <- paste0("path", seq_along(paths))[nm == ""]

  rows <- lapply(seq_along(paths), function(i) {
    x  <- validate_confoundsens(paths[[i]])
    m  <- x$meta
    rp <- robustness_points(x, alpha = alpha)
    fs <- .framework_summary(m)
    data.frame(
      path          = nm[i],
      level         = rp$level,
      framework     = if (is.null(m$framework)) NA_character_ else m$framework,
      source        = if (is.null(m$source)) NA_character_ else m$source,
      parameter     = if (is.null(m$lambda_label)) NA_character_ else m$lambda_label,
      summary_name  = fs$name,
      summary_value = if (is.null(fs$value)) NA_real_ else as.numeric(fs$value),
      theta_start   = rp$theta_start,
      lambda_null   = rp$lambda_null,
      lambda_sig    = rp$lambda_sig,
      sig_rule      = rp$sig_rule,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  class(out) <- c("confoundvis_table", "data.frame")
  out
}

#' @export
print.confoundvis_table <- function(x, digits = 4, ...) {
  d <- as.data.frame(unclass(x), stringsAsFactors = FALSE)
  if (all(d$level == "pooled")) d$level <- NULL
  num <- vapply(d, is.numeric, logical(1))
  d[num] <- lapply(d[num], function(v) signif(v, digits))
  cat("<confoundvis_table>\n")
  print(d, row.names = FALSE, right = FALSE)
  cat("\nEach row keeps its own sensitivity parameter (column `parameter`);",
      "lambda values are comparable only within a row.\n")
  invisible(x)
}
