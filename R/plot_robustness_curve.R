#' Plot a robustness curve
#'
#' Plots the sensitivity path stored in a `confoundsens` object as a function
#' of `lambda`. By default, plots the effect path `theta(lambda)`; optionally
#' plots the test-statistic path `t(lambda)`. Pointwise 95% confidence bands
#' are drawn when standard errors are available.
#'
#' @param x A `confoundsens` object created by [new_confoundsens()] or
#'   [as_confoundsens()].
#' @param what Character; which path to plot: `"theta"` (default) or `"t"`.
#' @param bands Logical; if `TRUE` and `x$se` is present, draw 95% pointwise
#'   confidence bands (applies only when `what = "theta"`).
#' @param points Logical; if `TRUE`, overlay points on the line.
#' @param reference Logical; if `TRUE` (default) and `x$meta` records a null
#'   value (`meta$null`) or a critical value (`meta$threshold`), draw them as
#'   horizontal reference lines (dotted and dashed respectively). Objects
#'   created by [new_confoundsens()] without `meta` are unaffected.
#' @param facet_level Logical; if `TRUE` and `x$level` is present, facet the
#'   plot by level.
#'
#' @return A [ggplot2::ggplot()] object.
#' @export
#'
#' @examples
#' x <- new_confoundsens(
#'   lambda = seq(0, 0.2, length.out = 25),
#'   theta  = 1 - 2 * seq(0, 0.2, length.out = 25),
#'   se     = rep(0.05, 25),
#'   level  = rep(c("within", "between"), length.out = 25)
#' )
#' plot_robustness_curve(x)
#' plot_robustness_curve(x, bands = FALSE, facet_level = FALSE)
plot_robustness_curve <- function(x,
                                  what        = c("theta", "t"),
                                  bands       = TRUE,
                                  points      = TRUE,
                                  facet_level = TRUE,
                                  reference   = TRUE) {
  what <- match.arg(what)

  if (!inherits(x, "confoundsens")) {
    stop("`x` must be an object of class 'confoundsens'.", call. = FALSE)
  }
  x <- validate_confoundsens(x)

  n <- length(x$lambda)

  level <- if (is.null(x$level)) rep("pooled", n) else as.character(x$level)
  se    <- if (is.null(x$se))    rep(NA_real_, n) else as.numeric(x$se)
  tval  <- if (is.null(x$t))     rep(NA_real_, n) else as.numeric(x$t)

  y <- if (what == "theta") as.numeric(x$theta) else tval
  if (what == "t" && all(is.na(y))) {
    stop("`what = 't'` requested but `x$t` is missing or all NA.", call. = FALSE)
  }

  df <- data.frame(
    lambda = as.numeric(x$lambda),
    level  = factor(level, levels = unique(level)),
    y      = y,
    se     = se,
    stringsAsFactors = FALSE
  )

  # Compute bands before plotting so lo/hi exist in the same df
  df$lo <- df$y - 1.96 * df$se
  df$hi <- df$y + 1.96 * df$se

  meta <- x$meta
  ylab <- if (what == "theta") {
    if (!is.null(meta$theta_label)) meta$theta_label else expression(theta(lambda))
  } else expression(t(lambda))
  xlab <- if (!is.null(meta$lambda_label)) meta$lambda_label else expression(lambda)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$lambda, y = .data$y))

  if (isTRUE(reference)) {
    ref <- if (what == "theta") meta$null else if (length(meta)) 0 else NULL
    if (!is.null(ref)) {
      p <- p + ggplot2::geom_hline(yintercept = ref, linetype = "dotted",
                                   colour = "grey40")
    }
    if (what == "theta" && !is.null(meta$threshold)) {
      p <- p + ggplot2::geom_hline(yintercept = meta$threshold,
                                   linetype = "dashed", colour = "grey40")
    }
  }

  p <- p + ggplot2::geom_line() +
    ggplot2::labs(x = xlab, y = ylab)

  if (isTRUE(bands) && what == "theta" && !all(is.na(df$se))) {
    p <- p + ggplot2::geom_ribbon(
      ggplot2::aes(ymin = .data$lo, ymax = .data$hi),
      alpha = 0.2,
      colour = NA
    )
  }

  if (isTRUE(points)) {
    p <- p + ggplot2::geom_point(size = 1.5)
  }

  if (isTRUE(facet_level) && !is.null(x$level)) {
    p <- p + ggplot2::facet_wrap(~ .data$level, scales = "free_y")
  }

  p
}
