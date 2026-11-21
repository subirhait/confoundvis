#' Run every applicable sensitivity analysis for one estimate
#'
#' One call that computes the partial-\eqn{R^2}, impact-threshold (ITCV), and
#' replacement (RIR) analyses for the coefficient of `treatment` in a fitted
#' linear model, adds any further paths the user supplies (for example an
#' E-value path from [from_evalue()] or a Rosenbaum path from
#' [sens_path_pairs()]), and collects them in one object. The object prints
#' as a side-by-side table ([sens_table()]) and plots as one figure with one
#' panel per analysis. Each analysis keeps its own sensitivity parameter, its
#' own scale, and its own decision rule; the panels do not share an axis, and
#' the values in different rows are not to be compared or ranked.
#'
#' The default grids are chosen from the closed-form values (robustness
#' value, ITCV, RIR proportion), so that each crossing falls well inside the
#' grid. They use 2,001 points, which places the interpolated crossings
#' within about \eqn{10^{-6}} of the closed-form values in typical cases.
#'
#' @param model A fitted `lm` object.
#' @param treatment Character string naming the treatment coefficient.
#' @param alpha Two-sided significance level for all three regression
#'   analyses.
#' @param extra Optional named list of further `confoundsens` paths to add to
#'   the profile.
#' @param n_grid Number of grid points for each default path.
#'
#' @return An object of class `confoundvis_profile`: a list with `paths` (a
#'   named list of `confoundsens` objects), `table` (the [sens_table()]
#'   output), `itcv` (the full [itcv_lm()] result), and `rir` (the full
#'   [rir_lm()] result).
#' @seealso [sens_table()], [sens_report()], [plot_robustness_curve()]
#' @export
#'
#' @examples
#' fit <- lm(mpg ~ am + wt + hp, data = mtcars)
#' pr <- sens_profile(fit, "wt")
#' pr
#' plot(pr)
sens_profile <- function(model, treatment, alpha = 0.05, extra = NULL,
                         n_grid = 2001) {
  pr <- .lm_parts(model, treatment)
  if (!is.numeric(alpha) || length(alpha) != 1L || alpha <= 0 || alpha >= 1) {
    stop("`alpha` must be a single number in (0, 1).", call. = FALSE)
  }
  if (!is.null(extra)) {
    if (!is.list(extra) || inherits(extra, "confoundsens") ||
        !all(vapply(extra, inherits, logical(1), what = "confoundsens"))) {
      stop("`extra` must be a named list of `confoundsens` objects.", call. = FALSE)
    }
    if (is.null(names(extra)) || any(names(extra) == "")) {
      stop("Every element of `extra` must be named.", call. = FALSE)
    }
  }
  rv <- .rv(pr$t, pr$dof)
  r  <- pr$t / sqrt(pr$t^2 + pr$dof)
  up_r2   <- min(0.999, max(0.01, 2 * rv))
  up_itcv <- min(0.999, max(0.01, 1.5 * abs(r)))
  p_r2 <- sens_path_lm(model, treatment, lambda = seq(0, up_r2, length.out = n_grid))
  it   <- itcv_lm(model, treatment, alpha = alpha,
                  lambda = seq(0, up_itcv, length.out = n_grid))
  rr   <- rir_lm(model, treatment, alpha = alpha,
                 lambda = seq(0, 1, length.out = n_grid))
  paths <- c(list(partial_r2 = p_r2, itcv = it$path, rir = rr$path), extra)
  tb <- sens_table(paths, alpha = alpha)
  structure(list(paths = paths, table = tb, itcv = it, rir = rr,
                 treatment = treatment, alpha = alpha),
            class = "confoundvis_profile")
}

#' @export
print.confoundvis_profile <- function(x, ...) {
  cat(sprintf("<confoundvis_profile> treatment: %s; %d analyses\n\n",
              x$treatment, length(x$paths)))
  tb <- x$table
  keep <- c("path", "summary_name", "summary_value", "lambda_null",
            "lambda_sig", "alpha")
  out <- tb[, keep]
  class(out) <- "data.frame"
  out$summary_value <- signif(out$summary_value, 4)
  out$lambda_null <- signif(out$lambda_null, 4)
  out$lambda_sig <- signif(out$lambda_sig, 4)
  print(out, row.names = FALSE, ...)
  cat("\nlambda (the sensitivity parameter) in each row:\n")
  for (i in seq_len(nrow(tb))) {
    cat(sprintf("  %-11s %s\n", tb$path[i],
                if (is.na(tb$parameter[i])) "-" else tb$parameter[i]))
  }
  cat("Values are comparable only within a row; do not rank the rows.\n")
  invisible(x)
}

#' Plot a sensitivity profile
#'
#' One panel per analysis, each with its own axes. The panel title names the
#' framework and states what the horizontal axis (x, the sensitivity
#' parameter) and the vertical axis (y, the path value) measure in that
#' panel. The dotted line is the null value, the dashed line is the stored
#' threshold (when the analysis has one), and the vertical long-dashed line
#' marks the confounder strength at which significance is lost.
#'
#' @param x A `confoundvis_profile` object from [sens_profile()].
#' @param ... Unused.
#' @return A [ggplot2::ggplot()] object.
#' @export
plot.confoundvis_profile <- function(x, ...) {
  nm <- names(x$paths)
  pts <- x$table
  fw_name <- c(partial_r2 = "Partial R2", itcv = "ITCV", rir = "RIR",
               evalue = "E-value", rosenbaum = "Rosenbaum")
  lab <- function(i) {
    m <- x$paths[[i]]$meta
    f <- if (!is.null(m$framework) && m$framework %in% names(fw_name))
      fw_name[[m$framework]] else nm[i]
    wrap <- function(z) paste(strwrap(z, width = 48, exdent = 3), collapse = "\n")
    paste(sprintf("%s. %s", LETTERS[i], f),
          wrap(paste("x:", if (is.null(m$lambda_label)) "lambda" else m$lambda_label)),
          wrap(paste("y:", if (is.null(m$theta_label)) "theta" else m$theta_label)),
          sep = "\n")
  }
  panels <- vapply(seq_along(nm), lab, character(1))
  dat <- do.call(rbind, lapply(seq_along(nm), function(i) {
    p <- x$paths[[i]]
    data.frame(panel = panels[i], lambda = p$lambda, theta = p$theta,
               stringsAsFactors = FALSE)
  }))
  ref <- do.call(rbind, lapply(seq_along(nm), function(i) {
    m <- x$paths[[i]]$meta
    data.frame(panel = panels[i],
               null = if (is.null(m$null)) NA_real_ else m$null,
               threshold = if (is.null(m$threshold)) NA_real_ else m$threshold,
               lambda_sig = pts$lambda_sig[pts$path == nm[i]][1L],
               stringsAsFactors = FALSE)
  }))
  dat$panel <- factor(dat$panel, levels = panels)
  ref$panel <- factor(ref$panel, levels = panels)
  ggplot2::ggplot(dat, ggplot2::aes(x = .data$lambda, y = .data$theta)) +
    ggplot2::geom_hline(data = ref[!is.na(ref$null), ],
                        ggplot2::aes(yintercept = .data$null),
                        linetype = "dotted", colour = "grey40") +
    ggplot2::geom_hline(data = ref[!is.na(ref$threshold), ],
                        ggplot2::aes(yintercept = .data$threshold),
                        linetype = "dashed", colour = "grey40") +
    ggplot2::geom_vline(data = ref[!is.na(ref$lambda_sig), ],
                        ggplot2::aes(xintercept = .data$lambda_sig),
                        linetype = "longdash", colour = "grey60") +
    ggplot2::geom_line() +
    ggplot2::facet_wrap(~ .data$panel, scales = "free",
                        ncol = min(2L, length(panels))) +
    ggplot2::labs(x = "Sensitivity parameter (own units in each panel)",
                  y = "Path value (own scale in each panel)")
}
