#' Three-panel Taylor approximation illustration
#'
#' Produces a three-panel figure (linear / concave-down / convex-up) showing
#' a simulated confounding path together with its first-order (tangent) and
#' second-order (local quadratic) approximations, for three curvature
#' regimes: linear, concave-down, and convex-up.
#'
#' @param delta_max Positive numeric scalar. Upper bound of the delta grid.
#' @param step Positive numeric scalar. Grid step size.
#' @param theta0 Numeric scalar. Baseline effect at delta = 0.
#' @param slope Numeric scalar. First-order slope at delta = 0.
#' @param kappa Non-negative numeric scalar. Curvature magnitude.
#' @param local_max_delta Positive numeric scalar <= `delta_max`. Width of
#'   the local window used for the quadratic approximation.
#'
#' @return A named list with ggplot objects `A` (linear), `B` (concave), and
#'   `C` (convex), returned invisibly. The three panels are also printed to
#'   the active graphics device via [gridExtra::grid.arrange()] if
#'   \pkg{gridExtra} is installed, or via [graphics::layout()] otherwise.
#' @export
#'
#' @details This is a didactic illustration built from simulated paths
#'   ([simulate_taylor_demo()]); it does not analyse data. It was called
#'   `plot_figure2_taylor_panels()` in version 0.1.0; that name still works
#'   but is deprecated.
#'
#' @examples
#' plots <- plot_taylor_panels()
#' # Access individual panels
#' plots$A
plot_taylor_panels <- function(delta_max       = 1.5,
                                       step            = 0.02,
                                       theta0          = 0.40,
                                       slope           = -0.70,
                                       kappa           = 0.40,
                                       local_max_delta = 0.20) {

  is_scalar_pos <- function(z) is.numeric(z) && length(z) == 1L && !is.na(z) && z > 0
  is_scalar_num <- function(z) is.numeric(z) && length(z) == 1L && !is.na(z)

  if (!is_scalar_pos(delta_max))       stop("`delta_max` must be > 0.",          call. = FALSE)
  if (!is_scalar_pos(step))            stop("`step` must be > 0.",               call. = FALSE)
  if (!is_scalar_num(theta0))          stop("`theta0` must be numeric.",          call. = FALSE)
  if (!is_scalar_num(slope))           stop("`slope` must be numeric.",           call. = FALSE)
  if (!is.numeric(kappa) || length(kappa) != 1L || is.na(kappa) || kappa < 0)
    stop("`kappa` must be >= 0.", call. = FALSE)
  if (!is_scalar_pos(local_max_delta)) stop("`local_max_delta` must be > 0.",    call. = FALSE)
  if (local_max_delta > delta_max)     stop("`local_max_delta` must be <= `delta_max`.", call. = FALSE)

  sims <- simulate_taylor_demo(
    delta_max = delta_max, step = step,
    theta0 = theta0, slope = slope, kappa = kappa
  )

  # Build the long-form df that plot_local_taylor() expects
  build_long <- function(cs_obj, kind = c("linear", "concave", "convex")) {
    kind  <- match.arg(kind)
    delta <- cs_obj$lambda
    theta <- cs_obj$theta

    K <- switch(kind,
                linear  = 0,
                concave = -abs(kappa),
                convex  =  abs(kappa))

    keep  <- delta <= local_max_delta
    dloc  <- delta[keep]
    n_loc <- length(dloc)

    data.frame(
      delta  = rep(dloc, 3L),
      series = rep(c("Observed path", "Tangent (1st order)", "Quadratic (2nd order)"),
                   each = n_loc),
      value  = c(
        theta[keep],
        theta0 + slope * dloc,
        theta0 + slope * dloc + 0.5 * K * dloc^2
      ),
      stringsAsFactors = FALSE
    )
  }

  make_panel <- function(cs_obj, title, kind) {
    df_long <- build_long(cs_obj, kind)
    plot_local_taylor(df_long, facet = FALSE) +
      ggplot2::ggtitle(title) +
      ggplot2::theme(legend.position = "bottom")
  }

  pA <- make_panel(sims$linear,  "Panel A: Linear (K = 0)",       "linear")
  pB <- make_panel(sims$concave, "Panel B: Concave-down (K < 0)", "concave")
  pC <- make_panel(sims$convex,  "Panel C: Convex-up (K > 0)",    "convex")

  # Arrange panels — use gridExtra if available, base graphics otherwise.
  # gridExtra is in Suggests; check before calling to satisfy R CMD CHECK.
  if (requireNamespace("gridExtra", quietly = TRUE)) {
    gridExtra::grid.arrange(pA, pB, pC, nrow = 1L)
  } else {
    message(
      "Install the 'gridExtra' package for side-by-side panel layout. ",
      "Falling back to sequential base-graphics output."
    )
    op <- graphics::par(mfrow = c(1L, 3L))
    on.exit(graphics::par(op), add = TRUE)
    print(pA)
    print(pB)
    print(pC)
  }

  invisible(list(A = pA, B = pB, C = pC))
}


#' @rdname plot_taylor_panels
#' @usage NULL
#' @export
plot_figure2_taylor_panels <- function(...) {
  .Deprecated("plot_taylor_panels", package = "confoundvis")
  plot_taylor_panels(...)
}
