#' Reversal cone geometry plot
#'
#' Visualizes the two-dimensional cross-section of the reversal cone
#' \eqn{C_l} at a fixed confounding effect magnitude \eqn{|\delta|}. The cone
#' partitions the \eqn{(q, p)} confounding parameter space into an
#' **attenuation zone** (where the effect shrinks but does not reverse sign)
#' and a **reversal zone** (where the sign changes). The boundary between
#' zones is the reversal curve.
#'
#' @details The display uses the stylized bilinear bias model
#' \eqn{\theta_{adj} = \theta_0 - |\delta|\, q\, p}, so the reversal
#' boundary is the hyperbola \eqn{p = \theta_0 / (|\delta| q)}. This is a
#' conceptual map of the parameter space, not an estimate: the annotated
#' "reversal share" is the fraction of the *displayed* grid that falls in the
#' reversal zone and therefore depends on `q_range` and `p_range`. It should
#' not be reported as a probability or as a property of the data.
#'
#' @param theta0 Numeric; baseline estimated effect at zero confounding (must
#'   be nonzero).
#' @param delta Numeric; the fixed confounding effect magnitude \eqn{|\delta|}
#'   at which the cross-section is evaluated (> 0).
#' @param q_range Numeric vector of length 2; range of the confounding
#'   prevalence parameter \eqn{q \in [0, 1]} (default `c(0, 1)`).
#' @param p_range Numeric vector of length 2; range of the confounding
#'   impact parameter \eqn{p} (default `c(-1, 1)`).
#' @param grid_n Integer; grid resolution for each axis (default 200).
#' @param show_boundary Logical; draw the reversal boundary curve.
#' @param show_volume Logical; annotate with the share of the displayed grid
#'   lying in the reversal zone.
#'
#' @return A ggplot object.
#' @export
#'
#' @examples
#' plot_reversal_cone(theta0 = 0.40, delta = 0.20)
#' plot_reversal_cone(theta0 = 0.40, delta = 0.50)
plot_reversal_cone <- function(theta0        = 0.40,
                               delta         = 0.20,
                               q_range       = c(0, 1),
                               p_range       = c(-1, 1),
                               grid_n        = 200L,
                               show_boundary = TRUE,
                               show_volume   = TRUE) {

  if (!is.numeric(theta0) || length(theta0) != 1L || is.na(theta0) || theta0 == 0)
    stop("`theta0` must be a single nonzero numeric value.", call. = FALSE)
  if (!is.numeric(delta) || length(delta) != 1L || is.na(delta) || delta <= 0)
    stop("`delta` must be a single positive numeric value.", call. = FALSE)
  if (!is.numeric(q_range) || length(q_range) != 2L || any(is.na(q_range)) ||
      q_range[1] >= q_range[2])
    stop("`q_range` must be a numeric vector of length 2 with q_range[1] < q_range[2].",
         call. = FALSE)
  if (!is.numeric(p_range) || length(p_range) != 2L || any(is.na(p_range)) ||
      p_range[1] >= p_range[2])
    stop("`p_range` must be a numeric vector of length 2 with p_range[1] < p_range[2].",
         call. = FALSE)
  grid_n <- as.integer(grid_n)
  if (is.na(grid_n) || grid_n < 20L)
    stop("`grid_n` must be an integer >= 20.", call. = FALSE)

  q_seq <- seq(q_range[1], q_range[2], length.out = grid_n)
  p_seq <- seq(p_range[1], p_range[2], length.out = grid_n)

  grid <- expand.grid(q = q_seq, p = p_seq)

  # The adjusted estimate under confounding:
  # theta_adj = theta0 - delta * q * p
  # Reversal occurs when sign(theta_adj) != sign(theta0)
  # Attenuation: |theta_adj| < |theta0| but same sign
  # Reversal:    theta_adj * theta0 < 0
  grid$theta_adj <- theta0 - delta * grid$q * grid$p
  grid$zone <- ifelse(
    grid$theta_adj * theta0 < 0,
    "Reversal",
    ifelse(abs(grid$theta_adj) < abs(theta0), "Attenuation", "Amplification")
  )
  grid$zone <- factor(grid$zone, levels = c("Amplification", "Attenuation", "Reversal"))

  zone_colors <- c(
    "Amplification" = "#AEC6CF",
    "Attenuation"   = "#FFF3B0",
    "Reversal"      = "#FF6B6B"
  )

  # Proportion of grid in reversal zone (proxy for cone volume)
  vol_prop <- mean(grid$zone == "Reversal")

  p <- ggplot2::ggplot(grid, ggplot2::aes(x = .data$q, y = .data$p,
                                           fill = .data$zone)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_manual(values = zone_colors, name = "Zone",
                               drop = FALSE) +
    ggplot2::labs(
      x       = expression(q ~ "(confounding prevalence)"),
      y       = expression(p ~ "(confounding impact)"),
      title   = "Reversal Cone Cross-Section",
      subtitle = sprintf(
        "|delta| = %.2f  |  theta[0] = %.2f", delta, theta0
      )
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::coord_cartesian(expand = FALSE)

  if (isTRUE(show_boundary)) {
    # Reversal boundary: theta0 = delta * q * p  =>  p = theta0 / (delta * q)
    q_bnd <- q_seq[q_seq > 0]
    bnd <- data.frame(
      q = q_bnd,
      p = theta0 / (delta * q_bnd)
    )
    bnd <- bnd[bnd$p >= p_range[1] & bnd$p <= p_range[2], ]

    if (nrow(bnd) > 1L) {
      p <- p + ggplot2::geom_line(
        data = bnd,
        ggplot2::aes(x = .data$q, y = .data$p),
        inherit.aes = FALSE,
        color = "black", linewidth = 1, linetype = "dashed"
      )
    }
  }

  if (isTRUE(show_volume)) {
    p <- p + ggplot2::annotate(
      "text",
      x = q_range[1] + 0.05 * diff(q_range),
      y = p_range[2] - 0.05 * diff(p_range),
      label = sprintf("Reversal share\nof displayed grid\n= %.1f%%", 100 * vol_prop),
      hjust = 0, vjust = 1, size = 3.2, color = "black"
    )
  }

  p
}


#' Compare reversal cones across multilevel components
#'
#' Produces a two-panel plot comparing the reversal cone cross-sections for
#' the within-cluster and between-cluster confounding components in multilevel
#' settings. Each panel calls [plot_reversal_cone()] and they are displayed
#' side by side using faceting.
#'
#' @param theta0_within Numeric; within-cluster baseline effect.
#' @param theta0_between Numeric; between-cluster baseline effect.
#' @param delta_within Numeric; within-cluster confounding magnitude (> 0).
#' @param delta_between Numeric; between-cluster confounding magnitude (> 0).
#' @param ... Additional arguments passed to [plot_reversal_cone()].
#'
#' @return Invisibly, a list with ggplot objects `within` and `between`.
#'   The plots are drawn side-by-side to the current device.
#' @export
#'
#' @examples
#' plot_cone_comparison(
#'   theta0_within  = 0.35, theta0_between = 0.50,
#'   delta_within   = 0.20, delta_between  = 0.30
#' )
plot_cone_comparison <- function(theta0_within  = 0.35,
                                 theta0_between = 0.50,
                                 delta_within   = 0.20,
                                 delta_between  = 0.30,
                                 ...) {

  pw <- plot_reversal_cone(theta0 = theta0_within,  delta = delta_within,  ...) +
          ggplot2::ggtitle("Within-cluster reversal cone")
  pb <- plot_reversal_cone(theta0 = theta0_between, delta = delta_between, ...) +
          ggplot2::ggtitle("Between-cluster reversal cone")

  grid::grid.newpage()
  lay <- grid::grid.layout(nrow = 1L, ncol = 2L)
  grid::pushViewport(grid::viewport(layout = lay))
  print(pw, vp = grid::viewport(layout.pos.row = 1L, layout.pos.col = 1L))
  print(pb, vp = grid::viewport(layout.pos.row = 1L, layout.pos.col = 2L))

  invisible(list(within = pw, between = pb))
}
