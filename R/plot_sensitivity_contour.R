#' Sensitivity contour plot
#'
#' Draws sensitivity boundaries in the plane of the two confounder
#' correlations \eqn{(r_{YU}, r_{DU})} and optionally overlays observed
#' covariate benchmarks as labelled points.
#'
#' Two kinds of boundary are available.
#'
#' * **Product reference** (always drawn, dashed when an exact boundary is
#'   also drawn): the four branches of \eqn{|r_{YU} r_{DU}| = }`threshold`.
#'   With `threshold` equal to the ITCV, this is the impact threshold of
#'   Frank (2000). It is exact for confounders whose two correlations have
#'   equal magnitude, but for unequal correlations the product alone does not
#'   determine whether significance is lost.
#' * **Exact boundary** (when `itcv` is supplied): the set of pairs at which
#'   the confounder-adjusted partial correlation
#'   \deqn{r_{\mathrm{adj}} = \frac{r - r_{DU} r_{YU}}
#'         {\sqrt{(1 - r_{DU}^2)(1 - r_{YU}^2)}}}{r_adj = (r - r_DU r_YU) / sqrt((1 - r_DU^2)(1 - r_YU^2))}
#'   has magnitude equal to the critical correlation. The solid curve is the
#'   boundary in the direction of the original estimate; beyond it the
#'   estimate is no longer significant ("Not significant" shading). The
#'   dotted curve is the boundary in the reversed direction; beyond it the
#'   adjusted estimate is significant with the opposite sign. Both curves are
#'   computed in closed form. Pairs for which
#'   the correlation matrix of \eqn{(D, Y, U)} given the covariates would not
#'   be positive definite, \eqn{1 - r^2 - r_{DU}^2 - r_{YU}^2 + 2 r r_{DU} r_{YU} \le 0}{1 - r^2 - r_DU^2 - r_YU^2 + 2 r r_DU r_YU <= 0},
#'   cannot occur and are shaded as infeasible.
#'
#' @param threshold Positive numeric scalar; the product threshold. The
#'   product reference satisfies \eqn{|r_{YU} \cdot r_{DU}| = }`threshold`.
#'   May be omitted when `itcv` is supplied, in which case the conditional
#'   ITCV is used.
#' @param grid_n Integer >= 50; resolution of the curves (points per branch
#'   for the product reference; grid points per axis for the exact boundary
#'   and shading).
#' @param benchmarks Optional data.frame with columns `r_yu` and `r_du`
#'   (numeric) and an optional `label` column (character). Each row is plotted
#'   as a labelled benchmark point; labels repel one another when the
#'   \pkg{ggrepel} package is installed. The output of [covariate_impacts()] with
#'   `metric = "itcv"` can be passed directly.
#' @param itcv Optional output of [itcv_lm()] (a list with `r` and `r_crit`),
#'   or a list or named numeric vector with elements `r` and `r_crit`. When
#'   supplied, the exact boundary and the infeasible region are drawn.
#' @param limits Numeric vector of length two; the displayed range of both
#'   axes. Narrow it (for example `c(-0.3, 0.3)`) when benchmarks lie close to
#'   the origin.
#'
#' @return A [ggplot2::ggplot()] object.
#' @references Frank, K. A. (2000). Impact of a confounding variable on a
#'   regression coefficient. *Sociological Methods & Research*, 29(2),
#'   147--194. \doi{10.1177/0049124100029002001}
#' @export
#'
#' @examples
#' b <- data.frame(
#'   r_yu  = c(0.10, 0.15),
#'   r_du  = c(0.20, 0.12),
#'   label = c("SES", "BMI")
#' )
#' plot_sensitivity_contour(threshold = 0.02, benchmarks = b)
#'
#' # Exact boundary and infeasible region from a fitted model
#' fit <- lm(mpg ~ am + wt + hp, data = mtcars)
#' plot_sensitivity_contour(itcv = itcv_lm(fit, "wt"))
plot_sensitivity_contour <- function(threshold  = NULL,
                                     grid_n     = 200,
                                     benchmarks = NULL,
                                     itcv       = NULL,
                                     limits     = c(-1, 1)) {
  exact <- !is.null(itcv)
  if (exact) {
    r      <- as.numeric(itcv[["r"]])
    r_crit <- abs(as.numeric(itcv[["r_crit"]]))
    if (length(r) != 1L || length(r_crit) != 1L || !is.finite(r) ||
        !is.finite(r_crit) || abs(r) >= 1 || r_crit >= 1) {
      stop("`itcv` must provide single values `r` and `r_crit` in (-1, 1).",
           call. = FALSE)
    }
    if (is.null(threshold)) threshold <- (abs(r) - r_crit) / (1 - r_crit)
  }
  if (!is.numeric(threshold) || length(threshold) != 1L ||
      is.na(threshold) || threshold <= 0) {
    stop("`threshold` must be a single positive numeric value ",
         "(or supply `itcv` for an initially significant estimate).", call. = FALSE)
  }
  if (!is.numeric(grid_n) || length(grid_n) != 1L || grid_n < 50) {
    stop("`grid_n` must be a single integer >= 50.", call. = FALSE)
  }
  grid_n <- as.integer(grid_n)
  if (!is.numeric(limits) || length(limits) != 2L || anyNA(limits) ||
      limits[1] >= limits[2] || any(abs(limits) > 1)) {
    stop("`limits` must be two increasing numbers in [-1, 1].", call. = FALSE)
  }

  # Product reference: four branches of r_du = +/- threshold / r_yu.
  eps <- 1e-6
  rr  <- seq(eps, 1, length.out = grid_n)
  make_branch <- function(r_yu_vals, r_du_vals, grp) {
    df <- data.frame(r_yu = r_yu_vals, r_du = r_du_vals, group = grp,
                     stringsAsFactors = FALSE)
    df[abs(df$r_du) <= 1 & abs(df$r_yu) <= 1, , drop = FALSE]
  }
  boundary <- rbind(
    make_branch( rr,  threshold / rr,  "Q1"),
    make_branch(-rr, -threshold / rr,  "Q2"),
    make_branch( rr, -threshold / rr,  "Q3"),
    make_branch(-rr,  threshold / rr,  "Q4")
  )

  p <- ggplot2::ggplot(
    boundary,
    ggplot2::aes(x = .data$r_yu, y = .data$r_du, group = .data$group))

  if (exact) {
    step <- 2 / grid_n
    g  <- (seq_len(grid_n) - 0.5) * step - 1
    gr <- expand.grid(r_yu = g, r_du = g)
    det <- 1 - r^2 - gr$r_yu^2 - gr$r_du^2 + 2 * r * gr$r_yu * gr$r_du
    radj <- (r - gr$r_yu * gr$r_du) /
      sqrt((1 - gr$r_yu^2) * (1 - gr$r_du^2))
    f <- sign(r) * radj
    gr$zone <- ifelse(det <= 0, "Infeasible",
               ifelse(f <= -r_crit, "Significant, reversed direction",
               ifelse(f <  r_crit,  "Not significant", NA_character_)))
    p <- p + ggplot2::geom_raster(
      data = gr,
      ggplot2::aes(x = .data$r_yu, y = .data$r_du, fill = .data$zone),
      inherit.aes = FALSE, na.rm = TRUE) +
      ggplot2::scale_fill_manual(
        values = c("Not significant" = "#bdd7e7",
                   "Significant, reversed direction" = "#fcbba1",
                   "Infeasible" = "grey60"),
        breaks = c("Not significant", "Significant, reversed direction",
                   "Infeasible"),
        na.value = NA, name = NULL)

    ex <- .exact_boundary(r, r_crit, n = 8L * grid_n)
    if (nrow(ex)) {
      p <- p + ggplot2::geom_path(
        data = ex,
        ggplot2::aes(x = .data$r_yu, y = .data$r_du, group = .data$group,
                     linetype = .data$kind),
        inherit.aes = FALSE, linewidth = 0.7) +
        ggplot2::scale_linetype_manual(
          values = c("Original direction" = "solid",
                     "Reversed direction" = "dotted"),
          name = "Boundary |r_adj| = critical r")
    }
  }

  p <- p +
    ggplot2::geom_path(linetype = if (exact) "dashed" else "solid") +
    ggplot2::coord_equal(xlim = limits, ylim = limits) +
    ggplot2::labs(
      x = expression(r[YU]),
      y = expression(r[DU]),
      caption = if (exact)
        sprintf("Exact boundaries for r = %s, critical r = %s\nDashed: product reference |r_YU r_DU| = %s",
                trimws(formatC(r, digits = 3, format = "g")),
                trimws(formatC(r_crit, digits = 3, format = "g")),
                trimws(formatC(threshold, digits = 3, format = "g")))
      else NULL
    )

  if (!is.null(benchmarks)) {
    if (!is.data.frame(benchmarks) ||
        !all(c("r_yu", "r_du") %in% names(benchmarks))) {
      stop("`benchmarks` must be a data.frame with columns `r_yu` and `r_du`.",
           call. = FALSE)
    }
    if (!is.numeric(benchmarks$r_yu) || !is.numeric(benchmarks$r_du)) {
      stop("`benchmarks$r_yu` and `benchmarks$r_du` must be numeric.",
           call. = FALSE)
    }
    benchmarks <- benchmarks[!is.na(benchmarks$r_yu) & !is.na(benchmarks$r_du), ,
                             drop = FALSE]

    p <- p + ggplot2::geom_point(
      data = benchmarks,
      ggplot2::aes(x = .data$r_yu, y = .data$r_du),
      inherit.aes = FALSE,
      size = 2.5
    )

    if ("label" %in% names(benchmarks)) {
      if (requireNamespace("ggrepel", quietly = TRUE)) {
        p <- p + ggrepel::geom_text_repel(
          data = benchmarks,
          ggplot2::aes(x = .data$r_yu, y = .data$r_du, label = .data$label),
          inherit.aes = FALSE,
          size = 3.5,
          min.segment.length = 0,
          seed = 1
        )
      } else {
        p <- p + ggplot2::geom_text(
          data = benchmarks,
          ggplot2::aes(x = .data$r_yu, y = .data$r_du, label = .data$label),
          inherit.aes = FALSE,
          nudge_y = 0.04,
          size = 3.5
        )
      }
    }
  }

  p
}


# Internal: exact boundaries |r_adj| = r_crit in closed form.
# For r > 0, write a = r_DU, b = r_YU, R = |r|, c = r_crit. Squaring
# (R - a b) = +/- c sqrt((1 - a^2)(1 - b^2)) gives, for each b, the quadratic
#   (b^2 + c^2 (1 - b^2)) a^2 - 2 R b a + (R^2 - c^2 (1 - b^2)) = 0.
# Roots with R - a b > 0 lie on the boundary in the original direction;
# roots with R - a b < 0 lie on the reversed-direction boundary. For r < 0
# the figure is reflected (a -> -a).
.exact_boundary <- function(r, r_crit, n = 1600L) {
  R  <- abs(r); cc <- r_crit
  b  <- seq(-1, 1, length.out = n + 2L)[-c(1L, n + 2L)]
  q  <- cc^2 * (1 - b^2)
  A  <- b^2 + q
  disc <- (R * b)^2 - A * (R^2 - q)
  out <- list()
  for (sgn in c(-1, 1)) {
    a  <- (R * b + sgn * sqrt(pmax(disc, 0))) / A
    ok <- disc >= 0 & abs(a) < 1
    kind <- ifelse(R - a * b > 0, "Original direction", "Reversed direction")
    key  <- ifelse(ok, kind, NA_character_)
    k2   <- ifelse(is.na(key), "", key)
    run  <- cumsum(c(TRUE, k2[-1L] != k2[-length(k2)]))
    run[is.na(key)] <- NA
    keep <- !is.na(run)
    if (any(keep)) {
      out[[length(out) + 1L]] <- data.frame(
        r_yu = b[keep], r_du = sign(r) * a[keep],
        kind = kind[keep], group = paste0(sgn, "_", run[keep]),
        stringsAsFactors = FALSE)
    }
  }
  if (!length(out)) {
    return(data.frame(r_yu = numeric(0), r_du = numeric(0),
                      kind = character(0), group = character(0)))
  }
  do.call(rbind, out)
}
