# ---------------------------------------------------------------------------
# Two-parameter sensitivity surfaces
#
# Most sensitivity frameworks describe an omitted confounder with two
# parameters, one for its association with treatment and one for its
# association with the outcome. A confoundsens path is one slice through
# that plane (for example, the equal-strength line). sens_surface() returns
# the whole plane for a path, with the boundaries at which the conclusion
# changes computed in closed form or by one-dimensional root finding from
# the published formula of each framework:
#   partial R2  Cinelli & Hazlett (2020), bias and adjusted SE
#   ITCV        Frank (2000), adjusted partial correlation
#   E-value     Ding & VanderWeele (2016), joint bounding factor
#   Rosenbaum   Rosenbaum & Silber (2009), amplification of Gamma
# ---------------------------------------------------------------------------

#' Two-parameter sensitivity surface
#'
#' Returns the two-parameter sensitivity surface of which a `confoundsens`
#' path is one slice. The horizontal axis is always the confounder's
#' association with treatment and the vertical axis its association with the
#' outcome, each in the units of the framework:
#'
#' | Framework | `a` (treatment side) | `b` (outcome side) | Value |
#' |---|---|---|---|
#' | `partial_r2` | \eqn{R^2_{D \sim U \mid X}} | \eqn{R^2_{Y \sim U \mid D, X}} | adjusted estimate |
#' | `itcv` | \eqn{r_{DU \mid X}} | \eqn{r_{YU \mid X}} | adjusted partial correlation |
#' | `evalue` | \eqn{RR_{UD}} | \eqn{RR_{UY}} | bound on the risk ratio |
#' | `rosenbaum` | \eqn{\Lambda} | \eqn{\Delta} | implied \eqn{\Gamma} |
#'
#' Each surface carries boundaries computed from the published formula of the
#' framework, in closed form where one exists and otherwise by root finding
#' (the partial-\eqn{R^2} significance boundaries, and the Rosenbaum
#' \eqn{\Gamma} at which the upper-bound P-value is 0.5): the
#' **null boundary**, where
#' the point estimate (or, for Rosenbaum, the upper-bound P-value of 0.5) is
#' reached, and the **significance boundary**, where the conclusion at level
#' `alpha` is lost. For the partial-\eqn{R^2} framework, the bias and
#' adjusted standard error of Cinelli and Hazlett (2020) are used, with
#' \eqn{d - 1} degrees of freedom for the significance boundary. For the ITCV,
#' the boundaries are those of [plot_sensitivity_contour()]. For E-values,
#' the joint bounding factor \eqn{B = RR_{UD} RR_{UY} / (RR_{UD} + RR_{UY} - 1)}
#' of Ding and VanderWeele (2016) is used, so the null boundary passes through
#' the point (E-value, E-value). For Rosenbaum's model, \eqn{\Gamma =
#' (\Lambda\Delta + 1) / (\Lambda + \Delta)} (Rosenbaum and Silber, 2009), so
#' the significance boundary is the set of \eqn{(\Lambda, \Delta)} pairs that
#' amplify the sensitivity value.
#'
#' The surface only evaluates each framework's formula. It does not
#' summarize the shape of the path, and it reports no area or share of the
#' plane: the size of a shaded zone depends on the chosen axis ranges and has
#' no inferential meaning.
#'
#' @param x A `confoundsens` object produced by [sens_path_lm()],
#'   [from_sensemakr()], [itcv_lm()] (its `$path`), [from_konfound()] with
#'   `index = "IT"`, [from_evalue()], or [sens_path_pairs()].
#' @param a,b Optional numeric grids for the two axes. Defaults depend on
#'   the framework and extend beyond the relevant boundary.
#' @param n Number of grid points per axis when `a` or `b` is not supplied.
#' @param alpha Significance level for the partial-\eqn{R^2} boundary.
#'   ITCV and Rosenbaum paths use the level stored in the path (a warning is
#'   given if `alpha` is supplied); for E-values the boundary is defined by
#'   the stored confidence limit.
#'
#' @return An object of class `confoundsurf`: a list with `grid` (a data
#'   frame with columns `a`, `b`, `value`, and `zone`), `boundary` (a data
#'   frame with columns `a`, `b`, `kind`, and `group`), and `meta` (framework,
#'   axis labels, and the quantities used). Plot it with
#'   [plot_sensitivity_surface()].
#'
#' @references
#' Cinelli, C., & Hazlett, C. (2020). Making sense of sensitivity: Extending
#' omitted variable bias. *Journal of the Royal Statistical Society: Series
#' B*, 82(1), 39--67. \doi{10.1111/rssb.12348}
#'
#' Ding, P., & VanderWeele, T. J. (2016). Sensitivity analysis without
#' assumptions. *Epidemiology*, 27(3), 368--377.
#' \doi{10.1097/EDE.0000000000000457}
#'
#' Frank, K. A. (2000). Impact of a confounding variable on a regression
#' coefficient. *Sociological Methods & Research*, 29(2), 147--194.
#' \doi{10.1177/0049124100029002001}
#'
#' Rosenbaum, P. R., & Silber, J. H. (2009). Amplification of sensitivity
#' analysis in matched observational studies. *Journal of the American
#' Statistical Association*, 104(488), 1398--1405.
#' \doi{10.1198/jasa.2009.tm08470}
#' @seealso [plot_sensitivity_surface()]
#' @export
#'
#' @examples
#' fit <- lm(mpg ~ am + wt + hp, data = mtcars)
#' s <- sens_surface(sens_path_lm(fit, "wt"))
#' s
#' plot_sensitivity_surface(s)
#'
#' # E-value: every pair of risk ratios that explains away the estimate
#' if (requireNamespace("EValue", quietly = TRUE)) {
#'   e <- EValue::evalues.RR(1.8, lo = 1.3, hi = 2.5)
#'   plot_sensitivity_surface(sens_surface(from_evalue(e)))
#' }
sens_surface <- function(x, a = NULL, b = NULL, n = 201, alpha = 0.05) {
  if (!inherits(x, "confoundsens")) {
    stop("`x` must be a `confoundsens` object.", call. = FALSE)
  }
  if (!is.null(x$level)) {
    stop("`x` has several levels; pass one path at a time.", call. = FALSE)
  }
  if (!is.numeric(n) || length(n) != 1L || n < 20) {
    stop("`n` must be a single number >= 20.", call. = FALSE)
  }
  n <- as.integer(n)
  fw <- if (is.null(x$meta$framework)) "" else x$meta$framework
  if (fw %in% c("itcv", "rosenbaum") && !missing(alpha)) {
    warning("`alpha` is ignored for ", fw, " paths; the level stored in the ",
            "path is used.", call. = FALSE)
  }
  out <- switch(fw,
    partial_r2 = .surf_pr2(x$meta, a, b, n, alpha),
    itcv       = .surf_itcv(x$meta, a, b, n),
    evalue     = .surf_evalue(x$meta, a, b, n),
    rosenbaum  = .surf_rosenbaum(x, a, b, n),
    rir        = stop("RIR has a single sensitivity parameter; there is no ",
                      "two-parameter surface.", call. = FALSE),
    stop("sens_surface() needs a path from sens_path_lm(), from_sensemakr(), ",
         "itcv_lm(), from_konfound(), from_evalue(), or sens_path_pairs().",
         call. = FALSE))
  out$meta$framework <- fw
  out$meta$treatment <- x$meta$treatment
  if (!is.null(x$meta$bounds)) out$meta$bounds <- x$meta$bounds
  class(out) <- "confoundsurf"
  out
}

.check_grid <- function(v, lo, hi, nm) {
  if (!is.numeric(v) || length(v) < 2L || anyNA(v) || any(v < lo) || any(v > hi)) {
    stop(sprintf("`%s` must be a numeric vector (length >= 2) with values in [%s, %s].",
                 nm, lo, hi), call. = FALSE)
  }
  sort(unique(v))
}

.run_groups <- function(ok) {
  g <- cumsum(c(TRUE, ok[-1L] != ok[-length(ok)]))
  g[!ok] <- NA
  g
}

.empty_bnd <- function() {
  data.frame(a = numeric(0), b = numeric(0), kind = character(0),
             group = character(0), style = character(0), stringsAsFactors = FALSE)
}

# Boundary pieces from a logical mask along a fine a-grid.
.bnd_piece <- function(a, b, ok, kind, tag, style) {
  if (!any(ok)) return(.empty_bnd())
  data.frame(a = a[ok], b = b[ok], kind = kind,
             group = paste0(tag, .run_groups(ok)[ok]), style = style,
             stringsAsFactors = FALSE)
}

# Trace the zero set of a smooth function g(a, b) on [ga] x [gb]: contour
# lines on a fine grid, then each vertex moved onto the curve by root finding
# along the local gradient. Used where the boundary has no closed form.
.trace_zero <- function(g, ga, gb, kind, tag, style, dom = c(0, 1, 0, 1)) {
  z  <- outer(ga, gb, g)
  if (all(is.na(z))) return(.empty_bnd())
  cl <- grDevices::contourLines(ga, gb, z, levels = 0)
  if (!length(cl)) return(.empty_bnd())
  h  <- 1e-7
  step <- 3 * max(diff(ga[1:2]), diff(gb[1:2]))
  lim <- function(v, u, lo, hi) {
    if (u > 0) c((lo - v) / u, (hi - v) / u)
    else if (u < 0) c((hi - v) / u, (lo - v) / u) else c(-Inf, Inf)
  }
  refine <- function(x, y) {
    gx <- (g(x + h, y) - g(x - h, y)) / (2 * h)
    gy <- (g(x, y + h) - g(x, y - h)) / (2 * h)
    nn <- sqrt(gx^2 + gy^2)
    if (!is.finite(nn) || nn == 0) return(c(x, y))
    ux <- gx / nn; uy <- gy / nn
    f1 <- function(t) g(x + t * ux, y + t * uy)
    bx <- lim(x, ux, dom[1] + 1e-9, dom[2] - 1e-9)
    by <- lim(y, uy, dom[3], dom[4] - 1e-9)
    ts <- c(max(-step, bx[1], by[1]), min(step, bx[2], by[2]))
    if (ts[1] >= ts[2]) return(c(x, y))
    f0 <- f1(ts[1]); f2 <- f1(ts[2])
    if (!is.finite(f0) || !is.finite(f2) || sign(f0) == sign(f2)) return(c(x, y))
    t0 <- stats::uniroot(f1, ts, tol = 1e-13)$root
    c(x + t0 * ux, y + t0 * uy)
  }
  out <- lapply(seq_along(cl), function(j) {
    pts <- t(mapply(refine, cl[[j]]$x, cl[[j]]$y))
    data.frame(a = pts[, 1], b = pts[, 2], kind = kind,
               group = paste0(tag, j), style = style, stringsAsFactors = FALSE)
  })
  do.call(rbind, out)
}

.zone_cols <- c("Not significant"                  = "#bdd7e7",
                "Significant, reversed direction"  = "#fcbba1",
                "Estimate at or past zero"         = "#fdd0a2",
                "Significant, original direction"  = "#c7e9c0",
                "Infeasible"                       = "grey60",
                "Confidence interval includes 1"   = "#bdd7e7",
                "Point estimate explained away"    = "#fcbba1",
                "Upper-bound P-value at least 0.5" = "#fcbba1")

# ---- partial R2 ------------------------------------------------------------
.surf_pr2 <- function(m, a, b, n, alpha) {
  est <- m$estimate; se <- m$se0; d <- m$dof
  f2 <- (est / se)^2 / d
  top <- min(0.99, max(3 * m$rv, 0.05))
  if (is.null(a)) a <- seq(0, top, length.out = n) else a <- .check_grid(a, 0, 0.999, "a")
  if (is.null(b)) b <- seq(0, top, length.out = n) else b <- .check_grid(b, 0, 0.999, "b")
  crit <- stats::qt(1 - alpha / 2, d - 1)
  adj <- function(aa, bb) sign(est) * (abs(est) - se * sqrt(d) * sqrt(aa * bb / (1 - aa)))
  tadj <- function(aa, bb) adj(aa, bb) /
    (se * sqrt((1 - bb) / (1 - aa)) * sqrt(d / (d - 1)))
  gr <- expand.grid(a = a, b = b)
  gr$value <- adj(gr$a, gr$b)
  tt <- sign(est) * tadj(gr$a, gr$b)          # t in the original direction
  sig0 <- abs(est / se) * sqrt((d - 1) / d) >= crit
  gr$zone <- ifelse(tt <= -crit, "Significant, reversed direction",
             ifelse(sign(est) * gr$value <= 0, "Estimate at or past zero",
             ifelse(tt < crit, "Not significant",
                    if (sig0) NA_character_ else "Significant, original direction")))

  # Null boundary in closed form: b = f2 (1 - a) / a.
  af <- seq(1e-6, max(a), length.out = 4L * n)
  bn <- f2 * (1 - af) / af
  nul <- .bnd_piece(af, bn, bn <= max(b), "Null (point estimate = 0)", "n", "dashed")
  # Significance boundaries |t| = crit. The adjusted SE shrinks as b
  # approaches 1, so these curves can turn back on themselves.
  gfun <- function(sgn) function(aa, bb) {
    out <- rep(NA_real_, max(length(aa), length(bb)))
    aa <- rep_len(aa, length(out)); bb <- rep_len(bb, length(out))
    ok <- aa > 0 & aa < 1 & bb >= 0 & bb < 1
    out[ok] <- sgn * sign(est) * tadj(aa[ok], bb[ok]) - crit
    out
  }
  K  <- 6L * n
  ga <- seq(1e-6, min(max(a), 0.999), length.out = K)
  gb <- seq(0, min(max(b), 0.999), length.out = K)
  sig <- .trace_zero(gfun(1), ga, gb, sprintf("Significance (alpha = %s)", alpha),
                     "s", "solid")
  rev <- .trace_zero(gfun(-1), ga, gb,
                     sprintf("Significance, reversed direction (alpha = %s)", alpha),
                     "r", "dotted")
  bnd <- rbind(nul, sig, rev)
  list(grid = gr, boundary = bnd,
       meta = list(a_label = "Partial R2 of U with treatment",
                   b_label = "Partial R2 of U with outcome",
                   value_label = "Adjusted estimate",
                   estimate = est, se = se, dof = d, rv = m$rv, alpha = alpha))
}

# ---- ITCV --------------------------------------------------------------------
.surf_itcv <- function(m, a, b, n) {
  r <- m$estimate; rc <- abs(m$r_crit)
  if (is.null(rc) || is.na(rc)) stop("The ITCV path has no critical correlation.", call. = FALSE)
  lim <- min(0.999, max(0.5, 3 * sqrt(max((abs(r) - rc) / (1 - rc), 1e-4))))
  if (is.null(a)) a <- seq(-lim, lim, length.out = n) else a <- .check_grid(a, -0.999, 0.999, "a")
  if (is.null(b)) b <- seq(-lim, lim, length.out = n) else b <- .check_grid(b, -0.999, 0.999, "b")
  gr <- expand.grid(a = a, b = b)
  det <- 1 - r^2 - gr$a^2 - gr$b^2 + 2 * r * gr$a * gr$b
  radj <- (r - gr$a * gr$b) / sqrt((1 - gr$a^2) * (1 - gr$b^2))
  f <- sign(r) * radj
  sig0 <- abs(r) >= rc
  gr$value <- radj
  gr$zone <- ifelse(det <= 0, "Infeasible",
             ifelse(f <= -rc, "Significant, reversed direction",
             ifelse(f < rc, "Not significant",
                    if (sig0) NA_character_ else "Significant, original direction")))
  ex <- .exact_boundary(r, rc, n = 8L * n)   # columns r_yu, r_du
  ex <- ex[ex$r_du >= min(a) & ex$r_du <= max(a) &
             ex$r_yu >= min(b) & ex$r_yu <= max(b), , drop = FALSE]
  bnd <- if (nrow(ex)) data.frame(
    a = ex$r_du, b = ex$r_yu,
    kind = ifelse(ex$kind == "Original direction",
                  "Significance (original direction)",
                  "Significance (reversed direction)"),
    group = ex$group,
    style = ifelse(ex$kind == "Original direction", "solid", "dotted"),
    stringsAsFactors = FALSE) else .empty_bnd()
  list(grid = gr, boundary = bnd,
       meta = list(a_label = "Partial correlation of U with treatment",
                   b_label = "Partial correlation of U with outcome",
                   value_label = "Adjusted partial correlation",
                   r = r, r_crit = rc, itcv = m$itcv))
}

# ---- E-value -------------------------------------------------------------------
.surf_evalue <- function(m, a, b, n) {
  rr <- m$estimate
  lo <- if (is.null(m$lower)) NA_real_ else m$lower
  hi <- if (is.null(m$upper)) NA_real_ else m$upper
  # Work on the scale where the estimate is above 1; the limit closest to 1.
  if (rr < 1) { rr_s <- 1 / rr; lim_s <- if (!is.na(hi)) 1 / hi else NA_real_ }
  else        { rr_s <- rr;     lim_s <- lo }
  ev <- rr_s + sqrt(rr_s * (rr_s - 1))
  top <- max(3, 1.5 * ev)
  if (is.null(a)) a <- seq(1, top, length.out = n) else a <- .check_grid(a, 1, Inf, "a")
  if (is.null(b)) b <- seq(1, top, length.out = n) else b <- .check_grid(b, 1, Inf, "b")
  gr <- expand.grid(a = a, b = b)
  B <- gr$a * gr$b / (gr$a + gr$b - 1)
  adj <- rr_s / B
  gr$value <- if (rr < 1) 1 / adj else adj
  gr$zone <- ifelse(adj <= 1, "Point estimate explained away", NA_character_)
  if (!is.na(lim_s)) {
    if (lim_s > 1) gr$zone[adj > 1 & (lim_s / B) <= 1] <- "Confidence interval includes 1"
    else           gr$zone[adj > 1] <- "Confidence interval includes 1"
  }
  af <- seq(1, max(a), length.out = 4L * n)
  curve <- function(target, kind, tag, style) {
    ok <- af > target
    bb <- rep(NA_real_, length(af))
    bb[ok] <- target * (af[ok] - 1) / (af[ok] - target)
    .bnd_piece(af, bb, ok & bb <= max(b), kind, tag, style)
  }
  bnd <- curve(rr_s, "Null (point estimate)", "n", "dashed")
  if (!is.na(lim_s) && lim_s > 1) {
    bnd <- rbind(bnd, curve(lim_s, "Confidence limit", "s", "solid"))
  }
  list(grid = gr, boundary = bnd,
       meta = list(a_label = "Risk ratio of U with treatment",
                   b_label = "Risk ratio of U with outcome",
                   value_label = "Risk ratio bound",
                   estimate = rr, lower = lo, upper = hi, evalue = ev))
}

# ---- Rosenbaum amplification -------------------------------------------------------
# Gamma at which the upper-bound P-value reaches 0.5 (deviate = 0), by root
# finding on the stored pair differences; falls back to the path grid.
.gamma_half <- function(x) {
  m <- x$meta
  if (is.null(m$d)) return(robustness_points(x)$lambda_null)
  dev <- if (identical(m$method, "binary")) function(g) .p_to_dev(.binary_p(m$d, g))
         else function(g) .wilcox_dev(m$d, g)
  if (dev(1) <= 0) return(1)
  hi <- 2
  while (dev(hi) > 0 && hi < 1e6) hi <- hi * 2
  if (dev(hi) > 0) return(Inf)
  stats::uniroot(dev, c(1, hi), tol = 1e-10)$root
}

.surf_rosenbaum <- function(x, a, b, n) {
  m <- x$meta
  sv <- m$sens_value
  g_half <- .gamma_half(x)
  ref <- if (!is.na(sv) && is.finite(sv)) sv else if (is.finite(g_half)) g_half else max(x$lambda)
  top <- max(4, 4 * ref)
  if (is.null(a)) a <- seq(1, top, length.out = n) else a <- .check_grid(a, 1, Inf, "a")
  if (is.null(b)) b <- seq(1, top, length.out = n) else b <- .check_grid(b, 1, Inf, "b")
  gr <- expand.grid(a = a, b = b)
  G <- (gr$a * gr$b + 1) / (gr$a + gr$b)
  gr$value <- G
  gr$zone <- if (is.na(sv)) "Not significant" else
    ifelse(G >= sv, "Not significant", NA_character_)
  if (!is.na(g_half) && is.finite(g_half)) {
    gr$zone[G >= g_half] <- "Upper-bound P-value at least 0.5"
  }
  af <- seq(1, max(a), length.out = 4L * n)
  curve <- function(target, kind, tag, style) {
    if (is.na(target) || !is.finite(target) || target <= 1) return(.empty_bnd())
    ok <- af > target
    bb <- rep(NA_real_, length(af))
    bb[ok] <- (target * af[ok] - 1) / (af[ok] - target)
    .bnd_piece(af, bb, ok & bb <= max(b), kind, tag, style)
  }
  bnd <- rbind(curve(g_half, "Upper-bound P-value = 0.5", "n", "dashed"),
               curve(sv, sprintf("Significance (one-sided alpha = %s)", m$alpha),
                     "s", "solid"))
  list(grid = gr, boundary = bnd,
       meta = list(a_label = "Lambda (U and treatment)",
                   b_label = "Delta (U and outcome)",
                   value_label = "Implied Gamma",
                   sens_value = sv, gamma_half = g_half, alpha = m$alpha))
}

#' @export
print.confoundsurf <- function(x, ...) {
  m <- x$meta
  cat("<confoundsurf>\n")
  cat("  framework :", m$framework, "\n")
  cat("  a (x-axis):", m$a_label, sprintf("[%g, %g]", min(x$grid$a), max(x$grid$a)), "\n")
  cat("  b (y-axis):", m$b_label, sprintf("[%g, %g]", min(x$grid$b), max(x$grid$b)), "\n")
  cat("  value     :", m$value_label, "\n")
  cat("  boundaries:", paste(unique(x$boundary$kind), collapse = "; "), "\n")
  invisible(x)
}


#' Plot a two-parameter sensitivity surface
#'
#' Draws a [sens_surface()] object: shaded zones where the conclusion
#' changes, the boundaries (solid for loss of significance, dotted for
#' significance in the reversed direction, dashed for the null), and optional
#' benchmark points. Unshaded cells are those where the original conclusion
#' holds.
#'
#' @param x A `confoundsurf` object.
#' @param benchmarks Optional data frame of reference points with columns
#'   `a` and `b` (and optionally `label`) on the axes of the surface. For a
#'   partial-\eqn{R^2} surface, the output of
#'   `covariate_impacts(..., metric = "partial_r2")` is accepted (columns
#'   `r2dz.x` and `r2yz.dx`); for an ITCV surface, the output of
#'   `covariate_impacts(..., metric = "itcv")` (columns `r_du` and `r_yu`).
#'   For a surface built from a `sensemakr` path, its formal bounds are
#'   added automatically unless `bounds = FALSE`.
#' @param bounds Logical; draw `sensemakr` formal bounds when present.
#'
#' @return A [ggplot2::ggplot()] object.
#' @seealso [sens_surface()]
#' @export
#'
#' @examples
#' fit <- lm(mpg ~ am + wt + hp, data = mtcars)
#' s <- sens_surface(sens_path_lm(fit, "wt"))
#' plot_sensitivity_surface(s, benchmarks = covariate_impacts(fit, "wt",
#'                                                            "partial_r2"))
plot_sensitivity_surface <- function(x, benchmarks = NULL, bounds = TRUE) {
  if (!inherits(x, "confoundsurf")) {
    stop("`x` must be a `confoundsurf` object from sens_surface().", call. = FALSE)
  }
  m <- x$meta
  gr <- x$grid
  present <- unique(stats::na.omit(gr$zone))
  zones <- .zone_cols[names(.zone_cols) %in% present]
  p <- ggplot2::ggplot() +
    ggplot2::geom_raster(
      data = gr, ggplot2::aes(x = .data$a, y = .data$b, fill = .data$zone),
      na.rm = TRUE, interpolate = FALSE) +
    ggplot2::scale_fill_manual(values = zones, breaks = names(zones),
                               na.value = NA, name = NULL, drop = FALSE)
  if (nrow(x$boundary)) {
    bd <- x$boundary
    lt <- stats::setNames(bd$style, bd$kind)
    lt <- lt[!duplicated(names(lt))]
    p <- p + ggplot2::geom_path(
      data = bd,
      ggplot2::aes(x = .data$a, y = .data$b, group = .data$group,
                   linetype = .data$kind), linewidth = 0.7) +
      ggplot2::scale_linetype_manual(values = lt, name = "Boundary")
  }
  pts <- NULL
  if (!is.null(benchmarks)) pts <- .surface_points(benchmarks, m$framework)
  if (isTRUE(bounds) && !is.null(m$bounds) && m$framework == "partial_r2") {
    bb <- m$bounds
    pts <- rbind(pts, data.frame(a = bb$r2dz.x, b = bb$r2yz.dx,
                                 label = bb$bound_label, stringsAsFactors = FALSE))
  }
  if (!is.null(pts) && nrow(pts)) {
    p <- p + ggplot2::geom_point(data = pts,
                                 ggplot2::aes(x = .data$a, y = .data$b), size = 2.3)
    if ("label" %in% names(pts)) {
      if (requireNamespace("ggrepel", quietly = TRUE)) {
        p <- p + ggrepel::geom_text_repel(
          data = pts, ggplot2::aes(x = .data$a, y = .data$b, label = .data$label),
          size = 3.3, min.segment.length = 0, seed = 1)
      } else {
        p <- p + ggplot2::geom_text(
          data = pts, ggplot2::aes(x = .data$a, y = .data$b, label = .data$label),
          size = 3.3, nudge_y = 0.02 * diff(range(gr$b)))
      }
    }
  }
  p + ggplot2::coord_cartesian(xlim = range(gr$a), ylim = range(gr$b), expand = FALSE) +
    ggplot2::labs(x = m$a_label, y = m$b_label)
}

.surface_points <- function(bm, fw) {
  if (!is.data.frame(bm)) stop("`benchmarks` must be a data frame.", call. = FALSE)
  lab <- if ("label" %in% names(bm)) bm$label else if ("covariate" %in% names(bm)) bm$covariate else NULL
  if (all(c("a", "b") %in% names(bm))) {
    out <- data.frame(a = bm$a, b = bm$b)
  } else if (fw == "partial_r2" && all(c("r2dz.x", "r2yz.dx") %in% names(bm))) {
    out <- data.frame(a = bm$r2dz.x, b = bm$r2yz.dx)
  } else if (fw == "itcv" && all(c("r_du", "r_yu") %in% names(bm))) {
    out <- data.frame(a = bm$r_du, b = bm$r_yu)
  } else {
    stop("`benchmarks` needs columns `a` and `b` on the axes of the surface.",
         call. = FALSE)
  }
  if (!is.null(lab)) out$label <- as.character(lab)
  out[stats::complete.cases(out[c("a", "b")]), , drop = FALSE]
}
