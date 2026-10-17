#' confoundvis: Visualization and Reporting Tools for Sensitivity Analysis
#'
#' @description
#' Draws and reports sensitivity analyses for unmeasured confounding. A
#' common `confoundsens` object stores a sensitivity path (the treatment
#' estimate as a function of hypothetical confounder strength) regardless of
#' the framework that produced it, so the same displays and reports work for
#' impact threshold (ITCV) analysis, partial R-squared omitted-variable bias
#' analysis, robustness of inference to replacement (RIR), E-values, and
#' Rosenbaum bounds for matched pairs.
#'
#' @details
#' Main entry points:
#' \itemize{
#'   \item Paths from fitted models: [sens_path_lm()], [itcv_lm()],
#'     [rir_lm()].
#'   \item Paths from other packages: [from_sensemakr()], [from_konfound()],
#'     [from_evalue()], or [as_confoundsens()] for data frames of
#'     precomputed (including stratified or multilevel) paths.
#'   \item Paths for matched pairs (Rosenbaum's sensitivity model):
#'     [sens_path_pairs()].
#'   \item Covariate benchmarks: [covariate_impacts()].
#'   \item Plots: [plot_robustness_curve()], [plot_sensitivity_love()],
#'     [plot_sensitivity_contour()], [plot_local_taylor()].
#'   \item Reporting: [robustness_points()], [sens_report()], and
#'     [sens_table()] for several paths side by side.
#'   \item Conceptual illustrations built on stylized models:
#'     [plot_reversal_cone()], [plot_taylor_panels()].
#' }
#'
#' A sensitivity display shows how strong confounding would have to be to
#' change a conclusion. It cannot show whether such a confounder exists, and
#' it inherits the assumptions of the framework that produced the path.
#'
#' @keywords internal
#'
#' @importFrom rlang .data
#' @importFrom ggplot2 ggplot aes geom_line geom_point geom_ribbon geom_vline geom_path geom_text coord_equal facet_wrap labs
#' @importFrom stats lm coef
#' @importFrom graphics par
"_PACKAGE"
