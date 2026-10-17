#' confoundvis: A Common Interface for Visualizing and Reporting Sensitivity Analyses for Unmeasured Confounding
#'
#' @description
#' A common interface to established sensitivity analyses for unmeasured
#' confounding. Results from the impact threshold (ITCV) and robustness of
#' inference to replacement (RIR) frameworks ('konfound'), partial
#' R-squared omitted-variable bias analysis ('sensemakr'), E-values
#' ('EValue'), and Rosenbaum bounds for matched pairs ('DOS2',
#' 'sensitivitymv') are stored in one `confoundsens` object, so the same
#' plots, tables, and plain-language reports apply to all of them. Each
#' result keeps the parameter and the published formula of its own
#' framework; confoundvis adds no new sensitivity parameter or measure.
#'
#' @details
#' Main entry points:
#' \itemize{
#'   \item Paths from fitted models: [sens_path_lm()], [itcv_lm()],
#'     [rir_lm()].
#'   \item Paths from other packages: [from_sensemakr()], [from_konfound()],
#'     [from_evalue()], or [as_confoundsens()] for data frames of
#'     precomputed (including subgroup-specific) paths.
#'   \item Paths for matched pairs (Rosenbaum's sensitivity model):
#'     [sens_path_pairs()].
#'   \item Covariate benchmarks: [covariate_impacts()].
#'   \item Two-parameter surfaces: [sens_surface()],
#'     [plot_sensitivity_surface()].
#'   \item Plots: [plot_robustness_curve()], [plot_sensitivity_love()],
#'     [plot_sensitivity_contour()].
#'   \item Reporting: [robustness_points()], [sens_report()], and
#'     [sens_table()] for several paths side by side.
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
