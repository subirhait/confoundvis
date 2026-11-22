#' Sensitivity analysis for unmeasured confounding in one call
#'
#' The simplest way to use the package. Give it a fitted model and the name
#' of the treatment, or a vector of matched-pair differences, and it runs
#' every analysis the package supports for that input, compares the result
#' with the measured covariates, and prints a short summary in plain
#' language. The returned object can be printed, plotted, and passed to the
#' other functions of the package.
#'
#' What is run depends on the input:
#' \itemize{
#'   \item A linear model fitted by [stats::lm()]: the partial-\eqn{R^2},
#'     ITCV, and RIR analyses ([sens_profile()]) and the measured-covariate
#'     benchmarks on both regression metrics ([covariate_impacts()]).
#'   \item A logistic regression or other [stats::glm()] fit: the RIR
#'     computed by 'konfound' (which must be installed) and imported with
#'     [from_konfound()]. konfound's RIR for such models is an
#'     approximation; see its documentation.
#'   \item A numeric vector: treated as treated-minus-control differences in
#'     matched pairs, analysed with Rosenbaum's sensitivity model
#'     ([sens_path_pairs()]; one-sided test).
#'   \item An instrumental-variable analysis from `iv.sensemakr`
#'     (an `iv_fit()` or `sensemakr()` result): the sensitivity of the
#'     Anderson--Rubin inference ([from_iv_sensemakr()]).
#'   \item A `DoubleMLPLR` or `DoubleMLIRM` fit from `DoubleML`, fitted with
#'     `store_predictions = TRUE`: the omitted-variable-bias analysis of
#'     Chernozhukov et al. ([sens_dml()]), with `rho = 1` and the one-sided
#'     level `1 - alpha / 2`, so that its significance crossing matches a
#'     two-sided test at `alpha`.
#' }
#'
#' Each analysis keeps its own sensitivity parameter and decision rule. The
#' numbers measure different things and must not be ranked against one
#' another. A sensitivity analysis shows how strong an unmeasured
#' confounder would need to be; it does not show whether one exists.
#'
#' @param x A fitted `lm` or `glm` object, a numeric vector of matched-pair
#'   differences, an `iv.sensemakr` `iv_fit` or `sensemakr` result, or a
#'   fitted `DoubleMLPLR` or `DoubleMLIRM` object.
#' @param treatment Name of the treatment coefficient (for a model).
#' @param alpha Significance level (two-sided for models, one-sided for
#'   matched pairs). For an `iv.sensemakr` `sensemakr()` result, the level of
#'   that object is used unless `alpha` is given, in which case the analysis
#'   is rerun at `alpha`.
#' @param labels Optional named character vector of readable covariate
#'   labels, passed to [covariate_impacts()].
#' @param extra Optional named list of further `confoundsens` paths (for
#'   example from [from_evalue()]) to add for a linear model.
#'
#' @return An object of class `confoundvis_result`: a list with `type`
#'   (`"lm"`, `"glm"`, `"pairs"`, `"iv"`, or `"dml"`), `treatment`, `paths` (named list of
#'   `confoundsens` objects), `table` ([sens_table()] output), `benchmarks`
#'   (for `"lm"`: the [covariate_impacts()] results for the two metrics),
#'   `profile` (for `"lm"`), and `summary` (the plain-language lines).
#' @seealso [sens_profile()], [sens_report()], [covariate_impacts()]
#' @export
#'
#' @examples
#' fit <- lm(mpg ~ am + wt + hp, data = mtcars)
#' res <- confoundvis(fit, "wt")
#' res
#' plot(res)
#' plot(res, which = "benchmarks")
#'
#' # matched pairs (Rosenbaum bounds)
#' d <- c(1.2, 0.8, 2.1, 0.4, 1.7, -0.3, 0.9, 1.5, 2.2, 0.6, 1.1, 0.7)
#' confoundvis(d)
confoundvis <- function(x, treatment = NULL, alpha = 0.05, labels = NULL,
                        extra = NULL) {
  if (is.numeric(x) && is.null(dim(x)) && !inherits(x, "lm")) {
    p <- sens_path_pairs(x, gamma = seq(1, 20, by = 0.01), alpha = alpha)
    paths <- list(rosenbaum = p)
    out <- list(type = "pairs", treatment = NA_character_, paths = paths,
                table = sens_table(paths), benchmarks = NULL, profile = NULL)
  } else if (inherits(x, c("iv.sensemakr", "iv_fit"))) {
    if (!requireNamespace("iv.sensemakr", quietly = TRUE)) {
      stop("Package 'iv.sensemakr' is needed for this input.", call. = FALSE)
    }
    if (inherits(x, "iv_fit")) {
      x <- iv.sensemakr::sensemakr(x, alpha = alpha)
    } else if (!missing(alpha) && abs(alpha - x$pars$alpha) > 1e-12) {
      pr <- x$pars
      x <- iv.sensemakr::sensemakr(x$unadjusted,
                                   benchmark_covariates = pr$benchmark_covariates,
                                   kz = pr$kz, ky = pr$ky, kd = pr$kd,
                                   q = pr$q, alpha = alpha)
    }
    alpha <- x$pars$alpha
    paths <- list(iv = from_iv_sensemakr(x))
    out <- list(type = "iv", treatment = paths$iv$meta$treatment, paths = paths,
                table = sens_table(paths), benchmarks = NULL, profile = NULL)
  } else if (inherits(x, "DoubleML")) {
    paths <- list(dml = sens_dml(x, level = 1 - alpha / 2,
                                 treatment = treatment))
    out <- list(type = "dml", treatment = paths$dml$meta$treatment, paths = paths,
                table = sens_table(paths), benchmarks = NULL, profile = NULL)
  } else if (inherits(x, "glm")) {
    if (is.null(treatment)) stop("Give the name of the treatment coefficient.", call. = FALSE)
    if (!requireNamespace("konfound", quietly = TRUE)) {
      stop("For a glm fit, confoundvis() uses the 'konfound' package. ",
           "Install it with install.packages(\"konfound\").", call. = FALSE)
    }
    if (!treatment %in% names(stats::coef(x))) {
      stop("`treatment` is not a coefficient of the model.", call. = FALSE)
    }
    k <- suppressMessages(eval(bquote(konfound::konfound(
      x, .(as.name(treatment)), alpha = .(alpha), index = "RIR",
      to_return = "raw_output"))))
    paths <- list(rir = from_konfound(k, index = "RIR", treatment = treatment,
                                      alpha = alpha))
    out <- list(type = "glm", treatment = treatment, paths = paths,
                table = sens_table(paths), benchmarks = NULL, profile = NULL)
  } else if (inherits(x, "lm")) {
    if (is.null(treatment)) stop("Give the name of the treatment coefficient.", call. = FALSE)
    prof <- sens_profile(x, treatment, alpha = alpha, extra = extra)
    bench <- list(
      itcv = covariate_impacts(x, treatment, metric = "itcv", alpha = alpha,
                               labels = labels),
      partial_r2 = covariate_impacts(x, treatment, metric = "partial_r2",
                                     alpha = alpha, labels = labels))
    out <- list(type = "lm", treatment = treatment, paths = prof$paths,
                table = prof$table, benchmarks = bench, profile = prof)
  } else {
    stop("`x` must be a fitted lm or glm model, a numeric vector of ",
         "matched-pair differences, an iv.sensemakr result, or a fitted ",
         "DoubleMLPLR or DoubleMLIRM object.", call. = FALSE)
  }
  out$alpha <- alpha
  out$summary <- .cv_summary(out)
  class(out) <- "confoundvis_result"
  out
}

# Internal: plain-language summary lines
.cv_summary <- function(r) {
  f <- function(v, d = 3) formatC(signif(v, d), format = "fg", digits = d)
  pct <- function(v) paste0(formatC(100 * v, format = "f", digits = 1), "%")
  tb <- r$table
  row <- function(nm) tb[tb$path == nm, , drop = FALSE][1, ]
  lines <- character()
  for (nm in names(r$paths)) {
    m <- r$paths[[nm]]$meta
    rw <- row(nm)
    fw <- if (is.null(m$framework)) nm else m$framework
    txt <- switch(fw,
      partial_r2 = sprintf(paste0(
        "Partial R2: an omitted confounder explaining %s of the residual ",
        "variance of both treatment and outcome would remove the estimate ",
        "(robustness value); %s would make it nonsignificant."),
        pct(rw$lambda_null), pct(rw$lambda_sig)),
      itcv = sprintf(paste0(
        "ITCV: a confounder whose partial correlations with treatment and ",
        "outcome each reach about %s (product %s) would make the estimate ",
        "nonsignificant."), f(sqrt(rw$lambda_sig)), f(rw$lambda_sig)),
      rir = sprintf(paste0(
        "RIR: %s%% of the estimate would have to be due to bias; %s cases ",
        "would have to be replaced with cases with no effect."),
        formatC(m$perc_bias, format = "f", digits = 1),
        if (is.null(m$rir) || is.na(m$rir)) "the number of" else
          format(m$rir, big.mark = ",")),
      evalue = sprintf(paste0(
        "E-value: a confounder with a risk ratio of %s with both treatment ",
        "and outcome could explain away the estimate; %s would move the ",
        "confidence limit to the null."), f(rw$lambda_null), f(rw$lambda_sig)),
      iv = sprintf(paste0(
        "Instrumental variables (Anderson-Rubin): an omitted variable explaining ",
        "%s of the residual variance of both the instrument and the outcome ",
        "would make the estimate nonsignificant (robustness value); %s would ",
        "move the bound on the estimate itself to the null."),
        pct(m$rv), pct(m$rv_point)),
      dml = sprintf(paste0(
        "DoubleML (rho = %s): a confounder with cf_y = cf_d = %s would move the ",
        "bias-adjusted bound to the null (RV); %s would make it ",
        "nonsignificant at alpha = %s (RVa)."),
        f(m$rho), pct(m$rv), pct(rw$lambda_sig), f(rw$alpha)),
      rosenbaum = sprintf(paste0(
        "Rosenbaum: the result stays significant (one-sided alpha = %s) ",
        "unless a hidden factor could make one unit in a pair %s times as ",
        "likely as the other, in odds, to be treated (Gamma)."),
        f(rw$alpha), f(m$sens_value)),
      sprintf("%s: significance is lost at lambda = %s.", nm, f(rw$lambda_sig)))
    lines <- c(lines, txt)
  }
  if (!is.null(r$benchmarks)) {
    for (mt in names(r$benchmarks)) {
      b <- r$benchmarks[[mt]]
      thr <- attr(b, "threshold")
      ok <- is.finite(b$impact)
      if (!any(ok)) next
      top <- b[ok, ][which.max(b$impact[ok]), ]
      n_over <- sum(b$impact[ok] >= thr)
      lines <- c(lines, sprintf(paste0(
        "Benchmark (%s): the strongest measured covariate is %s ",
        "(%s vs. threshold %s); %s."),
        if (mt == "itcv") "ITCV" else "partial R2", top$label,
        f(top$impact), f(thr),
        if (n_over == 0) "no measured covariate reaches the threshold" else
          sprintf("%d measured covariate%s reach the threshold", n_over,
                  if (n_over == 1) "" else "s")))
    }
  }
  c(lines, if (length(r$paths) > 1) paste0(
    "These numbers use different scales; do not compare or rank them. They ",
    "show how strong a confounder would need to be, not whether one exists.")
    else paste0("This value shows how strong a confounder would need to be, ",
                "not whether one exists."))
}

#' @export
print.confoundvis_result <- function(x, ...) {
  cat("confoundvis: sensitivity to unmeasured confounding\n")
  m1 <- x$paths[[1]]$meta
  if (x$type == "pairs") {
    cat(sprintf("Matched pairs: %s\n\n", m1$n_pairs))
  } else {
    est <- if (!is.null(m1$estimate)) m1$estimate else x$paths[[1]]$theta[1]
    cat(sprintf("Treatment: %s; estimate %s%s\n\n", x$treatment,
                formatC(signif(est, 4), format = "fg", digits = 4),
                if (x$type == "glm") " (log-odds scale)" else ""))
  }
  keep <- c("path", "summary_name", "summary_value", "lambda_null",
            "lambda_sig", "alpha")
  tb <- x$table[, keep]
  class(tb) <- "data.frame"
  for (v in c("summary_value", "lambda_null", "lambda_sig")) tb[[v]] <- signif(tb[[v]], 4)
  print(tb, row.names = FALSE)
  cat("\nIn plain language:\n")
  for (l in x$summary) {
    cat(paste(strwrap(l, width = 76, initial = "- ", prefix = "  "),
              collapse = "\n"), "\n", sep = "")
  }
  pl <- "Plots: plot(x) draws the paths"
  if (!is.null(x$benchmarks)) pl <- paste0(pl, "; plot(x, which = \"benchmarks\") the covariates")
  if (x$type != "glm") pl <- paste0(pl, "; plot(x, which = \"surface\") the surface")
  cat("\n", paste(strwrap(paste0(pl, "."), width = 78), collapse = "\n"), "\n", sep = "")
  invisible(x)
}

#' Plot a confoundvis result
#'
#' @param x A `confoundvis_result` from [confoundvis()].
#' @param which `"paths"` (default) for the sensitivity paths, `"benchmarks"`
#'   for the measured-covariate Love plot (partial-\eqn{R^2} metric),
#'   `"benchmarks_itcv"` for the same on the ITCV metric, or `"surface"` for
#'   the two-parameter surface of the first analysis that has one.
#' @param top Number of covariates shown in the benchmark plots.
#' @param ... Unused.
#' @return A [ggplot2::ggplot()] object.
#' @export
plot.confoundvis_result <- function(x, which = c("paths", "benchmarks",
                                                 "benchmarks_itcv", "surface"),
                                    top = 15, ...) {
  which <- match.arg(which)
  if (which == "paths") {
    if (!is.null(x$profile)) return(plot(x$profile))
    return(plot_robustness_curve(x$paths[[1]], points = FALSE))
  }
  if (which %in% c("benchmarks", "benchmarks_itcv")) {
    if (is.null(x$benchmarks)) {
      stop("Covariate benchmarks are available for linear models only.", call. = FALSE)
    }
    b <- x$benchmarks[[if (which == "benchmarks") "partial_r2" else "itcv"]]
    return(plot_sensitivity_love(b, top = top))
  }
  has_surface <- vapply(x$paths, function(p) !identical(p$meta$framework, "rir"),
                        logical(1))
  if (!any(has_surface)) {
    stop("RIR has one parameter and no two-parameter surface.", call. = FALSE)
  }
  p <- x$paths[[which(has_surface)[1]]]
  bm <- if (!is.null(x$benchmarks) && identical(p$meta$framework, "partial_r2")) {
    b <- x$benchmarks$partial_r2
    utils::head(b[order(-b$impact), ], 3)
  } else NULL
  plot_sensitivity_surface(sens_surface(p), benchmarks = bm)
}
