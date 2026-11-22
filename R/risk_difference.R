# ---------------------------------------------------------------------------
# Sensitivity of a standardized risk difference to a binary unmeasured
# confounder (Lin, Psaty, and Kronmal, 1998; VanderWeele and Arah, 2011).
#
# Model. U is binary with P(U = 1 | A = a, X) = p_a (it does not depend on X
# given A, the assumption of Lin et al.), and U acts on the outcome either
#   odds_ratio:      logit P(Y = 1 | a, X, U) = h_a(X) + log(gamma) U
#   risk_difference: P(Y = 1 | a, X, U)       = h_a(X) + delta U
# with no U-by-treatment interaction. The fitted model gives the observed
# risk r_a(X) = P(Y = 1 | A = a, X), which mixes the two U strata with
# weights p_a. Solving for h_a(X) and re-weighting the strata by the U
# distribution of the target population gives the counterfactual risks.
# For the odds-ratio model the solution is the positive root of a quadratic
# in t = exp(h): G (1 - r) t^2 + (p G + 1 - p - r (1 + G)) t - r = 0.
# ---------------------------------------------------------------------------

# exp(h) solving p * G t / (1 + G t) + (1 - p) * t / (1 + t) = r
.rd_t <- function(r, p, G) {
  A <- G * (1 - r)
  B <- p * G + 1 - p - r * (1 + G)
  D <- sqrt(B^2 + 4 * A * r)
  ifelse(B >= 0, 2 * r / (B + D), (-B + D) / (2 * A))
}

# Counterfactual risk and its derivative with respect to r, for observed
# risk r in arm a (U prevalence p), target U prevalence q, and strength G
# (odds ratio) or delta (risk difference).
.rd_cf <- function(r, p, q, s, type) {
  if (type == "odds_ratio") {
    if (abs(s - 1) < 1e-15) return(list(P = r, dP = rep(1, length(r))))
    t <- .rd_t(r, p, s)
    e0 <- t / (1 + t); e1 <- s * t / (1 + s * t)
    d0 <- t / (1 + t)^2; d1 <- s * t / (1 + s * t)^2
    P <- q * e1 + (1 - q) * e0
    dP <- (q * d1 + (1 - q) * d0) / (p * d1 + (1 - p) * d0)
    list(P = P, dP = dP)
  } else {
    # stratum risks h (U = 0) and h + delta (U = 1) must be probabilities
    h <- r - s * p
    if (any(h < -1e-12 | h > 1 + 1e-12 | h + s < -1e-12 | h + s > 1 + 1e-12)) {
      return(list(P = rep(NA_real_, length(r)), dP = rep(NA_real_, length(r))))
    }
    list(P = r + s * (q - p), dP = rep(1, length(r)))
  }
}

# Largest additive confounder effect (in the direction of the path) for
# which every stratum risk lies in [0, 1].
.rd_max_delta <- function(ob, dir) {
  lim <- function(r, p) {
    # need r - s p in [0, 1] and r + s (1 - p) in [0, 1], s = dir * delta
    if (dir > 0) min(c(r / p, (1 - r) / (1 - p))[c(p > 0, p < 1)])
    else min(c((1 - r) / p, r / (1 - p))[c(p > 0, p < 1)])
  }
  min(lim(ob$r1, ob$p1), lim(ob$r0, ob$p0))
}

# Standardized risks, risk difference, and delta-method SE at strengths s
# (and, for surfaces, treated-group prevalences p1v), in chunks so that the
# gradients are matrix products.
.rd_eval <- function(ob, s, p1v = NULL, se = TRUE) {
  k <- length(s)
  if (is.null(p1v)) p1v <- rep(ob$p1, k)
  out <- matrix(NA_real_, k, 4L,
                dimnames = list(NULL, c("risk1", "risk0", "rd", "se")))
  for (ch in split(seq_len(k), ceiling(seq_len(k) / 400L))) {
    m <- length(ch)
    G1 <- G0 <- matrix(0, ob$n, m)
    for (j in seq_len(m)) {
      i <- ch[j]
      q <- if (ob$target == "ATE") p1v[i] * ob$propensity + ob$p0 * (1 - ob$propensity)
           else rep(p1v[i], ob$n)
      c1 <- .rd_cf(ob$r1, p1v[i], q, s[i], ob$type)
      c0 <- .rd_cf(ob$r0, ob$p0, q, s[i], ob$type)
      out[i, 1:2] <- c(sum(ob$w * c1$P), sum(ob$w * c0$P))
      G1[, j] <- ob$w * c1$dP * ob$mu1
      G0[, j] <- ob$w * c0$dP * ob$mu0
    }
    if (se) {
      g <- crossprod(ob$X1, G1) - crossprod(ob$X0, G0)
      out[ch, 4] <- sqrt(pmax(colSums(g * (ob$V %*% g)), 0))
    }
  }
  out[, 3] <- out[, 1] - out[, 2]
  as.data.frame(out)
}

# Strength at which the confidence limit nearest zero reaches `thr`.
.rd_ci_cross <- function(ob, to_s, lo, hi, s0, thr, alpha) {
  z <- stats::qnorm(1 - alpha / 2)
  g <- function(l) { e <- .rd_eval(ob, to_s(l)); s0 * (e$rd - thr) - z * e$se }
  g0 <- g(lo)
  if (!is.finite(g0)) return(NA_real_)
  if (g0 <= 0) return(lo)
  gh <- g(hi)
  if (!is.finite(gh) || gh > 0) return(NA_real_)
  stats::uniroot(g, c(lo, hi), tol = 1e-10)$root
}

# Threshold crossing of the confidence limit at level alpha, from a path.
.rd_ci_cross_meta <- function(m, alpha) {
  type <- m$confounder_effect; dir <- m$direction
  to_s <- function(l) if (type == "odds_ratio") l^dir else dir * l
  lo <- if (type == "odds_ratio") 1 else 0
  hi <- if (type == "odds_ratio") 1e4 else min(1 - 1e-9, m$max_feasible)
  .rd_ci_cross(m$rd_setup, to_s, lo, hi, if (m$estimate >= 0) 1 else -1,
               m$rd_threshold, alpha)
}

# Build everything that does not depend on the confounder strength.
.rd_setup <- function(model, treatment, p1, p0, target, propensity, type) {
  if (!inherits(model, "glm") || stats::family(model)$family != "binomial") {
    stop("`model` must be a binomial glm (for example logistic regression).",
         call. = FALSE)
  }
  pw <- stats::weights(model, type = "prior")
  if (!is.null(pw) && any(abs(pw - 1) > 1e-12)) {
    stop("The outcome must be binary (0/1) with unit weights.", call. = FALSE)
  }
  dat <- model$data
  if (is.null(dat) || !is.data.frame(dat)) {
    stop("Fit the model with a `data` argument so that risks can be ",
         "predicted with the treatment set to each level.", call. = FALSE)
  }
  if (!is.character(treatment) || length(treatment) != 1L ||
      !treatment %in% names(dat)) {
    stop("`treatment` must be the name of a variable in the model data.",
         call. = FALSE)
  }
  mf <- stats::model.frame(model)
  rows <- rownames(mf)
  if (!all(rows %in% rownames(dat))) {
    stop("Rows of the model frame could not be matched to the model data.",
         call. = FALSE)
  }
  dat <- dat[rows, , drop = FALSE]
  a <- dat[[treatment]]
  if (is.factor(a) || is.character(a)) {
    lv <- if (is.factor(a)) levels(a) else sort(unique(a))
    if (length(lv) != 2L) stop("The treatment must have exactly two levels.", call. = FALSE)
    a1 <- as.integer(a == lv[2L])
    set <- function(v) if (is.factor(a)) factor(lv[v + 1L], levels = lv) else lv[v + 1L]
  } else if (is.logical(a)) {
    a1 <- as.integer(a); set <- function(v) as.logical(v)
  } else {
    if (!all(a %in% c(0, 1))) stop("A numeric treatment must be coded 0/1.", call. = FALSE)
    a1 <- as.integer(a); set <- function(v) v
  }
  if (!treatment %in% all.vars(stats::formula(model)[[3L]])) {
    stop("`treatment` does not appear on the right-hand side of the model.",
         call. = FALSE)
  }
  for (nm in c("p1", "p0")) {
    v <- get(nm)
    if (!is.numeric(v) || length(v) != 1L || is.na(v) || v < 0 || v > 1) {
      stop(sprintf("`%s` must be a single proportion in [0, 1].", nm), call. = FALSE)
    }
  }
  if (p1 == p0) {
    stop("`p1` equals `p0`: a confounder equally common in both groups ",
         "cannot bias the comparison in this model.", call. = FALSE)
  }

  # Observed risks with the treatment set to each level, and the model
  # matrices for the delta method.
  d1 <- dat; d1[[treatment]] <- set(rep(1L, nrow(dat)))
  d0 <- dat; d0[[treatment]] <- set(rep(0L, nrow(dat)))
  tt <- stats::delete.response(stats::terms(model))
  mm <- function(d) stats::model.matrix(tt, stats::model.frame(tt, d, xlev = model$xlevels),
                                        contrasts.arg = model$contrasts)
  X1 <- mm(d1); X0 <- mm(d0)
  b <- stats::coef(model)
  ok <- !is.na(b)
  X1 <- X1[, ok, drop = FALSE]; X0 <- X0[, ok, drop = FALSE]; b <- b[ok]
  fam <- stats::family(model)
  # linear predictors from predict(), so that offsets are kept
  eta1 <- unname(stats::predict(model, newdata = d1, type = "link"))
  eta0 <- unname(stats::predict(model, newdata = d0, type = "link"))
  clamp <- function(r) pmin(pmax(r, 1e-12), 1 - 1e-12)
  r1 <- clamp(fam$linkinv(eta1)); r0 <- clamp(fam$linkinv(eta0))
  V <- stats::vcov(model)[ok, ok, drop = FALSE]

  # Distribution of U in the target population.
  n <- nrow(dat)
  if (target == "ATE") {
    if (is.null(propensity)) {
      keep <- attr(tt, "term.labels")
      keep <- keep[!vapply(keep, function(z)
        treatment %in% all.vars(stats::reformulate(z)), logical(1))]
      pdat <- dat; pdat$.cv_a <- a1
      fpr <- if (length(keep)) stats::reformulate(keep, response = ".cv_a")
             else stats::as.formula(".cv_a ~ 1")
      propensity <- suppressWarnings(stats::fitted(
        stats::glm(fpr, data = pdat, family = stats::binomial())))
    } else if (!is.numeric(propensity) || length(propensity) != n ||
               anyNA(propensity) || any(propensity < 0 | propensity > 1)) {
      stop(sprintf("`propensity` must be a vector of %d probabilities.", n),
           call. = FALSE)
    }
    q <- p1 * propensity + p0 * (1 - propensity)
    w <- rep(1 / n, n)
  } else {
    q <- rep(p1, n)
    w <- a1 / sum(a1)
  }
  list(r1 = r1, r0 = r0, mu1 = fam$mu.eta(eta1), mu0 = fam$mu.eta(eta0),
       X1 = X1, X0 = X0, V = V, w = w, q = q, p1 = p1, p0 = p0, type = type,
       n = n, n_treated = sum(a1), target = target,
       propensity = if (target == "ATE") propensity else NULL)
}

#' Sensitivity of a risk difference to a binary unmeasured confounder
#'
#' For a binary outcome and a binary treatment analysed with a logistic (or
#' other binomial) regression, computes the standardized risk difference
#' that would be obtained after adjusting for a hypothetical binary
#' confounder \eqn{U}, over a range of confounder strengths. The confounder
#' is described as in Lin, Psaty, and Kronmal (1998): its prevalence in the
#' treated (`p1`) and untreated (`p0`) groups, conditional on the measured
#' covariates, and its association with the outcome, either an odds ratio
#' (the default) or a risk difference (VanderWeele and Arah, 2011). The
#' answer is on the absolute-risk scale, so it can be compared with a
#' practically important difference given in `threshold`.
#'
#' The fitted model gives each person's risk with the treatment set to each
#' level. Under the confounder model these risks are mixtures over the two
#' strata of \eqn{U}; the function solves for the stratum risks (in closed
#' form) and averages them over the distribution of \eqn{U} in the target
#' population (g-computation). With no confounding the result is the
#' usual standardized risk difference, the same as
#' `marginaleffects::avg_comparisons()` for the treatment. Assumptions: the
#' prevalence of \eqn{U} depends on the treatment but not on the measured
#' covariates, \eqn{U} does not modify the effect of treatment, and for the
#' odds-ratio model the fitted model is the outcome model in each stratum of
#' \eqn{U} up to the term for \eqn{U}. For the risk-difference model the
#' adjustment is exact and equals \eqn{RD - \delta (p_1 - p_0)}.
#'
#' The sensitivity parameter \eqn{\lambda} is the strength of the
#' confounder-outcome association, oriented so that larger values move the
#' estimate toward zero: for the odds-ratio model the odds ratio is
#' \eqn{\lambda} or \eqn{1/\lambda}, and for the risk-difference model the
#' risk difference is \eqn{\pm\lambda}, with the sign chosen from the signs
#' of the estimate and of \eqn{p_1 - p_0}.
#'
#' Two kinds of uncertainty are kept apart. The path over \eqn{\lambda}
#' describes the range across confounding scenarios; it is not sampling
#' error. Standard errors (delta method over the model coefficients, with
#' the covariates held fixed, as in `marginaleffects`) describe sampling
#' error for each scenario. [robustness_points()] uses them for the
#' significance crossing, and `meta` also gives the strength at which the
#' estimate reaches `threshold` and at which the confidence limit does.
#'
#' @param model A fitted [stats::glm()] with `family = binomial`, fitted with
#'   a `data` argument and a binary 0/1 (or two-level factor) outcome.
#' @param treatment Name of the binary treatment variable in the model data
#'   (numeric 0/1, logical, or a two-level factor whose second level is the
#'   treated level).
#' @param p1,p0 Prevalence of the confounder in the treated and untreated
#'   groups, conditional on the measured covariates.
#' @param threshold A practically important risk difference, on the scale
#'   of the estimate and with its sign (for example `0.02` for two extra
#'   events per 100 people). Default 0.
#' @param confounder_effect `"odds_ratio"` (default) or `"risk_difference"`:
#'   the scale of the confounder-outcome association.
#' @param target `"ATE"` (default; standardized to everyone in the model
#'   data) or `"ATT"` (standardized to the treated).
#' @param lambda Optional grid of strengths: odds ratios (>= 1) or risk
#'   differences (in \[0, 1)). By default from no confounding to 1.5 times
#'   the strength that removes the estimate (at most an odds ratio of 50).
#' @param alpha Two-sided significance level for the confidence limits in
#'   `meta` (default 0.05).
#' @param propensity Optional vector of propensity scores for the target
#'   `"ATE"`, used for the distribution of \eqn{U} in the population. By
#'   default a logistic regression of the treatment on the other terms of
#'   the model.
#'
#' For the risk-difference model, a scenario is used only if every implied
#' stratum risk lies in \eqn{[0, 1]}; the largest such strength is
#' `meta$max_feasible`, the default grid stops there, and infeasible
#' scenarios in a user grid are `NA`. Offsets in the model are kept.
#' `threshold` may lie on either side of zero; the crossing is the strength
#' at which the adjusted risk difference reaches it.
#'
#' @return A `confoundsens` object with `meta$framework = "risk_difference"`.
#'   `theta` is the adjusted risk difference and `se` its standard error;
#'   `meta$risks` holds the adjusted risks in both groups; `meta$lambda_zero`
#'   and `meta$lambda_threshold` are the strengths at which the estimate
#'   reaches 0 and `threshold`, and `meta$lambda_threshold_ci` the strength
#'   at which the confidence limit nearest zero reaches `threshold`. `NA`
#'   means not reached.
#' @references Lin, D. Y., Psaty, B. M., & Kronmal, R. A. (1998). Assessing
#'   the sensitivity of regression results to unmeasured confounders in
#'   observational studies. *Biometrics*, 54(3), 948--963.
#'   \doi{10.2307/2533848}
#'
#'   VanderWeele, T. J., & Arah, O. A. (2011). Bias formulas for sensitivity
#'   analysis of unmeasured confounding for general outcomes, treatments,
#'   and confounders. *Epidemiology*, 22(1), 42--52.
#'   \doi{10.1097/EDE.0b013e3181f74493}
#' @seealso [sens_surface()] for the plane of `p1` and the confounder
#'   strength.
#' @export
#'
#' @examples
#' set.seed(1)
#' n <- 2000
#' x <- rnorm(n)
#' a <- rbinom(n, 1, plogis(0.5 * x))
#' y <- rbinom(n, 1, plogis(-1.5 + 0.4 * a + 0.6 * x))
#' fit <- glm(y ~ a + x, family = binomial, data = data.frame(y, a, x))
#' rd <- sens_risk_difference(fit, "a", p1 = 0.4, p0 = 0.2, threshold = 0.02)
#' rd
#' sens_report(rd)
#' plot_robustness_curve(rd)
sens_risk_difference <- function(model, treatment, p1, p0, threshold = 0,
                                 confounder_effect = c("odds_ratio", "risk_difference"),
                                 target = c("ATE", "ATT"), lambda = NULL,
                                 alpha = 0.05, propensity = NULL) {
  type <- match.arg(confounder_effect)
  target <- match.arg(target)
  if (!is.numeric(alpha) || length(alpha) != 1L || is.na(alpha) ||
      alpha <= 0 || alpha >= 1) {
    stop("`alpha` must be a single number in (0, 1).", call. = FALSE)
  }
  if (!is.numeric(threshold) || length(threshold) != 1L || !is.finite(threshold)) {
    stop("`threshold` must be a single finite number.", call. = FALSE)
  }
  ob <- .rd_setup(model, treatment, p1, p0, target, propensity, type)
  e0 <- .rd_eval(ob, if (type == "odds_ratio") 1 else 0)
  rd0 <- e0$rd
  dir <- sign(rd0) * sign(p1 - p0)
  if (dir == 0) dir <- 1
  # strength lambda -> model parameter, oriented toward zero
  to_s <- function(l) if (type == "odds_ratio") l^dir else dir * l
  f <- function(l, what = "rd") .rd_eval(ob, to_s(l), se = what == "se")[[what]]
  lo <- if (type == "odds_ratio") 1 else 0
  max_feasible <- if (type == "odds_ratio") Inf else .rd_max_delta(ob, dir)
  hi <- if (type == "odds_ratio") 1e4 else min(1 - 1e-9, max_feasible)
  s0 <- if (rd0 >= 0) 1 else -1
  cross <- function(g) {
    g0 <- g(lo)
    if (!is.finite(g0)) return(NA_real_)
    if (g0 <= 0) return(lo)
    gh <- g(hi)
    if (!is.finite(gh) || gh > 0) return(NA_real_)
    stats::uniroot(g, c(lo, hi), tol = 1e-10)$root
  }
  l_zero <- cross(function(l) s0 * f(l))
  l_thr <- cross(function(l) s0 * (f(l) - threshold))
  l_thr_ci <- .rd_ci_cross(ob, to_s, lo, hi, s0, threshold, alpha)
  if (is.null(lambda)) {
    ref <- if (is.finite(l_zero)) l_zero else NA_real_
    top <- if (type == "odds_ratio") {
      if (is.na(ref)) 20 else min(50, max(3, 1 + 1.5 * (ref - 1)))
    } else {
      min(max_feasible, if (is.na(ref)) 0.5 else max(0.05, 1.5 * ref))
    }
    lambda <- seq(lo, top, length.out = 401L)
  } else {
    if (!is.numeric(lambda) || length(lambda) < 2L || anyNA(lambda) ||
        any(lambda < lo) || (type == "risk_difference" && any(lambda >= 1))) {
      stop(if (type == "odds_ratio") "`lambda` must be odds ratios >= 1."
           else "`lambda` must be risk differences in [0, 1).", call. = FALSE)
    }
    lambda <- sort(unique(lambda))
  }
  ev <- .rd_eval(ob, to_s(lambda))
  if (type == "risk_difference" && any(lambda > max_feasible)) {
    warning(sprintf(paste0(
      "Confounder-outcome risk differences above %s imply risks outside ",
      "[0, 1] in some confounder stratum; those scenarios are set to NA."),
      signif(max_feasible, 4)), call. = FALSE)
  }
  eff_lab <- if (type == "odds_ratio")
    "Odds ratio of U with the outcome (toward the null)" else
    "Risk difference of U with the outcome (toward the null)"
  new_confoundsens(
    lambda = lambda, theta = ev$rd, se = ev$se,
    meta = list(
      framework          = "risk_difference",
      source             = "glm",
      treatment          = treatment,
      estimate           = rd0,
      se0                = e0$se,
      risk1_obs          = e0$risk1,
      risk0_obs          = e0$risk0,
      null               = 0,
      rd_threshold       = threshold,
      p1                 = p1,
      p0                 = p0,
      confounder_effect  = type,
      direction          = dir,
      target             = target,
      default_alpha      = alpha,
      max_feasible       = max_feasible,
      lambda_zero        = l_zero,
      lambda_threshold   = l_thr,
      lambda_threshold_ci = l_thr_ci,
      risks              = data.frame(lambda = lambda, risk1 = ev$risk1,
                                      risk0 = ev$risk0),
      rd_setup           = ob,
      n                  = ob$n,
      lambda_label       = eff_lab,
      theta_label        = sprintf("Adjusted risk difference (%s)", target)
    )
  )
}

# ---- surface: p1 (a) by confounder strength (b), p0 fixed ---------------------
.surf_rd <- function(m, a, b, n, alpha) {
  ob <- m$rd_setup; type <- m$confounder_effect; dir <- m$direction
  if (is.null(alpha)) alpha <- m$default_alpha
  z <- stats::qnorm(1 - alpha / 2)
  thr <- m$rd_threshold; rd0 <- m$estimate
  s0 <- if (rd0 >= 0) 1 else -1
  lo <- if (type == "odds_ratio") 1 else 0
  top <- if (type == "odds_ratio") {
    if (is.na(m$lambda_zero)) 20 else min(50, max(3, 1 + 2 * (m$lambda_zero - 1)))
  } else min(0.99, if (is.na(m$lambda_zero)) 0.5 else 2 * m$lambda_zero)
  n <- min(n, 61L)   # each point is a g-computation with its gradient
  if (is.null(a)) a <- seq(0, 1, length.out = n) else a <- .check_grid(a, 0, 1, "a")
  if (is.null(b)) b <- seq(lo, top, length.out = n) else
    b <- .check_grid(b, lo, if (type == "odds_ratio") Inf else 0.999, "b")
  to_s <- function(l) if (type == "odds_ratio") l^dir else dir * l
  gr <- expand.grid(a = a, b = b)
  v <- .rd_eval(ob, to_s(gr$b), p1v = gr$a)
  gr$value <- v$rd
  past0 <- s0 * v$rd <= 0
  below <- s0 * (v$rd - thr) <= 0
  nsig <- s0 * v$rd - z * v$se <= 0
  gr$zone <- ifelse(is.na(v$rd), "Infeasible",
             ifelse(past0, "Estimate at or past zero",
             ifelse(below, "Below practical threshold",
             ifelse(nsig, "Not significant", NA_character_))))
  # Boundaries from contour lines of the grid values (each grid value is a
  # full g-computation, so no further root finding is done).
  cl_bnd <- function(v, kind, tag, style) {
    zm <- matrix(v, nrow = length(a), ncol = length(b))
    cl <- grDevices::contourLines(a, b, zm, levels = 0)
    if (!length(cl)) return(.empty_bnd())
    do.call(rbind, lapply(seq_along(cl), function(j)
      data.frame(a = cl[[j]]$x, b = cl[[j]]$y, kind = kind,
                 group = paste0(tag, j), style = style, stringsAsFactors = FALSE)))
  }
  bnd <- rbind(
    cl_bnd(s0 * v$rd, "Null (risk difference = 0)", "n", "dashed"),
    if (s0 * thr > 0) cl_bnd(s0 * (v$rd - thr),
                             sprintf("Practical threshold (%s)", thr), "t", "dotdash"),
    cl_bnd(s0 * v$rd - z * v$se, sprintf("Significance (alpha = %s)", alpha),
           "s", "solid"))
  list(grid = gr, boundary = bnd,
       meta = list(a_label = sprintf("Prevalence of U among the treated (p0 = %s)", ob$p0),
                   b_label = m$lambda_label,
                   value_label = sprintf("Adjusted risk difference (%s)", ob$target),
                   estimate = rd0, threshold = thr, p1 = ob$p1, p0 = ob$p0,
                   alpha = alpha))
}
