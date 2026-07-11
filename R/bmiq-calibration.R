#' Internal Helpers for BMIQ Calibration
#'
#' Non-exported utilities that implement the corrected BMIQ pipeline used
#' by [BMIQcalibration()]. Policies include adaptive endpoint clipping for
#' EM, weighted Beta-density intersection thresholds (mixture-posterior
#' semantics), lower-class equality, hard failures for unusable classes,
#' continuous H stitching, dual-criterion EM convergence, and optional H
#' fallback. See [BMIQcalibration()] for the public API.
#'
#' @keywords internal
#' @name bmiq-calibration-helpers
#' @noRd
NULL

#' Clip Beta Values for Mixture Fitting
#'
#' Map raw beta values from the closed unit interval \eqn{[0, 1]} into the
#' open interval \eqn{(0, 1)} required by the Beta likelihood. Endpoint
#' handling matches historical `blc()` / `blc2()` exactly:
#' \preformatted{
#'   Ymn <- min(Y[Y > 0]);  Y <- pmax(Y, Ymn / 2)
#'   Ymx <- max(Y[Y < 1]);  Y <- pmin(Y, 1 - (1 - Ymx) / 2)
#' }
#' This preserves a data-driven scale for endpoint likelihood contributions.
#'
#' The raw domain is assumed to have already been checked by
#' `scan_finite_unit_interval_cpp()`. Callers must supply at least one
#' strictly positive and one strictly sub-unit value (as real methylation
#' profiles do); there is no machine-epsilon fallback.
#'
#' @param y Numeric vector of beta values in \eqn{[0, 1]}.
#'
#' @return Numeric vector of the same length as `y`, strictly inside
#'   \eqn{(0, 1)} when the blc extremes exist.
#'
#' @keywords internal
#' @noRd
clipBetaForFit <- function(y) {
  y <- as.numeric(y)
  # Exact blc/blc2 adaptive shrink (no fixed-eps fallback).
  y_min_pos <- min(y[y > 0])
  y_max_lt1 <- max(y[y < 1])
  y <- pmax(y, y_min_pos / 2)
  y <- pmin(y, 1 - (1 - y_max_lt1) / 2)
  y
}

#' Validate Ordered Mixture Thresholds
#'
#' Lightweight structural checks only (length, optional open unit interval,
#' strict increase). Domain scanning of beta matrices happens once at the
#' [BMIQcalibration()] entry point.
#'
#' @param thresholds Numeric vector of candidate boundaries.
#' @param nL Integer number of mixture components.
#' @param name Character label for error messages.
#' @param require.unit.interval Logical; if `TRUE`, thresholds must lie
#'   strictly inside \eqn{(0, 1)}.
#'
#' @return `NULL`, invisibly.
#'
#' @keywords internal
#' @noRd
validateThresholds <- function(
  thresholds,
  nL,
  name,
  require.unit.interval = FALSE
) {
  if (length(thresholds) != nL - 1L) {
    stop(name, " must have length nL - 1 (", nL - 1L, ").", call. = FALSE)
  }
  if (require.unit.interval &&
        any(!is.finite(thresholds) | thresholds <= 0 | thresholds >= 1)) {
    stop(name, " must lie strictly inside (0, 1).", call. = FALSE)
  }
  if (is.unsorted(thresholds, strictly = TRUE)) {
    stop(name, " must be strictly increasing.", call. = FALSE)
  }
  invisible(NULL)
}

#' Classify Values by Ordered Thresholds
#'
#' Assign each value to an ordered mixture class. Equality stays in the
#' lower class.
#'
#' For `nL = 3`:
#' * class 1: `x <= threshold[1]`
#' * class 2: `threshold[1] < x <= threshold[2]`
#' * class 3: `x > threshold[2]`
#'
#' For `nL = 2`:
#' * class 1: `x <= threshold[1]`
#' * class 2: `x > threshold[1]`
#'
#' @param beta Numeric vector of beta values.
#' @param thresholds Strictly increasing boundaries of length `nL - 1`.
#' @param nL Integer number of classes (`2` or `3`).
#'
#' @return Integer vector of class labels in `1:nL`.
#'
#' @keywords internal
#' @noRd
classifyByThresholds <- function(beta, thresholds, nL) {
  # Thresholds are validated at the API boundary / when constructed.
  class <- rep.int(1L, length(beta))
  for (boundary in seq_along(thresholds)) {
    class[beta > thresholds[boundary]] <- boundary + 1L
  }
  class
}

#' Build Hard One-Hot Responsibilities
#'
#' Construct an \eqn{n \times nL} responsibility matrix from threshold
#' classification. Each row is a one-hot indicator of the assigned class.
#'
#' @param beta Numeric vector of beta values.
#' @param thresholds Strictly increasing boundaries of length `nL - 1`.
#' @param nL Integer number of mixture components.
#'
#' @return Numeric matrix of hard responsibilities.
#'
#' @keywords internal
#' @noRd
initialResponsibilities <- function(beta, thresholds, nL) {
  class <- classifyByThresholds(beta, thresholds, nL)

  responsibility <- matrix(
    0,
    nrow = length(beta),
    ncol = nL
  )

  responsibility[
    cbind(seq_along(beta), class)
  ] <- 1

  responsibility
}

#' Require Nonempty Mixture Classes
#'
#' Fail if any class has fewer than `min.count` observations. Domain
#' scanners cannot guarantee nonempty mixture classes.
#'
#' @param class Integer class labels.
#' @param nL Integer number of classes.
#' @param context Character label for error messages.
#' @param min.count Minimum allowed count per class.
#'
#' @return Integer vector of class counts of length `nL`.
#'
#' @keywords internal
#' @noRd
requireAllClasses <- function(class, nL, context, min.count = 1L) {
  counts <- tabulate(class, nbins = nL)
  if (any(counts < min.count)) {
    stop(
      context,
      " has insufficient class counts [",
      paste(counts, collapse = ", "),
      "]; need >= ",
      min.count,
      " per class.",
      call. = FALSE
    )
  }
  counts
}

#' Thresholds from Weighted Beta-Density Intersections
#'
#' Build class boundaries by solving adjacent mixture-weighted density
#' crossings. For components \eqn{k} and \eqn{k+1}, find \eqn{t} with
#' \preformatted{
#'   log(eta_k) + log f_k(t) = log(eta_{k+1}) + log f_{k+1}(t)
#' }
#' where \eqn{f_j} is the Beta density of component \eqn{j}. The root is
#' sought between the adjacent component means (interval
#' \eqn{(\min(\mu_k, \mu_{k+1}), \max(\mu_k, \mu_{k+1}))}).
#'
#' This is the mixture-posterior / MAP decision boundary between two
#' components. It does **not** guarantee that each component mean lies
#' inside its MAP region (unlike component-mean midpoints).
#'
#' @param a,b Numeric component shape vectors of length `nL`.
#' @param eta Numeric mixture weights of length `nL` (positive, sum to 1).
#' @param component.means Numeric component means of length `nL`.
#' @param context Character label for error messages.
#'
#' @return Numeric vector of `nL - 1` strictly increasing thresholds.
#'
#' @keywords internal
#' @noRd
thresholdsFromDensityCrossings <- function(
  a,
  b,
  eta,
  component.means,
  context
) {
  a <- as.numeric(a)
  b <- as.numeric(b)
  eta <- as.numeric(eta)
  means <- as.numeric(component.means)
  nL <- length(means)

  if (nL < 2L) {
    stop(context, " needs at least two components.", call. = FALSE)
  }
  if (length(a) != nL || length(b) != nL || length(eta) != nL) {
    stop(
      context,
      " a, b, eta, and component.means must share length nL.",
      call. = FALSE
    )
  }
  if (any(!is.finite(a) | !is.finite(b) | !is.finite(eta) | !is.finite(means))) {
    stop(context, " mixture parameters must be finite.", call. = FALSE)
  }
  if (any(a <= 0) || any(b <= 0) || any(eta <= 0)) {
    stop(
      context,
      " shapes and mixture weights must be strictly positive.",
      call. = FALSE
    )
  }

  log_score_diff <- function(t, k) {
    # log(eta_k f_k) - log(eta_{k+1} f_{k+1})
    (log(eta[k]) + stats::dbeta(t, a[k], b[k], log = TRUE)) -
      (log(eta[k + 1L]) + stats::dbeta(t, a[k + 1L], b[k + 1L], log = TRUE))
  }

  find_crossing <- function(k) {
    lo <- min(means[k], means[k + 1L])
    hi <- max(means[k], means[k + 1L])
    if (!(hi > lo)) {
      stop(
        context,
        " adjacent component means are not separated for boundary ",
        k,
        " (means: ",
        paste(signif(means, 8), collapse = ", "),
        ").",
        call. = FALSE
      )
    }

    # Nudge off exact means so endpoint density quirks are less likely.
    span <- hi - lo
    left <- lo + 1e-8 * span
    right <- hi - 1e-8 * span
    if (!(right > left)) {
      left <- lo
      right <- hi
    }

    f_left <- log_score_diff(left, k)
    f_right <- log_score_diff(right, k)
    if (!is.finite(f_left) || !is.finite(f_right)) {
      stop(
        context,
        " non-finite weighted log-density at boundary ",
        k,
        " endpoints.",
        call. = FALSE
      )
    }
    if (f_left == 0) {
      return(left)
    }
    if (f_right == 0) {
      return(right)
    }

    if (f_left * f_right > 0) {
      # No sign change at mean endpoints: scan for a bracket inside (lo, hi).
      grid <- seq(left, right, length.out = 257L)
      vals <- vapply(grid, log_score_diff, numeric(1L), k = k)
      if (any(!is.finite(vals))) {
        stop(
          context,
          " non-finite weighted log-density on search grid for boundary ",
          k,
          ".",
          call. = FALSE
        )
      }
      sign_change <- which(vals[-length(vals)] * vals[-1L] <= 0)
      if (!length(sign_change)) {
        stop(
          context,
          " no weighted-density crossing between component means for ",
          "boundary ",
          k,
          " (means ",
          signif(lo, 8),
          ", ",
          signif(hi, 8),
          "; a = ",
          paste(signif(a[k:(k + 1L)], 6), collapse = "/"),
          "; b = ",
          paste(signif(b[k:(k + 1L)], 6), collapse = "/"),
          "; eta = ",
          paste(signif(eta[k:(k + 1L)], 6), collapse = "/"),
          ").",
          call. = FALSE
        )
      }
      i <- sign_change[[1L]]
      left <- grid[i]
      right <- grid[i + 1L]
      f_left <- vals[i]
      f_right <- vals[i + 1L]
      if (f_left == 0) {
        return(left)
      }
      if (f_right == 0) {
        return(right)
      }
    }

    root <- tryCatch(
      stats::uniroot(
        log_score_diff,
        interval = c(left, right),
        k = k,
        tol = .Machine$double.eps^0.5
      )$root,
      error = function(e) {
        stop(
          context,
          " failed to solve weighted-density crossing for boundary ",
          k,
          ": ",
          conditionMessage(e),
          call. = FALSE
        )
      }
    )
    root
  }

  thresholds <- vapply(seq_len(nL - 1L), find_crossing, numeric(1L))

  if (is.unsorted(thresholds, strictly = TRUE)) {
    stop(
      context,
      " density-crossing thresholds are not strictly increasing: ",
      paste(signif(thresholds, 8), collapse = ", "),
      " (means: ",
      paste(signif(means, 8), collapse = ", "),
      ").",
      call. = FALSE
    )
  }

  thresholds
}

#' Truncated Two-Class Beta Quantile Map (nL = 2)
#'
#' Continuous class-wise quantile normalization joined at the gold threshold.
#' Sample cut \(t_s\) and gold cut \(t_g\) are the U/M density-crossing
#' thresholds. For \(x \le t_s\) (U):
#' \preformatted{
#'   u = F_sU(x) / F_sU(t_s)
#'   g(x) = F_gU^{-1}(u * F_gU(t_g))
#' }
#' For \(x > t_s\) (M), upper tails for numerical stability:
#' \preformatted{
#'   r = (1 - F_sM(x)) / (1 - F_sM(t_s))
#'   g(x) = F_gM^{-1}(1 - r * (1 - F_gM(t_g)))
#' }
#' Both sides map \(t_s\) to \(t_g\), so the piecewise map is continuous at
#' the cut (slope kinks remain possible). Conditional quantiles within each
#' hard class are preserved.
#'
#' @param beta Numeric beta vector.
#' @param class Integer labels in `{1, 2}` (1 = U, 2 = M).
#' @param sample.a,sample.b Length-2 sample component shapes.
#' @param gold.a,gold.b Length-2 gold component shapes.
#' @param sample.threshold Scalar sample U/M cut \(t_s\).
#' @param gold.threshold Scalar gold U/M cut \(t_g\).
#' @param context Character label for error messages.
#'
#' @return Numeric vector of calibrated betas, same length as `beta`.
#'
#' @keywords internal
#' @noRd
normalizeNL2Truncated <- function(
  beta,
  class,
  sample.a,
  sample.b,
  gold.a,
  gold.b,
  sample.threshold,
  gold.threshold,
  context = "nL=2 truncated map"
) {
  beta <- as.numeric(beta)
  class <- as.integer(class)
  sample.a <- as.numeric(sample.a)
  sample.b <- as.numeric(sample.b)
  gold.a <- as.numeric(gold.a)
  gold.b <- as.numeric(gold.b)
  sample.threshold <- as.numeric(sample.threshold)[1L]
  gold.threshold <- as.numeric(gold.threshold)[1L]

  if (length(sample.a) != 2L || length(sample.b) != 2L ||
        length(gold.a) != 2L || length(gold.b) != 2L) {
    stop(context, " expects length-2 component shapes.", call. = FALSE)
  }
  if (!is.finite(sample.threshold) || !is.finite(gold.threshold) ||
        sample.threshold <= 0 || sample.threshold >= 1 ||
        gold.threshold <= 0 || gold.threshold >= 1) {
    stop(
      context,
      " thresholds must lie strictly in (0, 1); got sample = ",
      signif(sample.threshold, 8),
      ", gold = ",
      signif(gold.threshold, 8),
      ".",
      call. = FALSE
    )
  }

  # Sample / gold truncated-component denominators at the cuts.
  FsU_ts <- stats::pbeta(
    sample.threshold,
    sample.a[1L],
    sample.b[1L],
    lower.tail = TRUE
  )
  FsM_ts_upper <- stats::pbeta(
    sample.threshold,
    sample.a[2L],
    sample.b[2L],
    lower.tail = FALSE
  )
  FgU_tg <- stats::pbeta(
    gold.threshold,
    gold.a[1L],
    gold.b[1L],
    lower.tail = TRUE
  )
  FgM_tg_upper <- stats::pbeta(
    gold.threshold,
    gold.a[2L],
    gold.b[2L],
    lower.tail = FALSE
  )

  if (!(is.finite(FsU_ts) && FsU_ts > 0)) {
    stop(
      context,
      " sample U CDF at threshold is not positive (",
      signif(FsU_ts, 8),
      ").",
      call. = FALSE
    )
  }
  if (!(is.finite(FsM_ts_upper) && FsM_ts_upper > 0)) {
    stop(
      context,
      " sample M upper-tail CDF at threshold is not positive (",
      signif(FsM_ts_upper, 8),
      ").",
      call. = FALSE
    )
  }
  if (!(is.finite(FgU_tg) && FgU_tg > 0 && FgU_tg <= 1)) {
    stop(
      context,
      " gold U CDF at threshold is unusable (",
      signif(FgU_tg, 8),
      ").",
      call. = FALSE
    )
  }
  if (!(is.finite(FgM_tg_upper) && FgM_tg_upper > 0 && FgM_tg_upper <= 1)) {
    stop(
      context,
      " gold M upper-tail CDF at threshold is unusable (",
      signif(FgM_tg_upper, 8),
      ").",
      call. = FALSE
    )
  }

  out <- beta
  u_idx <- which(class == 1L)
  m_idx <- which(class == 2L)

  if (length(u_idx)) {
    # Conditional quantile of F_sU on (0, t_s], mapped into (0, t_g].
    u <- stats::pbeta(
      beta[u_idx],
      sample.a[1L],
      sample.b[1L],
      lower.tail = TRUE
    ) / FsU_ts
    u <- pmin(1, pmax(0, u))
    out[u_idx] <- stats::qbeta(
      u * FgU_tg,
      gold.a[1L],
      gold.b[1L],
      lower.tail = TRUE
    )
  }

  if (length(m_idx)) {
    # Conditional upper-tail quantile of F_sM on (t_s, 1], mapped into (t_g, 1].
    r <- stats::pbeta(
      beta[m_idx],
      sample.a[2L],
      sample.b[2L],
      lower.tail = FALSE
    ) / FsM_ts_upper
    r <- pmin(1, pmax(0, r))
    out[m_idx] <- stats::qbeta(
      r * FgM_tg_upper,
      gold.a[2L],
      gold.b[2L],
      lower.tail = FALSE
    )
  }

  if (any(!is.finite(out))) {
    stop(context, " produced non-finite calibrated values.", call. = FALSE)
  }
  out
}

#' Estimate a Univariate Mode
#'
#' Mode via `stats::density()`. Singletons and constant samples return the
#' constant value.
#'
#' @param x Finite numeric vector.
#' @param context Character label for error messages.
#'
#' @return Scalar mode estimate.
#'
#' @keywords internal
#' @noRd
estimateMode <- function(x, context) {
  # Domain already bulk-scanned at the API boundary; only emptiness remains
  # (e.g. no probes in a mode window).
  if (!length(x)) {
    stop(context, " is empty; cannot estimate mode.", call. = FALSE)
  }
  if (length(x) == 1L || all(x == x[1L])) {
    return(x[1L])
  }
  estimate <- density(x)
  estimate$x[which.max(estimate$y)]
}

#' Require Strictly Ordered Component Means
#'
#' Enforce complete ordering \eqn{\mu_U < \cdots < \mu_M}. Components are
#' not reordered; identities are fixed by initialization.
#'
#' Used for gold-standard fits and for H normalization. Optional H
#' fallback does **not** assert that successful sample fits always have
#' strictly ordered means. See [BMIQcalibration()] for the full policy.
#'
#' @param mu Numeric component means in component order.
#' @param context Character label for error messages.
#'
#' @return `mu` (as numeric), invisibly.
#'
#' @keywords internal
#' @noRd
requireOrderedComponentMeans <- function(mu, context) {
  mu <- as.numeric(mu)
  if (is.unsorted(mu, strictly = TRUE)) {
    stop(
      context,
      " component means are not strictly increasing: ",
      paste(signif(mu, 8), collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  invisible(mu)
}

#' Require Separated U and M Anchors
#'
#' Check only that the unmethylated and methylated component means are
#' separated (\eqn{\mu_U < \mu_M}). Enough for U/M quantile calibration;
#' full ordering is required only when applying H.
#'
#' Does not assert a valid ordered U/H/M interpretation of any interior
#' component. MAP labels may still depend on the interior through the
#' mixture fit.
#'
#' @param mu Numeric component means in component order.
#' @param context Character label for error messages.
#'
#' @return `mu` (as numeric), invisibly.
#'
#' @keywords internal
#' @noRd
requireOrderedAnchors <- function(mu, context) {
  mu <- as.numeric(mu)
  if (!length(mu) || any(!is.finite(mu))) {
    stop(context, " has missing or non-finite component means.", call. = FALSE)
  }
  if (mu[1L] >= mu[length(mu)]) {
    stop(
      context,
      " U and M component means are not separated (",
      signif(mu[1L], 8),
      " vs ",
      signif(mu[length(mu)], 8),
      "); U/M anchoring is impossible.",
      call. = FALSE
    )
  }
  invisible(mu)
}

#' Enforce EM Fit Policy
#'
#' Apply `fit.policy` to a Beta-mixture EM result.
#'
#' * `"usable"`: accept finite converged, max-iteration, or stalled
#'   component estimates. Unusable components already hard-fail in C++.
#' * `"converged"`: require the tightened outer-EM dual criterion
#'   (parameter state on \eqn{\log a}, \eqn{\log b}, \eqn{\eta} and
#'   relative log-likelihood, after at least two iterations) **and**
#'   every most-recent component optimizer to report convergence.
#'
#' @param em List returned by the C++ EM fitter.
#' @param fit.policy Either `"usable"` or `"converged"`.
#' @param context Character label for error messages.
#'
#' @return `NULL`, invisibly.
#'
#' @keywords internal
#' @noRd
enforceFitPolicy <- function(em, fit.policy, context) {
  if (fit.policy == "usable") {
    return(invisible(NULL))
  }
  if (!isTRUE(em$converged) || any(em$fit_status != "converged")) {
    stop(
      context,
      " did not fully converge under fit.policy = \"converged\" ",
      "(EM: ",
      isTRUE(em$converged),
      "; parameter_criterion: ",
      if (!is.null(em$parameter_criterion)) {
        signif(em$parameter_criterion, 8)
      } else {
        "NA"
      },
      "; loglik_criterion: ",
      if (!is.null(em$loglik_criterion)) {
        signif(em$loglik_criterion, 8)
      } else {
        "NA"
      },
      "; components: ",
      paste(em$fit_status, collapse = ", "),
      ").",
      call. = FALSE
    )
  }
  invisible(NULL)
}

#' Fit a Beta Mixture via C++
#'
#' Thin R wrapper around `beta_mixture_em_cpp()`. Clips endpoints into
#' \eqn{(0, 1)} with `clipBetaForFit()` (adaptive blc-style) before calling
#' C++.
#'
#' Pre-conditions (enforced at the API boundary; trusted by C++):
#' * `beta` already in \eqn{[0, 1]} (bulk-scanned);
#' * `initial.responsibility` hard one-hot, \eqn{n \times nL}, nonempty
#'   classes;
#' * hyperparameters already validated.
#'
#' @param beta Numeric beta vector.
#' @param initial.responsibility Hard responsibility matrix.
#' @param nL Number of components.
#' @param maxiter Maximum outer EM iterations.
#' @param tol Dual EM convergence tolerance.
#' @param beta.maxit Newton iterations per component.
#' @param beta.score.tol Newton score tolerance.
#' @param debug If `TRUE`, retain EM criterion traces.
#'
#' @return List from `beta_mixture_em_cpp()` (component shapes, means,
#'   responsibilities, dual criteria, fit statuses, and optional traces).
#'
#' @keywords internal
#' @noRd
fitBetaMixture <- function(
  beta,
  initial.responsibility,
  nL = 3L,
  maxiter = 5L,
  tol = 0.001,
  beta.maxit = 50L,
  beta.score.tol = 1e-10,
  debug = FALSE
) {
  beta_mixture_em_cpp(
    y = clipBetaForFit(beta),
    initial_responsibility = initial.responsibility,
    nL = nL,
    maxiter = maxiter,
    tol = tol,
    beta_maxit = beta.maxit,
    beta_score_tol = beta.score.tol,
    debug = debug
  )
}

#' Shared Gold / Sample Mixture Pipeline
#'
#' Run hard initialization, EM subsample fitting, fit-policy enforcement,
#' component-mean ordering checks, weighted-density threshold construction,
#' and full-vector classification.
#'
#' @param beta Numeric beta vector for one profile (gold or sample).
#' @param thresholds Initial class boundaries of length `nL - 1`.
#' @param nL Number of components (`2` or `3`).
#' @param nfit Maximum probes used in the EM subsample.
#' @param niter Maximum outer EM iterations.
#' @param tol Dual EM convergence tolerance.
#' @param beta.maxit Newton iterations per component.
#' @param beta.score.tol Newton score tolerance.
#' @param fit.policy `"usable"` or `"converged"`.
#' @param context Character label for messages and errors.
#' @param debug If `TRUE`, retain EM diagnostics.
#' @param seed Integer seed for the EM subsample.
#' @param mean.order `"strict"` requires full U \eqn{<} \eqn{\cdots} \eqn{<} M
#'   (gold-standard). `"anchors"` requires only \eqn{\mu_U < \mu_M}
#'   (sample pre-H / U/M-only).
#'
#' @return Named list with EM fit, indices, counts, component means,
#'   posterior thresholds, and full-vector class labels.
#'
#' @keywords internal
#' @noRd
fitMixturePipeline <- function(
  beta,
  thresholds,
  nL,
  nfit,
  niter,
  tol,
  beta.maxit,
  beta.score.tol,
  fit.policy,
  context,
  debug = FALSE,
  seed = 1L,
  # "strict": full U < ... < M (gold-standard). "anchors": U < M only
  # (sample pre-H / U/M-only; does not assert a valid interior H).
  mean.order = c("strict", "anchors")
) {
  mean.order <- match.arg(mean.order)
  w0 <- initialResponsibilities(
    beta = beta,
    thresholds = thresholds,
    nL = nL
  )

  set.seed(seed)
  rand.idx <- sample.int(
    length(beta),
    min(nfit, length(beta)),
    replace = FALSE
  )

  initial.counts <- requireAllClasses(
    class = max.col(w0[rand.idx, , drop = FALSE], ties.method = "first"),
    nL = nL,
    context = paste0(context, " initial mixture"),
    min.count = 2L
  )

  em <- fitBetaMixture(
    beta = beta[rand.idx],
    initial.responsibility = w0[rand.idx, , drop = FALSE],
    nL = nL,
    maxiter = niter,
    tol = tol,
    beta.maxit = beta.maxit,
    beta.score.tol = beta.score.tol,
    debug = debug
  )

  enforceFitPolicy(
    em,
    fit.policy = fit.policy,
    context = paste0(context, " mixture")
  )

  component.means <- as.numeric(em$mu[, 1L])
  if (mean.order == "strict") {
    requireOrderedComponentMeans(
      component.means,
      paste0(context, " mixture")
    )
  } else {
    requireOrderedAnchors(
      component.means,
      paste0(context, " mixture")
    )
  }

  subset.class <- max.col(em$w, ties.method = "first")
  subset.counts <- requireAllClasses(
    class = subset.class,
    nL = nL,
    context = paste0(context, " posterior mixture"),
    min.count = 1L
  )

  # Mixture-posterior boundaries: weighted Beta-density intersections
  # between adjacent components (not MAP-extrema or mean midpoints).
  posterior.thresholds <- thresholdsFromDensityCrossings(
    a = as.numeric(em$a[, 1L]),
    b = as.numeric(em$b[, 1L]),
    eta = as.numeric(em$eta),
    component.means = component.means,
    context = paste0(context, " posterior mixture")
  )

  full.class <- classifyByThresholds(
    beta = beta,
    thresholds = posterior.thresholds,
    nL = nL
  )

  full.counts <- requireAllClasses(
    class = full.class,
    nL = nL,
    context = paste0(context, " complete mixture"),
    min.count = 1L
  )

  list(
    em = em,
    random_indices = rand.idx,
    initial_class_counts = initial.counts,
    component_means = component.means,
    subset_map_counts = subset.counts,
    thresholds = posterior.thresholds,
    full_class = full.class,
    complete_class_counts = full.counts
  )
}

#' Collect EM Diagnostics
#'
#' Build a diagnostics list from a `fitMixturePipeline()` result. Used only
#' when `debug = TRUE`. Reports the dual EM criteria (parameter + relative
#' log-likelihood), not a mean-only criterion.
#'
#' @param fit Result of `fitMixturePipeline()`.
#' @param extra Optional named list of extra fields to append.
#'
#' @return Named list of diagnostic fields.
#'
#' @keywords internal
#' @noRd
collectEmDiagnostics <- function(fit, extra = NULL) {
  em <- fit$em
  out <- list(
    random_indices = fit$random_indices,
    initial_class_counts = fit$initial_class_counts,
    eta = em$eta,
    component_means = fit$component_means,
    component_a = as.numeric(em$a[, 1L]),
    component_b = as.numeric(em$b[, 1L]),
    component_fit_status = em$fit_status,
    component_fit_reason = em$fit_reason,
    em_iterations = em$iterations,
    em_converged = em$converged,
    parameter_criterion = em$parameter_criterion,
    loglik_criterion = em$loglik_criterion,
    log_likelihood = em$llike,
    subset_map_counts = fit$subset_map_counts,
    thresholds = fit$thresholds,
    complete_class_counts = fit$complete_class_counts
  )
  if (!is.null(em$parameter_criterion_trace)) {
    out$parameter_criterion_trace <- em$parameter_criterion_trace
  }
  if (!is.null(em$loglik_criterion_trace)) {
    out$loglik_criterion_trace <- em$loglik_criterion_trace
  }
  if (!is.null(extra)) {
    out <- c(out, extra)
  }
  out
}

#' Construct a Sample-Level Error Condition
#'
#' Build a `bmiq_sample_error` condition for a single-sample failure.
#' Sample-specific errors can be caught when
#' `on.sample.error = "continue"`. Gold-standard and global input errors
#' are never converted to this condition.
#'
#' @param sample.index Integer sample index (row of `datM`).
#' @param sample.name Character sample name (may be empty).
#' @param stage Character pipeline stage at failure.
#' @param error Underlying error condition.
#' @param diagnostics Optional diagnostics list when `debug = TRUE`.
#'
#' @return An object of class `c("bmiq_sample_error", "error", "condition")`.
#'
#' @keywords internal
#' @noRd
newBMIQSampleError <- function(
  sample.index,
  sample.name,
  stage,
  error,
  diagnostics = NULL
) {
  original.message <- conditionMessage(error)

  structure(
    list(
      message = paste0(
        "Sample ",
        sample.index,
        if (nzchar(sample.name)) {
          paste0(" (", sample.name, ")")
        } else {
          ""
        },
        " failed during ",
        stage,
        ": ",
        original.message
      ),
      call = NULL,
      sample.index = sample.index,
      sample.name = sample.name,
      stage = stage,
      original.message = original.message,
      diagnostics = diagnostics,
      parent = error
    ),
    class = c(
      "bmiq_sample_error",
      "error",
      "condition"
    )
  )
}

#' Empty Sample Failure Table
#'
#' Create a zero-row data frame with the failure-table schema used by
#' [BMIQcalibration()].
#'
#' @return Data frame with columns `sample_index`, `sample_name`, `stage`,
#'   and `message`.
#'
#' @keywords internal
#' @noRd
emptyBMIQFailureTable <- function() {
  data.frame(
    sample_index = integer(),
    sample_name = character(),
    stage = character(),
    message = character(),
    stringsAsFactors = FALSE
  )
}

#' Calibrate Methylation Beta Values Against a Gold Standard
#'
#' Adjust DNA methylation beta values so each sample matches a gold-standard
#' beta profile using beta-mixture quantile (BMIQ) normalization.
#'
#' @param datM Numeric matrix of beta values: samples in rows, CpGs in
#'   columns. Values must be finite, non-missing, and in \eqn{[0, 1]}.
#' @param goldstandard.beta Numeric vector of gold-standard betas, one per
#'   column of `datM`.
#' @param nL Number of mixture components: `3` for unmethylated /
#'   intermediate / methylated (default), or `2` for unmethylated /
#'   methylated only.
#' @param doH Whether to normalize the intermediate (H) component.
#'   Default is `TRUE` when `nL = 3` and `FALSE` when `nL = 2`.
#' @param nfit Maximum number of probes used when fitting each mixture.
#' @param th1.v Initial gold-standard class boundaries (length `nL - 1`).
#'   Defaults to `c(0.2, 0.75)` for `nL = 3` and `0.5` for `nL = 2`.
#' @param niter Maximum outer EM iterations.
#' @param tol Convergence tolerance for the mixture fit.
#' @param beta.maxit Maximum iterations for each Beta-component fit.
#' @param beta.score.tol Convergence tolerance for each Beta-component fit.
#' @param fit.policy How strict mixture convergence must be.
#'   * `"usable"` (default): accept finite component fits even if not fully
#'     converged.
#'   * `"converged"`: require full mixture and component convergence.
#' @param h.policy What to do if intermediate (H) normalization fails when
#'   `doH = TRUE`.
#'   * `"optional"` (default): keep U/M calibration for that sample.
#'   * `"require"`: treat the sample as failed.
#' @param on.sample.error What to do when a sample fails.
#'   * `"stop"` (default): abort.
#'   * `"continue"`: record the failure and process remaining samples.
#' @param failed.sample For `on.sample.error = "continue"`, write either
#'   an all-`NA` row (`"NA"`, default) or keep the original betas
#'   (`"original"`).
#' @param debug If `TRUE`, attach detailed fit diagnostics to the result.
#' @param verbose If `TRUE`, print progress messages.
#'
#' @return An object of class `BMIQcalibration_result` with:
#' \describe{
#'   \item{`calibrated`}{Calibrated beta matrix (same size as `datM`).}
#'   \item{`success`}{Logical vector of per-sample success.}
#'   \item{`failures`}{Data frame of failed samples (empty if none).}
#'   \item{`h.applied`}{Whether H normalization was applied per sample
#'     (`NA` when H was not requested).}
#'   \item{`diagnostics`}{Detailed diagnostics if `debug = TRUE`, else
#'     `NULL`.}
#' }
#'
#' Use `as.matrix()` to extract `calibrated`, and `print()` for a short
#' summary.
#'
#' @details
#' With `nL = 3`, samples are fit as unmethylated (U), intermediate (H), and
#' methylated (M) components, then quantile-normalized to the gold standard
#' (with continuous H stitching when H runs).
#'
#' With `nL = 2`, only U and M are used (no H step). Normalization uses
#' **truncated component quantile maps** joined at the gold U/M threshold:
#' both sides send the sample cut to the gold cut, so the map is continuous
#' at the boundary (a slope kink is still possible).
#'
#' A failed gold-standard fit stops the whole call. Sample failures follow
#' `on.sample.error`.
#'
#' @seealso [horvath_goldstandard], [GPL21145_sample]
#'
#' @export
BMIQcalibration <- function(
  datM,
  goldstandard.beta,
  nL = 3L,
  doH = NULL,
  nfit = 20000L,
  th1.v = NULL,
  niter = 5L,
  tol = 0.001,
  beta.maxit = 50L,
  beta.score.tol = 1e-10,
  fit.policy = c("usable", "converged"),
  h.policy = c("optional", "require"),
  on.sample.error = c("stop", "continue"),
  failed.sample = c("NA", "original"),
  debug = FALSE,
  verbose = TRUE
) {
  call <- match.call()

  fit.policy <- match.arg(fit.policy)
  h.policy <- match.arg(h.policy)
  on.sample.error <- match.arg(on.sample.error)
  failed.sample <- match.arg(failed.sample)

  # ---------------------------------------------------------------------------
  # API-boundary validation only. Helpers trust this contract and do not
  # re-assert matrix domain properties.
  # ---------------------------------------------------------------------------
  checkmate::assert_flag(debug)
  checkmate::assert_flag(verbose)
  checkmate::assert_matrix(
    datM,
    mode = "numeric",
    any.missing = FALSE,
    min.rows = 1L,
    min.cols = 1L
  )
  storage.mode(datM) <- "double"

  goldstandard.beta <- as.numeric(goldstandard.beta)
  checkmate::assert_numeric(
    goldstandard.beta,
    any.missing = FALSE,
    len = ncol(datM),
    .var.name = "goldstandard.beta"
  )

  nL <- as.integer(checkmate::assert_int(nL, lower = 2L, upper = 3L))

  if (is.null(th1.v)) {
    th1.v <- if (nL == 2L) 0.5 else c(0.2, 0.75)
  }
  if (is.null(doH)) {
    doH <- nL == 3L
  } else {
    checkmate::assert_flag(doH)
    if (nL == 2L && isTRUE(doH)) {
      stop(
        "doH = TRUE is not valid when nL = 2 (no H component).",
        call. = FALSE
      )
    }
  }

  nfit <- as.integer(checkmate::assert_int(nfit, lower = 2L * nL))
  niter <- as.integer(checkmate::assert_int(niter, lower = 1L))
  beta.maxit <- as.integer(checkmate::assert_int(beta.maxit, lower = 1L))
  checkmate::assert_number(tol, lower = 0, finite = TRUE)
  checkmate::assert_true(tol > 0, .var.name = "tol")
  checkmate::assert_number(beta.score.tol, lower = 0, finite = TRUE)
  checkmate::assert_true(beta.score.tol > 0, .var.name = "beta.score.tol")

  validateThresholds(
    th1.v,
    nL = nL,
    name = "th1.v",
    require.unit.interval = TRUE
  )

  # Finite + [0, 1] scan in C++ (hot path). Endpoints 0/1 are accepted
  # here; EM uses adaptive blc-style clipping (clipBetaForFit) into (0, 1).
  scan_finite_unit_interval_cpp(datM, name = "datM", require_open = FALSE)
  scan_finite_unit_interval_cpp(
    matrix(goldstandard.beta, ncol = 1L),
    name = "goldstandard.beta",
    require_open = FALSE
  )

  number.of.samples <- nrow(datM)
  number.of.probes <- ncol(datM)

  if (number.of.probes < 2L * nL) {
    stop(
      "Too few probes (",
      number.of.probes,
      ") to fit nL = ",
      nL,
      ".",
      call. = FALSE
    )
  }

  sample.names <- rownames(datM)
  if (is.null(sample.names)) {
    sample.names <- rep.int("", number.of.samples)
  }

  # Preserve the original matrix so failed rows can explicitly remain
  # uncalibrated when failed.sample = "original".
  original.datM <- datM
  calibrated <- datM

  success <- rep.int(FALSE, number.of.samples)
  # NA = H not requested (doH=FALSE) or sample failed; TRUE/FALSE after success.
  h.applied.vec <- rep(NA, number.of.samples)
  failures <- list()

  sample.diagnostics <- if (debug) {
    vector("list", number.of.samples)
  } else {
    NULL
  }

  # Gold-standard failures are global and are never skipped.
  beta1.v <- goldstandard.beta

  if (verbose) {
    message("Fitting EM beta mixture to gold-standard probes")
  }

  gold.fit <- fitMixturePipeline(
    beta = beta1.v,
    thresholds = th1.v,
    nL = nL,
    nfit = nfit,
    niter = niter,
    tol = tol,
    beta.maxit = beta.maxit,
    beta.score.tol = beta.score.tol,
    fit.policy = fit.policy,
    context = "Gold-standard",
    debug = debug
  )

  em1.o <- gold.fit$em
  classAV1.v <- gold.fit$component_means
  nth1.v <- gold.fit$thresholds

  mod1U <- estimateMode(
    beta1.v[gold.fit$full_class == 1L],
    "Gold-standard unmethylated class"
  )
  mod1M <- estimateMode(
    beta1.v[gold.fit$full_class == nL],
    "Gold-standard methylated class"
  )

  gold.diagnostics <- if (debug) {
    collectEmDiagnostics(
      gold.fit,
      extra = list(
        unmethylated_mode = mod1U,
        methylated_mode = mod1M
      )
    )
  } else {
    NULL
  }

  if (verbose) {
    message("Gold-standard mixture fit complete")
  }

  # ===========================================================================
  # One-sample worker
  #
  # All changes are made to a local beta vector. The output matrix is updated
  # only after the complete sample succeeds.
  # ===========================================================================
  processOneSample <- function(ii) {
    beta2.v <- as.numeric(original.datM[ii, ])
    sample.name <- sample.names[ii]
    stage <- "initialization"

    diagnostic <- if (debug) {
      list(
        sample_index = ii,
        sample_name = sample.name,
        input_range = range(beta2.v)
      )
    } else {
      NULL
    }

    tryCatch(
      {
        # ---------------------------------------------------------------------
        # Sample modes
        # ---------------------------------------------------------------------
        stage <- "sample mode estimation"

        low.mode.values <- beta2.v[beta2.v < 0.4]
        high.mode.values <- beta2.v[beta2.v > 0.6]

        mod2U <- estimateMode(
          low.mode.values,
          paste0("Sample ", ii, " values below 0.4")
        )

        mod2M <- estimateMode(
          high.mode.values,
          paste0("Sample ", ii, " values above 0.6")
        )

        if (debug) {
          diagnostic$low_mode_window_count <-
            length(low.mode.values)
          diagnostic$high_mode_window_count <-
            length(high.mode.values)
          diagnostic$unmethylated_mode <- mod2U
          diagnostic$methylated_mode <- mod2M
        }

        # ---------------------------------------------------------------------
        # Initial type-2 thresholds
        # ---------------------------------------------------------------------
        stage <- "initial threshold construction"

        unmethylated.shift <- mod2U - mod1U
        methylated.shift <- mod2M - mod1M

        # nL = 3: legacy end-mode shifts on each boundary.
        # nL = 2: average of the two anchor shifts on the single U/M boundary.
        if (nL == 3L) {
          th2.initial <- c(
            nth1.v[1L] + unmethylated.shift,
            nth1.v[2L] + methylated.shift
          )
        } else {
          th2.initial <- nth1.v[1L] +
            0.5 * (unmethylated.shift + methylated.shift)
        }

        validateThresholds(
          th2.initial,
          nL = nL,
          name = paste0("Sample ", ii, " initial thresholds"),
          require.unit.interval = FALSE
        )

        if (debug) {
          diagnostic$unmethylated_shift <- unmethylated.shift
          diagnostic$methylated_shift <- methylated.shift
          diagnostic$initial_thresholds <- th2.initial
        }

        stage <- "sample mixture fitting"

        sample.fit <- fitMixturePipeline(
          beta = beta2.v,
          thresholds = th2.initial,
          nL = nL,
          nfit = nfit,
          niter = niter,
          tol = tol,
          beta.maxit = beta.maxit,
          beta.score.tol = beta.score.tol,
          fit.policy = fit.policy,
          context = paste0("Sample ", ii),
          debug = debug,
          # U/M anchors only; full ordering required inside H when doH.
          mean.order = "anchors"
        )

        em2.o <- sample.fit$em
        classAV2.v <- sample.fit$component_means
        class2.v <- sample.fit$full_class

        if (debug) {
          diagnostic <- modifyList(
            diagnostic,
            collectEmDiagnostics(sample.fit)
          )
          diagnostic$posterior_thresholds <- sample.fit$thresholds
        }

        nbeta2.v <- beta2.v

        # U/M classes are nonempty (requireAllClasses on full MAP labels).
        U <- 1L
        M <- nL
        selU.idx <- which(class2.v == U)
        selUL.idx <- selU.idx[beta2.v[selU.idx] < classAV2.v[U]]
        selUR.idx <- selU.idx[beta2.v[selU.idx] > classAV2.v[U]]
        selM.idx <- which(class2.v == M)
        selML.idx <- selM.idx[beta2.v[selM.idx] < classAV2.v[M]]
        selMR.idx <- selM.idx[beta2.v[selM.idx] > classAV2.v[M]]

        if (nL == 2L) {
          # -----------------------------------------------------------------
          # nL = 2: truncated component quantile maps joined at gold cut.
          # Both sides send sample threshold t_s to gold threshold t_g
          # (continuous at the U/M boundary; slope kink still possible).
          # -----------------------------------------------------------------
          stage <- "nL=2 truncated U/M quantile normalization"
          nbeta2.v <- normalizeNL2Truncated(
            beta = beta2.v,
            class = class2.v,
            sample.a = as.numeric(em2.o$a[, 1L]),
            sample.b = as.numeric(em2.o$b[, 1L]),
            gold.a = as.numeric(em1.o$a[, 1L]),
            gold.b = as.numeric(em1.o$b[, 1L]),
            sample.threshold = sample.fit$thresholds[1L],
            gold.threshold = nth1.v[1L],
            context = paste0("Sample ", ii, " nL=2 map")
          )
          if (debug) {
            diagnostic$nl2_sample_threshold <- sample.fit$thresholds[1L]
            diagnostic$nl2_gold_threshold <- nth1.v[1L]
          }
        } else {
          # -----------------------------------------------------------------
          # nL = 3: classical left/right Beta quantile maps about component
          # means. M-left is left for H (or unchanged if H is skipped).
          # -----------------------------------------------------------------
          stage <- "unmethylated quantile normalization"

          if (length(selUL.idx)) {
            nbeta2.v[selUL.idx] <- qbeta(
              pbeta(
                beta2.v[selUL.idx],
                em2.o$a[U, 1L],
                em2.o$b[U, 1L],
                lower.tail = TRUE
              ),
              em1.o$a[U, 1L],
              em1.o$b[U, 1L],
              lower.tail = TRUE
            )
          }

          if (length(selUR.idx)) {
            nbeta2.v[selUR.idx] <- qbeta(
              pbeta(
                beta2.v[selUR.idx],
                em2.o$a[U, 1L],
                em2.o$b[U, 1L],
                lower.tail = FALSE
              ),
              em1.o$a[U, 1L],
              em1.o$b[U, 1L],
              lower.tail = FALSE
            )
          }

          stage <- "methylated quantile normalization"

          if (length(selMR.idx)) {
            nbeta2.v[selMR.idx] <- qbeta(
              pbeta(
                beta2.v[selMR.idx],
                em2.o$a[M, 1L],
                em2.o$b[M, 1L],
                lower.tail = FALSE
              ),
              em1.o$a[M, 1L],
              em1.o$b[M, 1L],
              lower.tail = FALSE
            )
          }
          # M-left remains raw for nL = 3 until H conformal map (or forever
          # if H is disabled / skipped).
        }

        if (debug) {
          diagnostic$tail_counts <- c(
            U_left = length(selUL.idx),
            U_right = length(selUR.idx),
            M_left = length(selML.idx),
            M_right = length(selMR.idx)
          )
        }

        # ---------------------------------------------------------------------
        # Intermediate/H conformal transformation (nL = 3 only)
        #
        # All stops in this body run before the single nbeta2.v[H] write, so a
        # failed attempt leaves U/M-only normalization intact. h.policy =
        # "optional" demotes to that U/M-only state; "require" fails the sample.
        # Optional fallback does not assert a valid ordered U/H/M fit.
        # ---------------------------------------------------------------------
        h.applied <- FALSE
        if (doH) {
          h.attempt <- tryCatch(
            {
              stage <- "intermediate/H normalization"

              requireOrderedComponentMeans(
                classAV2.v,
                paste0("Sample ", ii, " mixture (H)")
              )

              # Right-side M anchor is required for the conformal map.
              if (!length(selMR.idx)) {
                stop(
                  "H normalization needs methylated probes above the ",
                  "methylated-component mean.",
                  call. = FALSE
                )
              }

              # Single intermediate class (nL = 3): MAP-H plus legacy M-left.
              selH.idx <- unique(c(which(class2.v == 2L), selML.idx))
              if (!length(selH.idx)) {
                stop(
                  "Intermediate/H normalization set is empty.",
                  call. = FALSE
                )
              }

              minH <- min(beta2.v[selH.idx])
              maxH <- max(beta2.v[selH.idx])
              deltaH <- maxH - minH
              # Continuous piecewise map: glue H to the transformed U right
              # edge and M-right left edge (no raw-gap carry into output).
              # Slope kinks remain possible; jumps at the boundaries do not.
              nminH <- max(nbeta2.v[selU.idx])
              nmaxH <- min(nbeta2.v[selMR.idx])
              ndeltaH <- nmaxH - nminH

              if (!(deltaH > 0) || !(ndeltaH > 0)) {
                stop(
                  "H conformal map is degenerate ",
                  "(deltaH = ",
                  signif(deltaH, 8),
                  ", ndeltaH = ",
                  signif(ndeltaH, 8),
                  ").",
                  call. = FALSE
                )
              }

              hf <- ndeltaH / deltaH
              nbeta2.v[selH.idx] <-
                nminH + hf * (beta2.v[selH.idx] - minH)

              if (debug) {
                diagnostic$H_count <- length(selH.idx)
                diagnostic$H_input_range <- c(minH, maxH)
                diagnostic$H_output_anchors <- c(nminH, nmaxH)
                diagnostic$H_scale <- hf
              }

              TRUE
            },
            error = function(e) e
          )

          if (isTRUE(h.attempt)) {
            h.applied <- TRUE
          } else if (h.policy == "optional") {
            if (verbose) {
              message(
                "  H skipped for sample ",
                ii,
                " (",
                conditionMessage(h.attempt),
                "); U/M-only calibration."
              )
            }
            if (debug) {
              diagnostic$h_skip_reason <- conditionMessage(h.attempt)
            }
          } else {
            stop(h.attempt)
          }
        }

        stage <- "sample output validation"

        # Material out-of-range is an error; tiny roundoff is projected to [0, 1].
        if (any(!is.finite(nbeta2.v)) ||
          any(nbeta2.v < -1e-12 | nbeta2.v > 1 + 1e-12)) {
          stop(
            "Normalization produced invalid beta values; range: [",
            paste(signif(range(nbeta2.v, finite = TRUE), 8), collapse = ", "),
            "].",
            call. = FALSE
          )
        }
        nbeta2.v <- pmin(1, pmax(0, nbeta2.v))

        if (debug) {
          diagnostic$output_range <- range(nbeta2.v)
          diagnostic$h_applied <- if (doH) h.applied else NA
          diagnostic$success <- TRUE
        }

        list(
          beta = nbeta2.v,
          h_applied = h.applied,
          diagnostics = diagnostic
        )
      },
      error = function(error) {
        if (inherits(error, "bmiq_sample_error")) {
          stop(error)
        }

        if (debug) {
          diagnostic$success <- FALSE
          diagnostic$failure_stage <- stage
          diagnostic$failure_message <-
            conditionMessage(error)
        }

        stop(
          newBMIQSampleError(
            sample.index = ii,
            sample.name = sample.name,
            stage = stage,
            error = error,
            diagnostics = diagnostic
          )
        )
      }
    )
  }

  # ===========================================================================
  # Sample processing loop
  # ===========================================================================
  for (ii in seq_len(number.of.samples)) {
    if (verbose) {
      message(
        "Processing sample ",
        ii,
        " of ",
        number.of.samples,
        if (nzchar(sample.names[ii])) {
          paste0(" (", sample.names[ii], ")")
        } else {
          ""
        }
      )
    }

    attempt <- tryCatch(
      processOneSample(ii),
      bmiq_sample_error = function(error) error
    )

    if (inherits(attempt, "bmiq_sample_error")) {
      if (on.sample.error == "stop") {
        stop(attempt)
      }

      success[ii] <- FALSE

      if (failed.sample == "NA") {
        calibrated[ii, ] <- NA_real_
      } else {
        calibrated[ii, ] <- original.datM[ii, ]
      }

      failures[[length(failures) + 1L]] <- data.frame(
        sample_index = attempt$sample.index,
        sample_name = attempt$sample.name,
        stage = attempt$stage,
        message = attempt$original.message,
        stringsAsFactors = FALSE
      )

      if (debug) {
        sample.diagnostics[[ii]] <-
          attempt$diagnostics
      }

      next
    }

    calibrated[ii, ] <- attempt$beta
    success[ii] <- TRUE
    if (doH) {
      h.applied.vec[ii] <- attempt$h_applied
    }

    if (debug) {
      sample.diagnostics[[ii]] <-
        attempt$diagnostics
    }
  }

  if (length(failures)) {
    failures <- do.call(rbind, failures)
    rownames(failures) <- NULL

    warning(
      nrow(failures),
      " sample(s) failed BMIQ calibration and were ",
      if (failed.sample == "NA") {
        "replaced with NA rows"
      } else {
        "left uncalibrated"
      },
      ". See result$failures.",
      call. = FALSE
    )
  } else {
    failures <- emptyBMIQFailureTable()
  }

  skipped.h <- which(success & !is.na(h.applied.vec) & !h.applied.vec)
  if (length(skipped.h)) {
    warning(
      length(skipped.h),
      " sample(s) calibrated without H (U/M-only). ",
      "See result$h.applied.",
      call. = FALSE
    )
  }

  # Successful rows are validated locally before assignment; no full re-scan.

  diagnostics <- if (debug) {
    list(
      gold = gold.diagnostics,
      samples = sample.diagnostics
    )
  } else {
    NULL
  }

  result <- list(
    calibrated = calibrated,
    success = success,
    failures = failures,
    h.applied = h.applied.vec,
    diagnostics = diagnostics,
    call = call,
    settings = list(
      nL = nL,
      doH = doH,
      nfit = nfit,
      niter = niter,
      tol = tol,
      beta.maxit = beta.maxit,
      beta.score.tol = beta.score.tol,
      fit.policy = fit.policy,
      h.policy = h.policy,
      on.sample.error = on.sample.error,
      failed.sample = failed.sample,
      debug = debug
    )
  )

  class(result) <- "BMIQcalibration_result"
  result
}

#' @rdname BMIQcalibration
#' @param x A `BMIQcalibration_result` object.
#' @param ... Unused.
#' @export
print.BMIQcalibration_result <- function(x, ...) {
  number.of.samples <- length(x$success)
  number.succeeded <- sum(x$success)
  number.failed <- number.of.samples - number.succeeded
  number.h.skipped <- sum(x$success & !is.na(x$h.applied) & !x$h.applied)

  cat("BMIQ calibration result\n")
  cat("  Samples:   ", number.of.samples, "\n", sep = "")
  cat("  Succeeded: ", number.succeeded, "\n", sep = "")
  cat("  Failed:    ", number.failed, "\n", sep = "")
  if (isTRUE(x$settings$doH)) {
    cat("  H skipped: ", number.h.skipped, " (U/M-only)\n", sep = "")
  }

  if (number.failed) {
    cat("\nFailures:\n")
    print(x$failures, row.names = FALSE)
  }

  invisible(x)
}

#' @rdname BMIQcalibration
#' @export
as.matrix.BMIQcalibration_result <- function(x, ...) {
  x$calibrated
}
