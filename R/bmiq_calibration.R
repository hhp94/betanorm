# Legacy BMIQ caps the outer EM at five iterations. The niter default and the
# legacy-drift warning must stay in sync, so both read this constant.
LEGACY_BMIQ_NITER <- 5L

# Windows used to locate each sample's unmethylated / methylated density modes
# when constructing initial EM thresholds: modes are estimated from values
# below MODE_WINDOW_LOW and above MODE_WINDOW_HIGH respectively.
MODE_WINDOW_LOW <- 0.4
MODE_WINDOW_HIGH <- 0.6

check_thresholds <- function(
  thresholds,
  nL,
  name,
  require.unit.interval = FALSE
) {
  if (!is.numeric(thresholds) || length(thresholds) != nL - 1L) {
    stop(
      name,
      " must be a numeric vector of length ",
      nL - 1L,
      ".",
      call. = FALSE
    )
  }
  if (any(!is.finite(thresholds))) {
    stop(name, " must contain only finite values.", call. = FALSE)
  }
  if (require.unit.interval &&
    any(thresholds <= 0 | thresholds >= 1)) {
    stop(name, " must lie strictly inside (0, 1).", call. = FALSE)
  }
  if (any(diff(thresholds) <= 0)) {
    stop(name, " must be strictly increasing.", call. = FALSE)
  }
  invisible(thresholds)
}

class_by_thresh <- function(beta, thresholds) {
  class <- rep.int(1L, length(beta))
  for (boundary in seq_along(thresholds)) {
    class[beta > thresholds[boundary]] <- boundary + 1L
  }
  class
}

require_all_classes <- function(class, nL, context, min.count = 1L) {
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

density_thresholds <- function(
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

  if (nL < 2L ||
    length(a) != nL ||
    length(b) != nL ||
    length(eta) != nL) {
    stop(context, " has inconsistent mixture dimensions.", call. = FALSE)
  }

  if (any(!is.finite(c(a, b, eta, means))) ||
    any(a <= 0) || any(b <= 0) || any(eta <= 0) ||
    any(means <= 0 | means >= 1) ||
    any(diff(means) <= 0)) {
    stop(
      context,
      " has invalid or unordered mixture parameters.",
      call. = FALSE
    )
  }

  # The desired boundary between adjacent components k and k + 1 (ordered by
  # increasing mean) is where the lower-mean component stops dominating the
  # weighted density and the higher-mean component takes over. Two Beta
  # densities can cross twice, so we require that specific orientation.
  find_crossing <- function(k) {
    lo <- means[k]
    hi <- means[k + 1L]

    da <- a[k] - a[k + 1L]
    db <- b[k] - b[k + 1L]

    constant <-
      log(eta[k]) - lbeta(a[k], b[k]) -
      log(eta[k + 1L]) + lbeta(a[k + 1L], b[k + 1L])

    log_score_diff <- function(x) {
      constant + da * log(x) + db * log1p(-x)
    }

    slope <- function(x) {
      da / x - db / (1 - x)
    }

    # The log-density ratio has at most one interior stationary point.
    cuts <- c(lo, hi)
    denominator <- da + db

    if (denominator != 0) {
      turning.point <- da / denominator
      if (is.finite(turning.point) &&
        turning.point > lo &&
        turning.point < hi) {
        cuts <- sort(c(lo, turning.point, hi))
      }
    }

    values <- vapply(cuts, log_score_diff, numeric(1L))
    if (any(!is.finite(values))) {
      stop(
        context,
        " produced a non-finite density ratio for boundary ",
        k,
        ".",
        call. = FALSE
      )
    }

    for (i in seq_len(length(cuts) - 1L)) {
      left <- cuts[i]
      right <- cuts[i + 1L]
      f.left <- values[i]
      f.right <- values[i + 1L]

      # Exact endpoint root with the required lower-to-higher orientation.
      if (f.left == 0 && slope(left) < 0) {
        return(left)
      }
      if (f.right == 0 && slope(right) < 0) {
        return(right)
      }

      # Component k dominates on the left and k + 1 on the right.
      if (f.left > 0 && f.right < 0) {
        return(
          stats::uniroot(
            log_score_diff,
            interval = c(left, right),
            tol = sqrt(.Machine$double.eps)
          )$root
        )
      }
    }

    stop(
      context,
      " has no lower-to-higher weighted-density crossing between means ",
      signif(lo, 8),
      " and ",
      signif(hi, 8),
      " for boundary ",
      k,
      ".",
      call. = FALSE
    )
  }

  thresholds <- vapply(
    seq_len(nL - 1L),
    find_crossing,
    numeric(1L)
  )

  if (any(diff(thresholds) <= 0)) {
    stop(
      context,
      " density boundaries are not ordinal: ",
      paste(signif(thresholds, 8), collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  thresholds
}

# Continuous class-wise quantile map for nL = 2, joined at the gold cut.
# Sample cut t_s and gold cut t_g are the U/M density-crossing thresholds.
# U side (x <= t_s): conditional quantile of F_sU on (0, t_s] -> (0, t_g].
# M side (x > t_s): conditional upper-tail map of F_sM on (t_s, 1] -> (t_g, 1].
# Shapes and thresholds come from fit_mixture()/density_thresholds(), which
# already guarantee length-2 canonical shapes and cuts strictly inside (0, 1).
normalize_nl2 <- function(
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
  # Work in log-probability throughout: an extreme but nonzero tail mass can
  # underflow pbeta() to exactly 0 on the natural scale and force a spurious
  # rejection, whereas log.p returns a finite log-mass. A conditional CDF
  # (numerator log-mass minus threshold log-mass) is clamped at 0 (probability
  # 1), then re-inflated by the gold threshold log-mass before qbeta().
  threshold_log_mass <- function(q, a, b, lower.tail, label) {
    value <- stats::pbeta(q, a, b, lower.tail = lower.tail, log.p = TRUE)
    # A real probability has log-mass <= 0; a tiny positive value is rounding
    # noise at prob = 1 and is harmless, but -Inf means the threshold sits on
    # the boundary with no usable conditioning mass.
    if (!is.finite(value)) {
      stop(
        context,
        " ",
        label,
        " at threshold is unusable (log = ",
        signif(value, 8),
        ").",
        call. = FALSE
      )
    }
    min(0, value)
  }

  log.FsU.ts <- threshold_log_mass(
    sample.threshold, sample.a[1L], sample.b[1L],
    lower.tail = TRUE, "sample U CDF"
  )
  log.FsM.ts <- threshold_log_mass(
    sample.threshold, sample.a[2L], sample.b[2L],
    lower.tail = FALSE, "sample M upper-tail CDF"
  )
  log.FgU.tg <- threshold_log_mass(
    gold.threshold, gold.a[1L], gold.b[1L],
    lower.tail = TRUE, "gold U CDF"
  )
  log.FgM.tg <- threshold_log_mass(
    gold.threshold, gold.a[2L], gold.b[2L],
    lower.tail = FALSE, "gold M upper-tail CDF"
  )

  map_tail <- function(x, k, lower.tail, log.sample.cut, log.gold.cut) {
    log.conditional <- pmin(
      0,
      stats::pbeta(
        x, sample.a[k], sample.b[k],
        lower.tail = lower.tail, log.p = TRUE
      ) - log.sample.cut
    )
    stats::qbeta(
      log.conditional + log.gold.cut,
      gold.a[k],
      gold.b[k],
      lower.tail = lower.tail,
      log.p = TRUE
    )
  }

  out <- as.numeric(beta)
  u_idx <- which(class == 1L)
  m_idx <- which(class == 2L)

  if (length(u_idx)) {
    out[u_idx] <- map_tail(
      out[u_idx], 1L,
      lower.tail = TRUE, log.FsU.ts, log.FgU.tg
    )
  }

  if (length(m_idx)) {
    out[m_idx] <- map_tail(
      out[m_idx], 2L,
      lower.tail = FALSE, log.FsM.ts, log.FgM.tg
    )
  }

  if (any(!is.finite(out))) {
    stop(context, " produced non-finite calibrated values.", call. = FALSE)
  }
  out
}

estimate_mode <- function(x, context) {
  if (!length(x)) {
    stop(context, " is empty; cannot estimate mode.", call. = FALSE)
  }
  if (all(x == x[1L])) {
    return(x[1L])
  }
  estimate <- density(x)
  # density() is not boundary-corrected, so its Gaussian tails can place the
  # mode just outside the unit interval; clamp so it can only shift the
  # initialization thresholds toward, never past, a valid beta value.
  mode <- estimate$x[which.max(estimate$y)]
  min(1, max(0, mode))
}

# Mixture component labels are arbitrary, so canonicalize each fit by sorting
# complete component tuples (a, b, eta, mu, responsibilities, fit status) by
# increasing mean. Downstream code can then treat component 1 as U and
# component nL as M, and density_thresholds() can assume adjacent components
# are ordered. Never sort a, b, eta, or thresholds independently.
canonicalize_em_components <- function(em, context) {
  a <- as.numeric(em$a[, 1L])
  b <- as.numeric(em$b[, 1L])
  eta <- as.numeric(em$eta)
  mu <- as.numeric(em$mu[, 1L])

  if (any(!is.finite(c(a, b, eta, mu))) ||
    any(a <= 0) || any(b <= 0) || any(eta <= 0) ||
    any(mu <= 0 | mu >= 1)) {
    stop(context, " returned invalid mixture parameters.", call. = FALSE)
  }

  ord <- order(mu)

  em$a <- em$a[ord, , drop = FALSE]
  em$b <- em$b[ord, , drop = FALSE]
  em$mu <- em$mu[ord, , drop = FALSE]
  em$eta <- em$eta[ord]
  em$w <- em$w[, ord, drop = FALSE]
  em$fit_status <- em$fit_status[ord]
  em$fit_reason <- em$fit_reason[ord]

  mu <- as.numeric(em$mu[, 1L])
  if (any(diff(mu) <= 0)) {
    stop(
      context,
      " has indistinguishable component means: ",
      paste(signif(mu, 8), collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  em
}

# Draw fit indices reproducibly without disturbing the caller's RNG stream.
# fit_mixture() runs once per sample, so seeding the global generator in place
# would silently reset .Random.seed for the whole session on every call.
draw_fit_indices <- function(n, size, seed) {
  has.seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  if (has.seed) {
    saved.seed <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
    on.exit(
      assign(".Random.seed", saved.seed, envir = globalenv()),
      add = TRUE
    )
  } else {
    on.exit(
      if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
        rm(".Random.seed", envir = globalenv())
      },
      add = TRUE
    )
  }
  set.seed(seed)
  sample.int(n, size, replace = FALSE)
}

fit_mixture <- function(
  beta,
  thresholds,
  nL,
  nfit,
  niter,
  tol,
  beta.maxit,
  beta.score.tol,
  context,
  debug = FALSE,
  seed = 1L,
  fit.idx = NULL
) {
  # The fit subset depends only on (length(beta), nfit, seed), so callers
  # looping over same-length samples pass a precomputed fit.idx instead of
  # redrawing the identical permutation for every sample.
  if (is.null(fit.idx)) {
    fit.idx <- draw_fit_indices(
      length(beta),
      min(nfit, length(beta)),
      seed
    )
  }
  beta.fit <- as.numeric(beta[fit.idx])

  initial.class <- class_by_thresh(beta.fit, thresholds)

  initial.counts <- require_all_classes(
    class = initial.class,
    nL = nL,
    context = paste0(context, " initial mixture"),
    min.count = 2L
  )

  # One-hot responsibilities for the fit subset only; allocating a full
  # length(beta) x nL matrix and then keeping only sampled rows wastes memory
  # proportional to the whole probe set.
  w.init <- matrix(0, nrow = length(fit.idx), ncol = nL)
  w.init[cbind(seq_along(fit.idx), initial.class)] <- 1

  # Historical BMIQ clips endpoints to half the distance to the nearest
  # interior observation. Guard the degenerate case where every fit value is
  # 0 (or every value is 1) so min()/max() do not return +/-Inf, and floor the
  # clips at endpoint.eps so a subnormal input cannot round the bound back to
  # an exact 0 or 1 (which would feed log(0) into the mixture fit).
  y_fit <- beta.fit
  endpoint.eps <- sqrt(.Machine$double.eps)
  positive <- y_fit[y_fit > 0]
  below.one <- y_fit[y_fit < 1]
  lower.clip <- if (length(positive)) {
    max(endpoint.eps, min(positive) / 2)
  } else {
    endpoint.eps
  }
  upper.clip <- if (length(below.one)) {
    min(1 - endpoint.eps, 1 - (1 - max(below.one)) / 2)
  } else {
    1 - endpoint.eps
  }
  if (!(lower.clip < upper.clip)) {
    stop(
      context,
      " could not construct a valid open-interval clipping range.",
      call. = FALSE
    )
  }
  y_fit <- pmin(upper.clip, pmax(lower.clip, y_fit))

  em <- beta_mixture_em_cpp(
    y = y_fit,
    initial_responsibility = w.init,
    nL = nL,
    maxiter = niter,
    tol = tol,
    beta_maxit = beta.maxit,
    beta_score_tol = beta.score.tol,
    debug = debug
  )

  em <- canonicalize_em_components(
    em,
    paste0(context, " mixture")
  )

  # Flatten the legacy K x 1 matrix shapes once at this boundary so all
  # downstream R code sees plain numeric vectors.
  em$a <- as.numeric(em$a[, 1L])
  em$b <- as.numeric(em$b[, 1L])
  em$mu <- as.numeric(em$mu[, 1L])

  # Components are already canonicalized by increasing mean.
  component.means <- em$mu

  # Diagnostic only: a valid soft component may never win the hard posterior
  # assignment, so do not reject on it. The load-bearing class checks are the
  # initial-count check above and the complete threshold-class check below.
  subset.class <- max.col(em$w, ties.method = "first")
  subset.counts <- tabulate(subset.class, nbins = nL)

  # The responsibility matrix is not returned in diagnostics and is only used
  # for the counts above; drop it so it is not carried for the whole run.
  em$w <- NULL

  posterior.thresholds <- density_thresholds(
    a = em$a,
    b = em$b,
    eta = as.numeric(em$eta),
    component.means = component.means,
    context = paste0(context, " posterior mixture")
  )

  full.class <- class_by_thresh(
    beta = beta,
    thresholds = posterior.thresholds
  )

  full.counts <- require_all_classes(
    class = full.class,
    nL = nL,
    context = paste0(context, " complete mixture"),
    min.count = 1L
  )

  list(
    em = em,
    random_indices = fit.idx,
    initial_class_counts = initial.counts,
    component_means = component.means,
    subset_map_counts = subset.counts,
    thresholds = posterior.thresholds,
    full_class = full.class,
    complete_class_counts = full.counts
  )
}

em_diagnostics <- function(fit, extra = NULL) {
  em <- fit$em
  out <- list(
    random_indices = fit$random_indices,
    initial_class_counts = fit$initial_class_counts,
    eta = em$eta,
    component_means = fit$component_means,
    component_a = em$a,
    component_b = em$b,
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

map_beta_q <- function(x, a.sample, b.sample, a.gold, b.gold, lower.tail) {
  qbeta(
    pbeta(x, a.sample, b.sample, lower.tail = lower.tail),
    a.gold,
    b.gold,
    lower.tail = lower.tail
  )
}

# Shared validation of the EM settings used by both bmiq_gold_fit() and
# bmiq_calibration(). Owns the legacy-drift warning so it fires exactly once
# per user-facing call. Returns the coerced integer settings.
validate_em_settings <- function(nL, nfit, niter, tol, beta.maxit,
                                 beta.score.tol) {
  nL <- as.integer(checkmate::assert_int(nL, lower = 2L, upper = 3L))
  nfit <- as.integer(checkmate::assert_int(nfit, lower = 2L * nL))
  niter <- as.integer(checkmate::assert_int(niter, lower = 1L))
  beta.maxit <- as.integer(checkmate::assert_int(beta.maxit, lower = 1L))
  checkmate::assert_number(tol, finite = TRUE)
  checkmate::assert_true(tol > 0, .var.name = "tol")
  checkmate::assert_number(beta.score.tol, finite = TRUE)
  checkmate::assert_true(beta.score.tol > 0, .var.name = "beta.score.tol")

  if (nL == 3L && niter > LEGACY_BMIQ_NITER) {
    warning(
      "nL = 3 with niter > ", LEGACY_BMIQ_NITER,
      " is not exactly compatible with legacy ",
      "five-iteration BMIQ results. This is expected if you intentionally ",
      "want the three-component fit to run further toward convergence.",
      call. = FALSE
    )
  }

  list(nL = nL, nfit = nfit, niter = niter, beta.maxit = beta.maxit)
}

# Fit the gold-standard mixture and reduce it to the parameters the per-sample
# calibration actually consumes. Settings are assumed validated by the caller.
fit_gold_standard <- function(
  goldstandard.beta,
  nL,
  nfit,
  th1.v,
  niter,
  tol,
  beta.maxit,
  beta.score.tol,
  debug,
  verbose
) {
  goldstandard.beta <- as.numeric(goldstandard.beta)
  checkmate::assert_numeric(
    goldstandard.beta,
    any.missing = FALSE,
    min.len = 2L * nL,
    .var.name = "goldstandard.beta"
  )
  scan_finite_unit_interval_cpp(
    goldstandard.beta,
    name = "goldstandard.beta",
    require_open = FALSE
  )

  if (is.null(th1.v)) {
    th1.v <- if (nL == 2L) 0.5 else c(0.2, 0.75)
  }
  check_thresholds(
    th1.v,
    nL = nL,
    name = "th1.v",
    require.unit.interval = TRUE
  )

  if (verbose) {
    message("Fitting EM beta mixture to gold-standard probes")
  }

  gold.fit <- fit_mixture(
    beta = goldstandard.beta,
    thresholds = th1.v,
    nL = nL,
    nfit = nfit,
    niter = niter,
    tol = tol,
    beta.maxit = beta.maxit,
    beta.score.tol = beta.score.tol,
    context = "Gold-standard",
    debug = debug
  )

  unmethylated.mode <- estimate_mode(
    goldstandard.beta[gold.fit$full_class == 1L],
    "Gold-standard unmethylated class"
  )
  methylated.mode <- estimate_mode(
    goldstandard.beta[gold.fit$full_class == nL],
    "Gold-standard methylated class"
  )

  if (verbose) {
    message("Gold-standard mixture fit complete")
  }

  structure(
    list(
      a = gold.fit$em$a,
      b = gold.fit$em$b,
      thresholds = as.numeric(gold.fit$thresholds),
      unmethylated.mode = unmethylated.mode,
      methylated.mode = methylated.mode,
      nL = nL,
      diagnostics = if (debug) {
        em_diagnostics(
          gold.fit,
          extra = list(
            unmethylated_mode = unmethylated.mode,
            methylated_mode = methylated.mode
          )
        )
      } else {
        NULL
      },
      settings = list(
        nfit = nfit,
        th1.v = th1.v,
        niter = niter,
        tol = tol,
        beta.maxit = beta.maxit,
        beta.score.tol = beta.score.tol
      )
    ),
    class = "bmiq_gold_fit"
  )
}

#' Fit the BMIQ Gold-Standard Beta Mixture Once
#'
#' Fits the gold-standard beta mixture used by [bmiq_calibration()] and
#' returns the fitted parameters. The gold standard is consumed purely
#' distributionally (fitted component shapes, density-crossing thresholds,
#' and the two class density modes) -- there is no per-probe alignment with
#' the matrix being calibrated. Fit it once with this function and pass the
#' result as `goldstandard.beta` to any number of [bmiq_calibration()] calls
#' (for example when calibrating a large matrix in chunks) instead of
#' refitting the same gold vector each time.
#'
#' @inheritParams bmiq_calibration
#' @param goldstandard.beta Numeric vector of gold-standard betas. Values
#'   must be finite, non-missing, and in \eqn{[0, 1]}; at least `2 * nL`
#'   values are required.
#'
#' @return An object of class `bmiq_gold_fit` with the fitted component
#'   shapes (`a`, `b`), the density-crossing `thresholds`, the class density
#'   modes (`unmethylated.mode`, `methylated.mode`), the `nL` used, fit
#'   `diagnostics` when `debug = TRUE` (else `NULL`), and the fit `settings`.
#'
#' @export
bmiq_gold_fit <- function(
  goldstandard.beta,
  nL = 3L,
  nfit = 20000L,
  th1.v = NULL,
  niter = LEGACY_BMIQ_NITER,
  tol = 0.001,
  beta.maxit = 50L,
  beta.score.tol = 1e-10,
  debug = FALSE,
  verbose = TRUE
) {
  checkmate::assert_flag(debug)
  checkmate::assert_flag(verbose)
  settings <- validate_em_settings(
    nL, nfit, niter, tol, beta.maxit, beta.score.tol
  )

  fit_gold_standard(
    goldstandard.beta,
    nL = settings$nL,
    nfit = settings$nfit,
    th1.v = th1.v,
    niter = settings$niter,
    tol = tol,
    beta.maxit = settings$beta.maxit,
    beta.score.tol = beta.score.tol,
    debug = debug,
    verbose = verbose
  )
}

#' Calibrate Methylation Beta Values Against a Gold Standard
#'
#' BMIQ-style calibration of DNA methylation beta values to a gold-standard
#' beta profile (beta-mixture quantile mapping). Provided for pipelines that
#' require this procedure; defaults target legacy BMIQ compatibility.
#'
#' @param datM Numeric matrix of beta values: samples in rows, CpGs in
#'   columns. Values must be finite, non-missing, and in \eqn{[0, 1]}.
#' @param goldstandard.beta Numeric vector of gold-standard betas, or a
#'   prefitted [bmiq_gold_fit()] object. The gold standard is consumed purely
#'   distributionally, so the vector does not need to align with (or match
#'   the length of) the columns of `datM`; a prefitted object lets one gold
#'   fit be reused across multiple calls.
#' @param nL Number of mixture components: `3` for unmethylated /
#'   intermediate / methylated (default, legacy three-state BMIQ), or `2`
#'   for unmethylated / methylated only. Choose `nL` by whether the
#'   intermediate component is scientifically appropriate for your data,
#'   not as a convergence setting: either model can be run further toward
#'   convergence by raising `niter` (see below). Using two components is not
#'   a mathematical substitute for converging a three-component model.
#' @param doH Whether to normalize the intermediate (H) component.
#'   Default is `TRUE` when `nL = 3` and `FALSE` when `nL = 2`.
#' @param nfit Maximum number of probes used when fitting each mixture.
#' @param th1.v Initial gold-standard class boundaries (length `nL - 1`).
#'   Defaults to `c(0.2, 0.75)` for `nL = 3` and `0.5` for `nL = 2`. Only
#'   used when fitting the gold standard, so it is ignored when
#'   `goldstandard.beta` is a prefitted [bmiq_gold_fit()] object.
#' @param niter Maximum outer EM iterations for the gold-standard fit and
#'   each sample fit. Default `5` is a legacy-compatibility setting matching
#'   common three-state BMIQ pipelines. Raising `niter` above `5` with
#'   `nL = 3` warns only that results will no longer match legacy
#'   five-iteration output; running a three-component fit further toward
#'   convergence is otherwise a valid choice.
#' @param tol Convergence tolerance for the mixture fit.
#' @param beta.maxit Maximum iterations for each Beta-component fit.
#' @param beta.score.tol Convergence tolerance for each Beta-component fit.
#' @param h.policy What to do if intermediate (H) normalization fails when
#'   `doH = TRUE`.
#'   * `"optional"` (default): keep the U-plus-upper-M calibration for that
#'     sample (lower-M observations, normally absorbed into the H map, are
#'     left unchanged).
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
#' @return An object of class `bmiq_calibration_result` with:
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
#' **Choosing `nL` and `niter` are independent decisions.** `nL` is a
#' modeling choice (is an intermediate component appropriate?); `niter` is a
#' stopping choice (how far to run EM). The `nL = 3`, `niter = 5` default is
#' the intentional legacy configuration for clock-style and other pipelines
#' that expect historical BMIQ behavior: the five-iteration cap is a
#' compatibility choice, not a claim of statistical convergence. Either the
#' two- or three-component model can be run further toward convergence by
#' raising `niter`; with `nL = 3` this emits a warning only because the output
#' then departs from legacy five-iteration results. Component fits are always
#' accepted when finite (the mixture is fit with a generalized-EM ascent guard
#' so accepting an unconverged component cannot decrease the log-likelihood);
#' per-component convergence status is reported in `diagnostics` when
#' `debug = TRUE`.
#'
#' With `nL = 3`, samples are fit as unmethylated (U), intermediate (H), and
#' methylated (M) components, then quantile-normalized to the gold standard
#' (with continuous H stitching when H runs). When H is skipped
#' (`h.policy = "optional"`), U and the upper-M tail are calibrated and
#' lower-M observations are left unchanged, matching legacy BMIQ.
#'
#' With `nL = 2`, only U and M are used (no H step). Normalization uses
#' truncated component quantile maps joined at the gold U/M threshold so both
#' sides send the sample cut to the gold cut (continuous at the boundary; a
#' slope kink is still possible).
#'
#' A failed gold-standard fit stops the whole call; sample failures follow
#' `on.sample.error`.
#'
#' @references
#' Teschendorff AE, Marabita F, Lechner M, Bartlett T, Tegner J,
#' Gomez-Cabrero D, Beck S (2013).
#' A beta-mixture quantile normalization method for correcting probe design
#' bias in Illumina Infinium 450 k DNA methylation data.
#' *Bioinformatics* 29(2), 189–196.
#' \doi{10.1093/bioinformatics/bts680}
#'
#' Horvath S (2013).
#' DNA methylation age of human tissues and cell types.
#' *Genome Biology* 14(10), R115.
#' \doi{10.1186/gb-2013-14-10-r115}
#'
#' @export
bmiq_calibration <- function(
  datM,
  goldstandard.beta,
  nL = 3L,
  doH = NULL,
  nfit = 20000L,
  th1.v = NULL,
  niter = LEGACY_BMIQ_NITER,
  tol = 0.001,
  beta.maxit = 50L,
  beta.score.tol = 1e-10,
  h.policy = c("optional", "require"),
  on.sample.error = c("stop", "continue"),
  failed.sample = c("NA", "original"),
  debug = FALSE,
  verbose = TRUE
) {
  call <- match.call()

  h.policy <- match.arg(h.policy)
  on.sample.error <- match.arg(on.sample.error)
  failed.sample <- match.arg(failed.sample)

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

  settings <- validate_em_settings(
    nL, nfit, niter, tol, beta.maxit, beta.score.tol
  )
  nL <- settings$nL
  nfit <- settings$nfit
  niter <- settings$niter
  beta.maxit <- settings$beta.maxit

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

  scan_finite_unit_interval_cpp(datM, name = "datM", require_open = FALSE)

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

  # Freshly allocated output so the C++ block scatter can mutate it in place
  # without ever touching the caller's datM (copy-on-write would otherwise make
  # the first write allocate a full-size copy anyway).
  calibrated <- matrix(
    NA_real_,
    nrow = number.of.samples,
    ncol = number.of.probes,
    dimnames = dimnames(datM)
  )

  success <- rep.int(FALSE, number.of.samples)
  h.applied.vec <- rep(NA, number.of.samples)
  failures <- list()

  sample.diagnostics <- if (debug) {
    vector("list", number.of.samples)
  } else {
    NULL
  }

  # The gold standard is consumed purely distributionally, so a prefitted
  # bmiq_gold_fit object can stand in for the raw vector and be reused across
  # calls (e.g. chunked calibration of one large matrix).
  if (inherits(goldstandard.beta, "bmiq_gold_fit")) {
    if (goldstandard.beta$nL != nL) {
      stop(
        "goldstandard.beta was fitted with nL = ",
        goldstandard.beta$nL,
        " but this call requested nL = ",
        nL,
        ".",
        call. = FALSE
      )
    }
    gold <- goldstandard.beta
  } else {
    gold <- fit_gold_standard(
      goldstandard.beta,
      nL = nL,
      nfit = nfit,
      th1.v = th1.v,
      niter = niter,
      tol = tol,
      beta.maxit = beta.maxit,
      beta.score.tol = beta.score.tol,
      debug = debug,
      verbose = verbose
    )
  }

  gold.a <- gold$a
  gold.b <- gold$b
  gold.thresholds <- gold$thresholds
  mod1U <- gold$unmethylated.mode
  mod1M <- gold$methylated.mode
  gold.diagnostics <- if (debug) gold$diagnostics else NULL

  # The fit subset is identical for every sample (same probe count, nfit, and
  # seed), so draw it once for the whole loop.
  sample.fit.idx <- draw_fit_indices(
    number.of.probes,
    min(nfit, number.of.probes),
    seed = 1L
  )

  process_sample <- function(ii, beta2.v) {
    beta2.v <- as.numeric(beta2.v)
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
        stage <- "sample mode estimation"

        low.mode.values <- beta2.v[beta2.v < MODE_WINDOW_LOW]
        high.mode.values <- beta2.v[beta2.v > MODE_WINDOW_HIGH]

        mod2U <- estimate_mode(
          low.mode.values,
          paste0("Sample ", ii, " values below ", MODE_WINDOW_LOW)
        )

        mod2M <- estimate_mode(
          high.mode.values,
          paste0("Sample ", ii, " values above ", MODE_WINDOW_HIGH)
        )

        if (debug) {
          diagnostic$low_mode_window_count <-
            length(low.mode.values)
          diagnostic$high_mode_window_count <-
            length(high.mode.values)
          diagnostic$unmethylated_mode <- mod2U
          diagnostic$methylated_mode <- mod2M
        }

        stage <- "initial threshold construction"

        unmethylated.shift <- mod2U - mod1U
        methylated.shift <- mod2M - mod1M

        # nL = 3: end-mode shifts on each boundary.
        # nL = 2: average of the two anchor shifts on the single U/M cut.
        if (nL == 3L) {
          th2.initial <- c(
            gold.thresholds[1L] + unmethylated.shift,
            gold.thresholds[2L] + methylated.shift
          )
        } else {
          th2.initial <- gold.thresholds[1L] +
            0.5 * (unmethylated.shift + methylated.shift)
        }

        # These are only EM initialization cuts, not fitted component-specific
        # roots, so sorting mode-shifted boundaries that cross is safe. (Never
        # sort fitted posterior thresholds.)
        th2.initial <- sort(th2.initial)

        check_thresholds(
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

        sample.fit <- fit_mixture(
          beta = beta2.v,
          thresholds = th2.initial,
          nL = nL,
          nfit = nfit,
          niter = niter,
          tol = tol,
          beta.maxit = beta.maxit,
          beta.score.tol = beta.score.tol,
          context = paste0("Sample ", ii),
          debug = debug,
          fit.idx = sample.fit.idx
        )

        em2.o <- sample.fit$em
        classAV2.v <- sample.fit$component_means
        class2.v <- sample.fit$full_class

        if (debug) {
          diagnostic <- modifyList(
            diagnostic,
            em_diagnostics(sample.fit)
          )
          diagnostic$posterior_thresholds <- sample.fit$thresholds
        }

        nbeta2.v <- beta2.v

        U <- 1L
        M <- nL

        # The tail split drives the nL = 3 maps and the H step; for nL = 2 it
        # is needed only for the debug tail-count diagnostic, so skip the
        # full-length passes on the plain nL = 2 path.
        if (nL == 3L || debug) {
          # Assign every U/M observation to exactly one tail; values exactly
          # at a component mean must not be left unnormalized.
          selU.idx <- which(class2.v == U)
          selUL.idx <- selU.idx[beta2.v[selU.idx] <= classAV2.v[U]]
          selUR.idx <- selU.idx[beta2.v[selU.idx] > classAV2.v[U]]
          selM.idx <- which(class2.v == M)
          selML.idx <- selM.idx[beta2.v[selM.idx] < classAV2.v[M]]
          selMR.idx <- selM.idx[beta2.v[selM.idx] >= classAV2.v[M]]
        }

        if (nL == 2L) {
          stage <- "nL=2 truncated U/M quantile normalization"
          nbeta2.v <- normalize_nl2(
            beta = beta2.v,
            class = class2.v,
            sample.a = em2.o$a,
            sample.b = em2.o$b,
            gold.a = gold.a,
            gold.b = gold.b,
            sample.threshold = sample.fit$thresholds[1L],
            gold.threshold = gold.thresholds[1L],
            context = paste0("Sample ", ii, " nL=2 map")
          )
          if (debug) {
            diagnostic$nl2_sample_threshold <- sample.fit$thresholds[1L]
            diagnostic$nl2_gold_threshold <- gold.thresholds[1L]
          }
        } else {
          stage <- "unmethylated quantile normalization"

          if (length(selUL.idx)) {
            nbeta2.v[selUL.idx] <- map_beta_q(
              beta2.v[selUL.idx],
              em2.o$a[U], em2.o$b[U],
              gold.a[U], gold.b[U],
              lower.tail = TRUE
            )
          }

          if (length(selUR.idx)) {
            nbeta2.v[selUR.idx] <- map_beta_q(
              beta2.v[selUR.idx],
              em2.o$a[U], em2.o$b[U],
              gold.a[U], gold.b[U],
              lower.tail = FALSE
            )
          }

          stage <- "methylated quantile normalization"

          if (length(selMR.idx)) {
            nbeta2.v[selMR.idx] <- map_beta_q(
              beta2.v[selMR.idx],
              em2.o$a[M], em2.o$b[M],
              gold.a[M], gold.b[M],
              lower.tail = FALSE
            )
          }
        }

        if (debug) {
          diagnostic$tail_counts <- c(
            U_left = length(selUL.idx),
            U_right = length(selUR.idx),
            M_left = length(selML.idx),
            M_right = length(selMR.idx)
          )
        }

        h.applied <- FALSE
        if (doH) {
          h.attempt <- tryCatch(
            {
              stage <- "intermediate/H normalization"

              # Component means are canonicalized in fit_mixture(), so no
              # separate ordering check is needed here.
              if (!length(selMR.idx)) {
                stop(
                  "H normalization needs methylated probes above the ",
                  "methylated-component mean.",
                  call. = FALSE
                )
              }

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
                "); U plus upper-M calibration (lower-M left unchanged)."
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
        if (debug) {
          diagnostic$success <- FALSE
          diagnostic$failure_stage <- stage
          diagnostic$failure_message <-
            conditionMessage(error)
        }

        original.message <- conditionMessage(error)
        stop(
          structure(
            list(
              message = paste0(
                "Sample ",
                ii,
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
              sample.index = ii,
              sample.name = sample.name,
              stage = stage,
              original.message = original.message,
              diagnostics = diagnostic,
              parent = error
            ),
            class = c("bmiq_sample_error", "error", "condition")
          )
        )
      }
    )
  }

  # Samples are rows in a column-major matrix, so a single sample is strided
  # across memory. Gather a contiguous run of samples into a probes x block
  # matrix (each sample contiguous), process the block, then scatter it back
  # into the freshly allocated output. Eight doubles is one 64-byte cache line.
  sample.block.size <- 8L

  for (block.start in seq.int(1L, number.of.samples, by = sample.block.size)) {
    block.count <- min(
      sample.block.size,
      number.of.samples - block.start + 1L
    )

    block <- gather_sample_block_cpp(
      datM,
      first_sample = block.start,
      sample_count = block.count
    )

    for (local.sample in seq_len(block.count)) {
      ii <- block.start + local.sample - 1L

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
        process_sample(ii, block[, local.sample]),
        bmiq_sample_error = function(error) error
      )

      if (inherits(attempt, "bmiq_sample_error")) {
        if (on.sample.error == "stop") {
          stop(attempt)
        }

        success[ii] <- FALSE

        if (failed.sample == "NA") {
          block[, local.sample] <- NA_real_
        }
        # failed.sample == "original": leave the gathered originals in place.

        failures[[length(failures) + 1L]] <- data.frame(
          sample_index = attempt$sample.index,
          sample_name = attempt$sample.name,
          stage = attempt$stage,
          message = attempt$original.message,
          stringsAsFactors = FALSE
        )

        if (debug) {
          sample.diagnostics[[ii]] <- attempt$diagnostics
        }

        next
      }

      block[, local.sample] <- attempt$beta
      success[ii] <- TRUE
      if (doH) {
        h.applied.vec[ii] <- attempt$h_applied
      }

      if (debug) {
        sample.diagnostics[[ii]] <- attempt$diagnostics
      }
    }

    scatter_sample_block_cpp(
      destination = calibrated,
      block = block,
      first_sample = block.start
    )

    rm(block)
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
    failures <- data.frame(
      sample_index = integer(),
      sample_name = character(),
      stage = character(),
      message = character(),
      stringsAsFactors = FALSE
    )
  }

  skipped.h <- which(success & !is.na(h.applied.vec) & !h.applied.vec)
  if (length(skipped.h)) {
    warning(
      length(skipped.h),
      " sample(s) calibrated without H (U plus upper-M). ",
      "See result$h.applied.",
      call. = FALSE
    )
  }

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
      h.policy = h.policy,
      on.sample.error = on.sample.error,
      failed.sample = failed.sample,
      debug = debug
    )
  )

  class(result) <- "bmiq_calibration_result"
  result
}

#' @rdname bmiq_calibration
#' @param x A `bmiq_calibration_result` object.
#' @param ... Unused.
#' @export
print.bmiq_calibration_result <- function(x, ...) {
  number.of.samples <- length(x$success)
  number.succeeded <- sum(x$success)
  number.failed <- number.of.samples - number.succeeded
  number.h.skipped <- sum(x$success & !is.na(x$h.applied) & !x$h.applied)

  cat("BMIQ calibration result\n")
  cat("  Samples:   ", number.of.samples, "\n", sep = "")
  cat("  Succeeded: ", number.succeeded, "\n", sep = "")
  cat("  Failed:    ", number.failed, "\n", sep = "")
  if (isTRUE(x$settings$doH)) {
    cat("  H skipped: ", number.h.skipped, " (U plus upper-M)\n", sep = "")
  }

  if (number.failed) {
    cat("\nFailures:\n")
    print(x$failures, row.names = FALSE)
  }

  invisible(x)
}

#' @rdname bmiq_calibration
#' @export
as.matrix.bmiq_calibration_result <- function(x, ...) {
  x$calibrated
}
