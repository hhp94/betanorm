# Derived from Steve Horvath's BMIQcalibration (his adaptation of Andrew
# Teschendorff's BMIQ v1.2), whose EM is modified from RPMM (Houseman and
# Koestler, GPL >= 2). GPL (>= 3) here; see inst/COPYRIGHTS.

# bmiq_calibration() at the legacy configuration against a built-in gold fit,
# fitting every value instead of a random subset. Rationale: betanorm's
# dev/DECISIONS.md.

FIXED_NITER <- 5L
FIXED_TOL <- 0.001
FIXED_BETA_MAXIT <- 50L
FIXED_BETA_SCORE_TOL <- 1e-10
FIXED_MODE_WINDOW_LOW <- 0.4
FIXED_MODE_WINDOW_HIGH <- 0.6

# bmiq_gold_fit() of each gold vector with nfit at its length, written by
# data-raw/bmiq_gold_fixed.R.
# horvath: horvath_goldstandard[["goldstandard2"]] (Horvath 2013).
FIXED_GOLD <- list(
  horvath = list(
    a = c(
      0x1.26c84d67fae6ap+1, # 2.302987743159794
      0x1.efe9694136734p+0, # 1.9371553215227779
      0x1.cb03b492cf27cp+3 # 14.344202315079421
    ),
    b = c(
      0x1.31b22bd7139bp+5, # 38.211997680936634
      0x1.1a89fe4786af5p+1, # 2.2073362206696125
      0x1.5c6ed041d7eadp+1 # 2.7221317597867114
    ),
    thresholds = c(
      0x1.6c71421f60fadp-3, # 0.17795039804707127
      0x1.66f53872391b2p-1 # 0.70108963388752543
    ),
    unmethylated.mode = 0x1.e2a7a87660432p-6, # 0.029458918108051964
    methylated.mode = 0x1.bde9e34ffc3ffp-1 # 0.87092504836448359
  )
)

fixed_class_by_thresh <- function(beta, thresholds) {
  class <- rep.int(1L, length(beta))
  class[beta > thresholds[1L]] <- 2L
  class[beta > thresholds[2L]] <- 3L
  class
}

fixed_require_all_classes <- function(class, min.count) {
  if (any(tabulate(class, nbins = 3L) < min.count)) {
    stop("Too few values in a class.")
  }
}

# Where component k stops dominating the weighted density and component
# k + 1 takes over, for a fit sorted by mean.
fixed_density_thresholds <- function(a, b, eta, means) {
  find_crossing <- function(k) {
    lo <- means[k]
    hi <- means[k + 1L]

    da <- a[k] - a[k + 1L]
    db <- b[k] - b[k + 1L]

    constant <-
      log(eta[k]) -
      lbeta(a[k], b[k]) -
      log(eta[k + 1L]) +
      lbeta(a[k + 1L], b[k + 1L])

    log_score_diff <- function(x) {
      constant + da * log(x) + db * log1p(-x)
    }

    slope <- function(x) {
      da / x - db / (1 - x)
    }

    cuts <- c(lo, hi)
    denominator <- da + db

    if (denominator != 0) {
      turning.point <- da / denominator
      if (
        is.finite(turning.point) &&
          turning.point > lo &&
          turning.point < hi
      ) {
        cuts <- sort(c(lo, turning.point, hi))
      }
    }

    values <- vapply(cuts, log_score_diff, numeric(1L))
    if (any(!is.finite(values))) {
      stop("Non-finite density ratio.")
    }

    for (i in seq_len(length(cuts) - 1L)) {
      left <- cuts[i]
      right <- cuts[i + 1L]
      f.left <- values[i]
      f.right <- values[i + 1L]

      if (f.left == 0 && slope(left) < 0) {
        return(left)
      }
      if (f.right == 0 && slope(right) < 0) {
        return(right)
      }

      if (f.left > 0 && f.right < 0) {
        return(
          stats::uniroot(
            log_score_diff,
            interval = c(left, right),
            tol = sqrt(.Machine[["double.eps"]])
          )[["root"]]
        )
      }
    }

    stop("No density crossing.")
  }

  thresholds <- vapply(1:2, find_crossing, numeric(1L))

  if (any(diff(thresholds) <= 0)) {
    stop("Density crossings out of order.")
  }

  thresholds
}

fixed_estimate_mode <- function(x) {
  if (!length(x)) {
    stop("No values to estimate a mode from.")
  }
  if (all(x == x[1L])) {
    return(x[1L])
  }
  estimate <- stats::density(x)
  mode <- estimate[["x"]][which.max(estimate[["y"]])]
  min(1, max(0, mode))
}

# Fits the mixture to every value of beta, in order, and classes them.
fixed_fit_mixture <- function(beta, thresholds) {
  initial.class <- fixed_class_by_thresh(beta, thresholds)
  fixed_require_all_classes(initial.class, min.count = 2L)

  w.init <- matrix(0, nrow = length(beta), ncol = 3L)
  w.init[cbind(seq_along(beta), initial.class)] <- 1

  y_fit <- beta
  endpoint.eps <- sqrt(.Machine[["double.eps"]])
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
    stop("No valid clipping range.")
  }
  y_fit <- pmin(upper.clip, pmax(lower.clip, y_fit))

  em <- beta_mixture_em_cpp(
    y = y_fit,
    initial_responsibility = w.init,
    nL = 3L,
    maxiter = FIXED_NITER,
    tol = FIXED_TOL,
    beta_maxit = FIXED_BETA_MAXIT,
    beta_score_tol = FIXED_BETA_SCORE_TOL
  )

  a <- as.numeric(em[["a"]][, 1L])
  b <- as.numeric(em[["b"]][, 1L])
  eta <- as.numeric(em[["eta"]])
  mu <- as.numeric(em[["mu"]][, 1L])

  if (
    any(!is.finite(c(a, b, eta, mu))) ||
      any(a <= 0) ||
      any(b <= 0) ||
      any(eta <= 0) ||
      any(mu <= 0 | mu >= 1)
  ) {
    stop("Invalid mixture parameters.")
  }

  # Sort the components by mean: 1 is U, 3 is M.
  ord <- order(mu)
  a <- a[ord]
  b <- b[ord]
  eta <- eta[ord]
  mu <- mu[ord]

  if (any(diff(mu) <= 0)) {
    stop("Indistinguishable component means.")
  }

  posterior.thresholds <- fixed_density_thresholds(
    a = a,
    b = b,
    eta = eta,
    means = mu
  )

  class <- fixed_class_by_thresh(beta, posterior.thresholds)
  fixed_require_all_classes(class, min.count = 1L)

  list(a = a, b = b, mu = mu, class = class)
}

# FIXED_QTABLE (R/sysdata.rda, built by data-raw/bmiq_gold_qtable.R) holds,
# per gold set, components U and M, each two half-tables, Beta(a, b) and
# Beta(b, a): log Q(p) as quintic Hermite pieces in log p, p in [1e-30, 1/2].
fixed_gold_qtable <- function(gold) {
  qtable <- FIXED_QTABLE[[gold]]
  shapes <- FIXED_GOLD[[gold]]
  built_for <- function(component, k) {
    ctab <- qtable[[component]]
    a <- shapes[["a"]][k]
    b <- shapes[["b"]][k]
    identical(ctab[["a"]], a) &&
      identical(ctab[["b"]], b) &&
      identical(ctab[["ab"]][["a"]], a) &&
      identical(ctab[["ab"]][["b"]], b) &&
      identical(ctab[["ba"]][["a"]], b) &&
      identical(ctab[["ba"]][["b"]], a)
  }
  if (!built_for("U", 1L) || !built_for("M", 3L)) {
    stop("FIXED_QTABLE does not match FIXED_GOLD.")
  }
  qtable
}

# Node i is rebuilt as u1 - (n - i) * h, exactly as the builder placed it.
fixed_eval_half <- function(u, half) {
  n <- half[["n"]]
  h <- half[["h"]]
  u1 <- half[["u1"]]
  j <- pmin(floor((u1 - u) / h), n - 2)
  i <- n - 1 - j
  t <- (u - (u1 - (n - i) * h)) / h
  v <- half[["c0"]][i] +
    t *
      (half[["c1"]][i] +
        t *
          (half[["c2"]][i] +
            t *
              (half[["c3"]][i] +
                t *
                  (half[["c4"]][i] +
                    t * half[["c5"]][i]))))
  exp(v)
}

# stats::qbeta(p, a, b, lower.tail) from a component table.
fixed_qbeta_table <- function(p, gold.table, lower.tail) {
  ab <- gold.table[["ab"]]
  ba <- gold.table[["ba"]]
  u0 <- ab[["u0"]]
  u1 <- ab[["u1"]]

  u <- log(p)
  out <- numeric(length(u))
  zero <- u == -Inf
  one <- u >= 0
  small <- u > u0 & u <= u1 # p in (p.min, 1/2]
  ubig <- log(-expm1(u)) # log(1 - p)
  big <- u > u1 & u < 0 & ubig > u0 # p in (1/2, 1 - p.min)
  if (lower.tail) {
    out[small] <- fixed_eval_half(u[small], ab)
    out[big] <- 1 - fixed_eval_half(ubig[big], ba)
    out[zero] <- 0
    out[one] <- 1
  } else {
    out[small] <- 1 - fixed_eval_half(u[small], ba)
    out[big] <- fixed_eval_half(ubig[big], ab)
    out[zero] <- 1
    out[one] <- 0
  }

  # Off the grid, or not placed in [0, 1]: Rmath.
  fallback <- !(zero | one | small | big) | !is.finite(out) | out < 0 | out > 1
  if (any(fallback)) {
    out[fallback] <- stats::qbeta(
      p[fallback],
      gold.table[["a"]],
      gold.table[["b"]],
      lower.tail = lower.tail
    )
  }
  out
}

fixed_map_beta_q <- function(x, a.sample, b.sample, gold.table, lower.tail) {
  fixed_qbeta_table(
    stats::pbeta(x, a.sample, b.sample, lower.tail = lower.tail),
    gold.table,
    lower.tail
  )
}

#' Calibrate beta values to a built-in gold standard
#'
#' `bmiq_calibration()` at the legacy configuration (`nL = 3` with H,
#' `niter = 5`) against a built-in gold-standard fit, with every value of
#' the gold and of each sample fitted instead of a random subset of 20000.
#' A sample that cannot be calibrated becomes a row of `NA`.
#'
#' @param datM Double matrix of beta values in \eqn{[0, 1]}, samples in rows.
#' @param gold Name of the built-in gold standard.
#'
#' @return The calibrated matrix, the size of `datM`.
#'
#' @keywords internal
#' @noRd
bmiq_calibration_fixed <- function(datM, gold = "horvath") {
  checkmate::assert_choice(gold, names(FIXED_GOLD))
  checkmate::assert_matrix(
    datM,
    mode = "double",
    min.rows = 1L,
    min.cols = 1L
  )
  scan_finite_unit_interval_cpp(datM, name = "datM", require_open = FALSE)

  number.of.samples <- nrow(datM)
  number.of.probes <- ncol(datM)

  if (number.of.probes < 6L) {
    stop("Too few probes.")
  }

  qtable <- fixed_gold_qtable(gold)

  # Above 200 samples, progress in steps of 20% on one line.
  progress.at <- integer()
  if (number.of.samples > 200L) {
    progress.at <- (number.of.samples * 1:5) %/% 5L
    message("Calibrating ", number.of.samples, " samples:", appendLF = FALSE)
  }

  gold <- FIXED_GOLD[[gold]]
  gold.thresholds <- gold[["thresholds"]]
  mod1U <- gold[["unmethylated.mode"]]
  mod1M <- gold[["methylated.mode"]]

  calibrated <- matrix(
    NA_real_,
    nrow = number.of.samples,
    ncol = number.of.probes,
    dimnames = dimnames(datM)
  )
  h.applied <- rep.int(FALSE, number.of.samples)
  failed <- rep.int(FALSE, number.of.samples)

  process_sample <- function(beta2.v) {
    mod2U <- fixed_estimate_mode(beta2.v[beta2.v < FIXED_MODE_WINDOW_LOW])
    mod2M <- fixed_estimate_mode(beta2.v[beta2.v > FIXED_MODE_WINDOW_HIGH])

    unmethylated.shift <- mod2U - mod1U
    methylated.shift <- mod2M - mod1M
    th2.initial <- sort(c(
      gold.thresholds[1L] + unmethylated.shift,
      gold.thresholds[2L] + methylated.shift
    ))
    if (any(diff(th2.initial) <= 0)) {
      stop("Initial thresholds not increasing.")
    }

    sample.fit <- fixed_fit_mixture(
      beta = beta2.v,
      thresholds = th2.initial
    )
    a2 <- sample.fit[["a"]]
    b2 <- sample.fit[["b"]]
    classAV2.v <- sample.fit[["mu"]]
    class2.v <- sample.fit[["class"]]

    nbeta2.v <- beta2.v

    selU.idx <- which(class2.v == 1L)
    selUL.idx <- selU.idx[beta2.v[selU.idx] <= classAV2.v[1L]]
    selUR.idx <- selU.idx[beta2.v[selU.idx] > classAV2.v[1L]]
    selM.idx <- which(class2.v == 3L)
    selML.idx <- selM.idx[beta2.v[selM.idx] < classAV2.v[3L]]
    selMR.idx <- selM.idx[beta2.v[selM.idx] >= classAV2.v[3L]]

    if (length(selUL.idx)) {
      nbeta2.v[selUL.idx] <- fixed_map_beta_q(
        beta2.v[selUL.idx],
        a2[1L],
        b2[1L],
        qtable[["U"]],
        lower.tail = TRUE
      )
    }
    if (length(selUR.idx)) {
      nbeta2.v[selUR.idx] <- fixed_map_beta_q(
        beta2.v[selUR.idx],
        a2[1L],
        b2[1L],
        qtable[["U"]],
        lower.tail = FALSE
      )
    }
    if (length(selMR.idx)) {
      nbeta2.v[selMR.idx] <- fixed_map_beta_q(
        beta2.v[selMR.idx],
        a2[3L],
        b2[3L],
        qtable[["M"]],
        lower.tail = FALSE
      )
    }

    # H (class 2 and lower-M) is mapped linearly between the calibrated U
    # maximum and upper-M minimum; if it cannot be, H and lower-M stay as
    # they are.
    h.ok <- tryCatch(
      {
        if (!length(selMR.idx)) {
          stop("No upper-M values.")
        }

        selH.idx <- c(which(class2.v == 2L), selML.idx)

        minH <- min(beta2.v[selH.idx])
        maxH <- max(beta2.v[selH.idx])
        deltaH <- maxH - minH
        nminH <- max(nbeta2.v[selU.idx])
        nmaxH <- min(nbeta2.v[selMR.idx])
        ndeltaH <- nmaxH - nminH

        if (!(deltaH > 0) || !(ndeltaH > 0)) {
          stop("Degenerate H map.")
        }

        hf <- ndeltaH / deltaH
        nbeta2.v[selH.idx] <-
          nminH + hf * (beta2.v[selH.idx] - minH)

        TRUE
      },
      error = function(e) FALSE
    )

    if (
      any(!is.finite(nbeta2.v)) ||
        any(nbeta2.v < -1e-12 | nbeta2.v > 1 + 1e-12)
    ) {
      stop("Calibrated values outside [0, 1].")
    }

    list(beta = pmin(1, pmax(0, nbeta2.v)), h_applied = h.ok)
  }

  # Samples are rows of a column-major matrix: gather eight at a time into
  # contiguous columns, calibrate them, and scatter them back.
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
      result <- tryCatch(
        process_sample(block[, local.sample]),
        error = function(e) NULL
      )
      if (is.null(result)) {
        failed[ii] <- TRUE
        block[, local.sample] <- NA_real_
      } else {
        block[, local.sample] <- result[["beta"]]
        h.applied[ii] <- result[["h_applied"]]
      }
      step <- match(ii, progress.at)
      if (!is.na(step)) {
        message(
          if (step == 1L) " " else "...",
          20L * step,
          "%",
          appendLF = step == 5L
        )
      }
    }

    scatter_sample_block_cpp(
      destination = calibrated,
      block = block,
      first_sample = block.start
    )
  }

  if (any(failed)) {
    warning("Some samples could not be calibrated and are NA.")
  }
  if (any(!h.applied & !failed)) {
    warning("Some samples were calibrated without H.")
  }

  calibrated
}
