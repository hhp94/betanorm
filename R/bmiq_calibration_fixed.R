# Derived from Steve Horvath's BMIQcalibration (his adaptation of Andrew
# Teschendorff's BMIQ v1.2), whose EM is modified from RPMM (Houseman and
# Koestler, GPL >= 2). GPL (>= 3) here; see inst/COPYRIGHTS.

# bmiq_calibration_fixed() is bmiq_calibration() cut down to the legacy
# Horvath configuration (nL = 3 with H, niter = 5, nfit = 20000, tol = 0.001)
# with the gold standard fitted ahead of time and carried here as numbers.
# It is the file a downstream package vendors in place of
# R/bmiq_calibration.R, so it shares no code with that file: every helper
# here carries the `fixed_` prefix so the two can live in one package. The
# arithmetic is bmiq_calibration()'s, call for call, except that the
# gold side of the quantile map reads a precomputed table (bmiq_gold_qtable())
# in place of stats::qbeta(); tests/testthat/test-bmiq-fixed.R pins the
# output to bmiq_calibration()'s within 1e-14.

# The legacy BMIQ configuration. bmiq_calibration() takes these as arguments;
# here they are the only values.
FIXED_NFIT <- 20000L
FIXED_NITER <- 5L
FIXED_TOL <- 0.001
FIXED_BETA_MAXIT <- 50L
FIXED_BETA_SCORE_TOL <- 1e-10

# Windows used to locate each sample's unmethylated / methylated density modes
# when constructing initial EM thresholds: modes are estimated from values
# below FIXED_MODE_WINDOW_LOW and above FIXED_MODE_WINDOW_HIGH respectively.
FIXED_MODE_WINDOW_LOW <- 0.4
FIXED_MODE_WINDOW_HIGH <- 0.6

# The gold quantile table (bmiq_gold_qtable()): one half-table is the
# lower-tail quantile of a Beta(a, b) for p in [FIXED_QTABLE_P_MIN, 1/2] on
# FIXED_QTABLE_N nodes uniform in u = log p. The version is a schema id for
# tables a caller built earlier; bump it whenever what a table holds
# changes: the layout, the grid, or how the coefficients are built.
FIXED_QTABLE_VERSION <- 1L
FIXED_QTABLE_N <- 8192L
FIXED_QTABLE_P_MIN <- 1e-30

# Gold-standard fits, by name. Each is bmiq_gold_fit(<vector>, nL = 3) at its
# defaults (nfit = 20000, th1.v = c(0.2, 0.75), niter = 5, tol = 0.001,
# beta.maxit = 50, beta.score.tol = 1e-10), reduced to what the sample path
# reads: component shapes, density-crossing thresholds and the two class
# density modes. Written by data-raw/bmiq_gold_fixed.R as hexadecimal
# floating-point literals, which R parses exactly on every platform; the
# decimal beside each is for reading.
#
# horvath: horvath_goldstandard[["goldstandard2"]], the 21368 gold-standard
# means published with Horvath (2013), doi:10.1186/gb-2013-14-10-r115.
FIXED_GOLD <- list(
  horvath = list(
    a = c(
      0x1.280b80329a07fp+1, # 2.3128509757116826
      0x1.f1ae1ff450641p+0, # 1.9440631839246552
      0x1.cfa405d4f181cp+3 # 14.488772312073756
    ),
    b = c(
      0x1.321f91f333b22p+5, # 38.265415096297389
      0x1.1c24ec6423668p+1, # 2.2198768128353557
      0x1.6077f690e4582p+1 # 2.7536609847526856
    ),
    thresholds = c(
      0x1.6c635fe8efd87p-3, # 0.17792391709916108
      0x1.66e53f652f122p-1 # 0.7009677706733195
    ),
    unmethylated.mode = 0x1.e295d00893456p-6, # 0.02945466342707833
    methylated.mode = 0x1.bde9e34ffc3ffp-1 # 0.87092504836448359
  )
)

fixed_class_by_thresh <- function(beta, thresholds) {
  class <- rep.int(1L, length(beta))
  class[beta > thresholds[1L]] <- 2L
  class[beta > thresholds[2L]] <- 3L
  class
}

fixed_require_all_classes <- function(class, context, min.count) {
  counts <- tabulate(class, nbins = 3L)
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
  invisible(counts)
}

# The class boundaries of a canonical three-component fit (components in
# increasing mean order): where the lower-mean component stops dominating the
# weighted density and the higher-mean one takes over. Two Beta densities can
# cross twice, so that specific orientation is required.
fixed_density_thresholds <- function(a, b, eta, means, context) {
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

    # The log-density ratio has at most one interior stationary point.
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
            tol = sqrt(.Machine[["double.eps"]])
          )[["root"]]
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

  thresholds <- vapply(1:2, find_crossing, numeric(1L))

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

fixed_estimate_mode <- function(x, context) {
  if (!length(x)) {
    stop(context, " is empty; cannot estimate mode.", call. = FALSE)
  }
  if (all(x == x[1L])) {
    return(x[1L])
  }
  estimate <- stats::density(x)
  # density() is not boundary-corrected, so its Gaussian tails can place the
  # mode just outside the unit interval; clamp so it can only shift the
  # initialization thresholds toward, never past, a valid beta value.
  mode <- estimate[["x"]][which.max(estimate[["y"]])]
  min(1, max(0, mode))
}

# Draw fit indices reproducibly without disturbing the caller's RNG stream:
# the caller's .Random.seed (which also records its RNGkind()) is restored on
# exit, and the draw itself runs under R's default generator.
fixed_draw_fit_indices <- function(n, size, seed) {
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
  set.seed(
    seed,
    kind = "Mersenne-Twister",
    normal.kind = "Inversion",
    sample.kind = "Rejection"
  )
  sample.int(n, size, replace = FALSE)
}

# Fit the three-component mixture to beta[fit.idx], canonicalize it by
# increasing mean, and class every value of beta by the fitted density
# crossings. Returns the canonical shapes and means and the full classing.
fixed_fit_mixture <- function(beta, thresholds, fit.idx, context) {
  beta.fit <- as.numeric(beta[fit.idx])

  initial.class <- fixed_class_by_thresh(beta.fit, thresholds)
  fixed_require_all_classes(
    initial.class,
    context = paste0(context, " initial mixture"),
    min.count = 2L
  )

  w.init <- matrix(0, nrow = length(fit.idx), ncol = 3L)
  w.init[cbind(seq_along(fit.idx), initial.class)] <- 1

  # Historical BMIQ clips endpoints to half the distance to the nearest
  # interior observation. Guard the degenerate case where every fit value is
  # 0 (or every value is 1) so min()/max() do not return +/-Inf, and floor the
  # clips at endpoint.eps so a subnormal input cannot round the bound back to
  # an exact 0 or 1 (which would feed log(0) into the mixture fit).
  y_fit <- beta.fit
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
    nL = 3L,
    maxiter = FIXED_NITER,
    tol = FIXED_TOL,
    beta_maxit = FIXED_BETA_MAXIT,
    beta_score_tol = FIXED_BETA_SCORE_TOL
  )

  # Component labels are arbitrary, so sort the complete tuples (a, b, eta,
  # mu) by increasing mean: component 1 is then U and component 3 is M.
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
    stop(context, " mixture returned invalid mixture parameters.", call. = FALSE)
  }

  ord <- order(mu)
  a <- a[ord]
  b <- b[ord]
  eta <- eta[ord]
  mu <- mu[ord]

  if (any(diff(mu) <= 0)) {
    stop(
      context,
      " mixture has indistinguishable component means: ",
      paste(signif(mu, 8), collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  posterior.thresholds <- fixed_density_thresholds(
    a = a,
    b = b,
    eta = eta,
    means = mu,
    context = paste0(context, " posterior mixture")
  )

  class <- fixed_class_by_thresh(beta, posterior.thresholds)
  fixed_require_all_classes(
    class,
    context = paste0(context, " complete mixture"),
    min.count = 1L
  )

  list(a = a, b = b, mu = mu, class = class)
}

# The quantile map sends x through the sample's fitted Beta(a, b) to a tail
# probability p and back through the gold component's Beta with stats::qbeta.
# The gold shapes are constants, so the gold side is one fixed function per
# component and tail, and a precomputed table stands in for stats::qbeta.
#
# One half-table is the lower-tail quantile x = Q(p; a, b) for p in
# [FIXED_QTABLE_P_MIN, 1/2] on FIXED_QTABLE_N nodes uniform in u = log p,
# stored as the per-interval coefficients of a quintic Hermite interpolant
# of v = log x in t = (u - u_i) / h, matching the node values and the
# analytic first and second derivatives of v in u. A component table is
# two halves, Beta(a, b) and Beta(b, a): the reflection x -> 1 - x swaps the
# shapes and the tails, so Q_upper(p; a, b) = 1 - Q(p; b, a), and p above
# 1/2 is read as 1 - p from the swapped half. p = 0 and 1 are exact; p at
# or below FIXED_QTABLE_P_MIN falls back to stats::qbeta. Measured against
# stats::qbeta on the EPIC inputs and a dense sweep: within 1e-14
# (dev/qtable-exp.R, 2026-10-20).
#
# The grid is anchored at u1 = log(1/2): node i is u1 - (n - i) * h, and
# fixed_eval_half() recomputes a node position with that same expression,
# so build and evaluation agree on the nodes bit for bit. Nodes from seq()
# reconstructed as u0 + i * h drifted by up to 2.5e-11 at n = 4096
# (dev/qtable-exp.R), which put a 1e-12 floor under the interpolant.
fixed_qtable_grid <- function() {
  n <- FIXED_QTABLE_N
  u1 <- log(0.5)
  h <- (u1 - log(FIXED_QTABLE_P_MIN)) / (n - 1)
  list(n = n, u0 = u1 - (n - 1) * h, u1 = u1, h = h)
}

fixed_make_half <- function(a, b) {
  grid <- fixed_qtable_grid()
  n <- grid[["n"]]
  h <- grid[["h"]]
  u <- grid[["u1"]] - ((n - 1):0) * h
  x <- stats::qbeta(u, a, b, log.p = TRUE)
  v <- log(x)
  # dv/du = p / (x f(x)), and its derivative in u, both analytic.
  m <- exp(u - v - stats::dbeta(x, a, b, log = TRUE))
  k <- m * (1 - m * (a - (b - 1) * x / (1 - x)))

  # Per interval, the quintic in t whose value, first and second derivative
  # match both ends (derivatives with respect to t are h m and h^2 k).
  i0 <- seq_len(n - 1L)
  i1 <- i0 + 1L
  v0 <- v[i0]
  v1 <- v[i1]
  d0 <- h * m[i0]
  d1 <- h * m[i1]
  e0 <- h^2 * k[i0]
  e1 <- h^2 * k[i1]
  A <- v1 - v0 - d0 - e0 / 2
  B <- d1 - d0 - e0
  C <- e1 - e0
  list(
    a = a,
    b = b,
    n = n,
    u0 = grid[["u0"]],
    u1 = grid[["u1"]],
    h = h,
    c0 = v0,
    c1 = d0,
    c2 = e0 / 2,
    c3 = 10 * A - 4 * B + C / 2,
    c4 = -15 * A + 7 * B - C,
    c5 = 6 * A - 3 * B + C / 2
  )
}

fixed_make_qtable <- function(a, b) {
  list(
    a = a,
    b = b,
    ab = fixed_make_half(a, b),
    ba = fixed_make_half(b, a)
  )
}

#' Build the Gold Quantile Table for the Legacy Calibration
#'
#' Precomputes the gold-standard side of the quantile map that
#' [bmiq_calibration_fixed()] applies: for each of its two gold components,
#' the Beta quantile function on a fixed grid, interpolated at evaluation
#' time in place of [stats::qbeta()]. `bmiq_calibration_fixed()` builds the
#' table itself when `qtable = NULL`; building it once here and passing it
#' saves that step when many matrices are calibrated in one session.
#'
#' @param gold Name of the built-in gold standard, as in
#'   [bmiq_calibration_fixed()].
#'
#' @return A plain list: `version` (an integer schema id), `gold` (the
#'   name), and one table per gold component, `U` and `M`. Treat it as
#'   opaque: build it in the session that uses it and pass it only to
#'   `bmiq_calibration_fixed()` with the same `gold`, which checks it
#'   against the built-in gold shapes and grid on every call.
#'
#' @details
#' Each component table holds two half-tables, for Beta(a, b) and
#' Beta(b, a): the lower-tail quantile for p in \eqn{[10^{-30}, 1/2]} on
#' 8192 nodes uniform in \eqn{\log p}, stored as per-interval quintic
#' Hermite coefficients of \eqn{\log x}. The reflection \eqn{x \to 1 - x}
#' gives the upper tail and p above 1/2, p = 0 and 1 are exact, and p at
#' or below \eqn{10^{-30}} falls back to [stats::qbeta()]. Against
#' `stats::qbeta()` the interpolant is within 1e-14 on real inputs and a
#' dense sweep. Building takes about 70 ms; a table is about 1.5 MB.
#'
#' A table passed to [bmiq_calibration_fixed()] is checked for its layout,
#' version, gold shapes and grid, and spot-checked against
#' `stats::qbeta()` at sixteen interval midpoints per half-table.
#'
#' @examples
#' qtable <- bmiq_gold_qtable()
#' # Reuse it across calls:
#' # lapply(matrices, bmiq_calibration_fixed, qtable = qtable)
#'
#' @export
bmiq_gold_qtable <- function(gold = "horvath") {
  checkmate::assert_choice(gold, names(FIXED_GOLD))
  shapes <- FIXED_GOLD[[gold]]
  list(
    version = FIXED_QTABLE_VERSION,
    gold = gold,
    U = fixed_make_qtable(shapes[["a"]][1L], shapes[["b"]][1L]),
    M = fixed_make_qtable(shapes[["a"]][3L], shapes[["b"]][3L])
  )
}

# A table the caller built earlier can be stale or belong to another gold,
# so a passed table is checked against what bmiq_gold_qtable() would build
# now: schema version, gold name, every shape identical to FIXED_GOLD, the
# grid constants, and finite coefficient vectors of the right length. The
# structure says nothing about the values, so each half is also evaluated
# at FIXED_QTABLE_SPOT_CHECKS interval midpoints spread over the grid and
# compared with stats::qbeta: a shifted, foreign or differently built set
# of coefficients misses by far more than the 1e-12 allowed, while a table
# built here (within 1e-14) or on another platform's Rmath passes.
FIXED_QTABLE_SPOT_CHECKS <- 16L
FIXED_QTABLE_SPOT_TOL <- 1e-12
fixed_validate_qtable <- function(qtable, gold) {
  refuse <- function(...) {
    stop(
      "qtable ",
      ...,
      "; build it with bmiq_gold_qtable(\"",
      gold,
      "\").",
      call. = FALSE
    )
  }
  is_list_named <- function(x, names) {
    is.list(x) && identical(names(x), names)
  }

  if (!is_list_named(qtable, c("version", "gold", "U", "M"))) {
    refuse("is not a list with elements version, gold, U and M")
  }
  if (!identical(qtable[["version"]], FIXED_QTABLE_VERSION)) {
    refuse("was built by another version of bmiq_gold_qtable()")
  }
  if (!identical(qtable[["gold"]], gold)) {
    built.for <- qtable[["gold"]]
    if (is.character(built.for) && length(built.for) == 1L) {
      refuse("was built for gold \"", built.for, "\", not \"", gold, "\"")
    }
    refuse("was not built for gold \"", gold, "\"")
  }

  grid <- fixed_qtable_grid()
  coefficient.names <- paste0("c", 0:5)
  check_half <- function(half, a, b, where) {
    if (
      !is_list_named(
        half,
        c("a", "b", "n", "u0", "u1", "h", coefficient.names)
      )
    ) {
      refuse(where, " is not a half-table")
    }
    if (!identical(half[["a"]], a) || !identical(half[["b"]], b)) {
      refuse(where, " was built for other shapes")
    }
    if (
      !identical(half[["n"]], grid[["n"]]) ||
        !identical(half[["u0"]], grid[["u0"]]) ||
        !identical(half[["u1"]], grid[["u1"]]) ||
        !identical(half[["h"]], grid[["h"]])
    ) {
      refuse(where, " was built on another grid")
    }
    for (name in coefficient.names) {
      coefficients <- half[[name]]
      if (
        !is.double(coefficients) ||
          length(coefficients) != grid[["n"]] - 1L ||
          !all(is.finite(coefficients))
      ) {
        refuse(
          where,
          "$",
          name,
          " is not a finite double vector of length n - 1"
        )
      }
    }
    n <- grid[["n"]]
    intervals <- round(seq(1, n - 1, length.out = FIXED_QTABLE_SPOT_CHECKS))
    u.mid <- grid[["u1"]] - (n - intervals) * grid[["h"]] + grid[["h"]] / 2
    x.table <- fixed_eval_half(u.mid, half)
    x.rmath <- stats::qbeta(u.mid, a, b, log.p = TRUE)
    if (!all(abs(x.table - x.rmath) <= FIXED_QTABLE_SPOT_TOL)) {
      refuse(where, " does not reproduce stats::qbeta() for its shapes")
    }
  }
  check_component <- function(name, k) {
    component <- qtable[[name]]
    a <- FIXED_GOLD[[gold]][["a"]][k]
    b <- FIXED_GOLD[[gold]][["b"]][k]
    where <- paste0("$", name)
    if (!is_list_named(component, c("a", "b", "ab", "ba"))) {
      refuse(where, " is not a component table")
    }
    if (!identical(component[["a"]], a) || !identical(component[["b"]], b)) {
      refuse(where, " was built for other gold shapes")
    }
    check_half(component[["ab"]], a, b, paste0(where, "$ab"))
    check_half(component[["ba"]], b, a, paste0(where, "$ba"))
  }
  check_component("U", 1L)
  check_component("M", 3L)
  invisible(TRUE)
}

# x = Q(p; a, b) for u = log p in (u0, u1]: find the interval from the top
# of the grid, rebuild its left node exactly as fixed_make_half() placed it,
# evaluate the quintic by Horner, and exponentiate.
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

# stats::qbeta(p, a, b, lower.tail) from the component table. Valid inputs
# only: p is stats::pbeta() of scanned beta values, so it is in [0, 1] and
# never NA.
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

  # p at or below p.min on either side is Rmath's, as is any element the
  # table did not place in [0, 1] (validation only spot-checks a passed
  # table, so this is a fallback rather than a clip).
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

# Sample side as in bmiq_calibration(): the tail probability of x under the
# sample's fitted component. Gold side: the table.
fixed_map_beta_q <- function(x, a.sample, b.sample, gold.table, lower.tail) {
  fixed_qbeta_table(
    stats::pbeta(x, a.sample, b.sample, lower.tail = lower.tail),
    gold.table,
    lower.tail
  )
}

# The end-of-run warnings list the samples concerned. R prints at most
# getOption("warning.length") characters of a warning (1000 by default),
# so a long list is shortened here, to FIXED_WARNING_SAMPLES entries and a
# count, rather than cut off there.
FIXED_WARNING_SAMPLES <- 10L
fixed_list_samples <- function(items, sep) {
  if (length(items) > FIXED_WARNING_SAMPLES) {
    items <- c(
      items[seq_len(FIXED_WARNING_SAMPLES)],
      paste("and", length(items) - FIXED_WARNING_SAMPLES, "more")
    )
  }
  paste(items, collapse = sep)
}

#' Calibrate Beta Values to a Built-in Gold Standard (Legacy Configuration)
#'
#' [bmiq_calibration()] reduced to the legacy three-state BMIQ
#' configuration against a gold standard that was fitted ahead of time and
#' is built into the package. Every setting is fixed: `nL = 3` with H
#' normalization, `niter = 5`, `nfit = 20000`, `tol = 0.001`,
#' `beta.maxit = 50`, `beta.score.tol = 1e-10`, a skipped H step keeps the
#' U-plus-upper-M calibration (`h.policy = "optional"`), and a failed sample
#' becomes a row of `NA` (`on.sample.error = "continue"`). The output
#' agrees with the `calibrated` element of
#' `bmiq_calibration(datM, <gold vector>, on.sample.error = "continue")`
#' to within 1e-14: the gold mixture is not refitted on each call, and the
#' gold side of the quantile map reads a precomputed table
#' ([bmiq_gold_qtable()]) in place of [stats::qbeta()].
#'
#' @param datM Double-precision numeric matrix of beta values: samples in
#'   rows, CpGs in columns. Values must be finite, non-missing, and in
#'   \eqn{[0, 1]}; integer matrices are rejected rather than coerced.
#' @param gold Name of the built-in gold standard. `"horvath"` is
#'   [bmiq_gold_fit()] of `horvath_goldstandard[["goldstandard2"]]`
#'   (Horvath 2013) at its defaults. As in [bmiq_calibration()], the gold
#'   standard is consumed purely distributionally, so `datM` need not hold
#'   its probes.
#' @param verbose If `TRUE`, print progress messages.
#' @param qtable A table from [bmiq_gold_qtable()] for the same `gold`, or
#'   `NULL` to build it in the call (about 70 ms). Build it once and pass it
#'   when calibrating many matrices in a session. A passed table is checked
#'   against the built-in gold shapes and grid and refused if it does not
#'   match.
#'
#' @return The calibrated beta matrix, the size of `datM` with its dimnames.
#'   A sample that cannot be calibrated is a row of `NA`, and a warning gives
#'   its index and the error; the other samples are unaffected. Another
#'   warning names the samples whose H step was skipped; their H class and
#'   lower-M observations are left unchanged, as in [bmiq_calibration()].
#'
#' @details
#' The function exists for pipelines that run one fixed configuration and
#' copy the code: it has no settings to choose, no gold vector to fit, and
#' no result object. See [bmiq_calibration()] for what the
#' procedure does and how it departs from legacy BMIQ code, and for the
#' two-component model, error policies and diagnostics it offers.
#'
#' @references
#' Teschendorff AE, Marabita F, Lechner M, Bartlett T, Tegner J,
#' Gomez-Cabrero D, Beck S (2013).
#' A beta-mixture quantile normalization method for correcting probe design
#' bias in Illumina Infinium 450 k DNA methylation data.
#' *Bioinformatics* 29(2), 189-196.
#' \doi{10.1093/bioinformatics/bts680}
#'
#' Horvath S (2013).
#' DNA methylation age of human tissues and cell types.
#' *Genome Biology* 14(10), R115.
#' \doi{10.1186/gb-2013-14-10-r115}
#'
#' @export
bmiq_calibration_fixed <- function(
  datM,
  gold = "horvath",
  verbose = TRUE,
  qtable = NULL
) {
  checkmate::assert_choice(gold, names(FIXED_GOLD))
  checkmate::assert_flag(verbose)
  # Storage mode double is asserted (not coerced) so the C++ layer never
  # pays a per-call coercion copy; missing values are caught by the single
  # finite/range scan below, which owns that invariant.
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
    stop(
      "Too few probes (",
      number.of.probes,
      ") to fit nL = 3.",
      call. = FALSE
    )
  }

  if (is.null(qtable)) {
    qtable <- bmiq_gold_qtable(gold)
  } else {
    fixed_validate_qtable(qtable, gold)
  }

  sample.names <- rownames(datM)
  if (is.null(sample.names)) {
    sample.names <- rep.int("", number.of.samples)
  }

  gold <- FIXED_GOLD[[gold]]
  gold.thresholds <- gold[["thresholds"]]
  mod1U <- gold[["unmethylated.mode"]]
  mod1M <- gold[["methylated.mode"]]

  # Freshly allocated output so the C++ block scatter can mutate it in place
  # without ever touching the caller's datM.
  calibrated <- matrix(
    NA_real_,
    nrow = number.of.samples,
    ncol = number.of.probes,
    dimnames = dimnames(datM)
  )
  h.applied <- rep.int(FALSE, number.of.samples)
  failed <- rep.int(FALSE, number.of.samples)
  failure.messages <- character(number.of.samples)

  # The fit subset is identical for every sample (same probe count, nfit, and
  # seed), so draw it once for the whole loop.
  sample.fit.idx <- fixed_draw_fit_indices(
    number.of.probes,
    min(FIXED_NFIT, number.of.probes),
    seed = 1L
  )

  # beta2.v arrives as a bare double vector (a gathered block column).
  process_sample <- function(ii, beta2.v) {
    mod2U <- fixed_estimate_mode(
      beta2.v[beta2.v < FIXED_MODE_WINDOW_LOW],
      paste0("Sample ", ii, " values below ", FIXED_MODE_WINDOW_LOW)
    )
    mod2M <- fixed_estimate_mode(
      beta2.v[beta2.v > FIXED_MODE_WINDOW_HIGH],
      paste0("Sample ", ii, " values above ", FIXED_MODE_WINDOW_HIGH)
    )

    # End-mode shifts on each gold boundary. These are only EM initialization
    # cuts, not fitted component-specific roots, so sorting mode-shifted
    # boundaries that cross is safe.
    unmethylated.shift <- mod2U - mod1U
    methylated.shift <- mod2M - mod1M
    th2.initial <- sort(c(
      gold.thresholds[1L] + unmethylated.shift,
      gold.thresholds[2L] + methylated.shift
    ))
    if (any(diff(th2.initial) <= 0)) {
      stop(
        "Sample ",
        ii,
        " initial thresholds must be strictly increasing.",
        call. = FALSE
      )
    }

    sample.fit <- fixed_fit_mixture(
      beta = beta2.v,
      thresholds = th2.initial,
      fit.idx = sample.fit.idx,
      context = paste0("Sample ", ii)
    )
    a2 <- sample.fit[["a"]]
    b2 <- sample.fit[["b"]]
    classAV2.v <- sample.fit[["mu"]]
    class2.v <- sample.fit[["class"]]

    nbeta2.v <- beta2.v

    # Assign every U/M observation to exactly one tail; values exactly at a
    # component mean must not be left unnormalized.
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

    h.attempt <- tryCatch(
      {
        if (!length(selMR.idx)) {
          stop(
            "H normalization needs methylated probes above the ",
            "methylated-component mean.",
            call. = FALSE
          )
        }

        # Class 2 and lower-M (class 3) are disjoint, and
        # fixed_require_all_classes() already guarantees class 2 is not
        # empty, so the H set is never empty.
        selH.idx <- c(which(class2.v == 2L), selML.idx)

        minH <- min(beta2.v[selH.idx])
        maxH <- max(beta2.v[selH.idx])
        deltaH <- maxH - minH
        # Glue H straight to the calibrated U maximum and upper-M minimum, so
        # the map has no jump at either join. Legacy BMIQ also carries the
        # raw U-H and H-M gaps (deltaUH, deltaHM); dropping them is
        # deliberate (DECISIONS.md, "Legacy compatibility").
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

        TRUE
      },
      error = function(e) e
    )

    if (!isTRUE(h.attempt) && verbose) {
      message(
        "  H skipped for sample ",
        ii,
        " (",
        conditionMessage(h.attempt),
        "); U plus upper-M calibration (H and lower-M left unchanged)."
      )
    }

    if (
      any(!is.finite(nbeta2.v)) ||
        any(nbeta2.v < -1e-12 | nbeta2.v > 1 + 1e-12)
    ) {
      stop(
        "Sample ",
        ii,
        " normalization produced invalid beta values; range: [",
        paste(signif(range(nbeta2.v, finite = TRUE), 8), collapse = ", "),
        "].",
        call. = FALSE
      )
    }

    list(beta = pmin(1, pmax(0, nbeta2.v)), h_applied = isTRUE(h.attempt))
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

      # Samples are independent, so one that cannot be calibrated becomes an
      # NA row and the others proceed.
      result <- tryCatch(
        process_sample(ii, block[, local.sample]),
        error = function(error) error
      )
      if (inherits(result, "error")) {
        failed[ii] <- TRUE
        failure.messages[ii] <- conditionMessage(result)
        block[, local.sample] <- NA_real_
        if (verbose) {
          message(
            "  Sample ",
            ii,
            " failed (",
            failure.messages[ii],
            "); replaced with NA."
          )
        }
        next
      }
      block[, local.sample] <- result[["beta"]]
      h.applied[ii] <- result[["h_applied"]]
    }

    scatter_sample_block_cpp(
      destination = calibrated,
      block = block,
      first_sample = block.start
    )
  }

  failed.samples <- which(failed)
  if (length(failed.samples)) {
    warning(
      length(failed.samples),
      " sample(s) failed BMIQ calibration and were replaced with NA rows: ",
      fixed_list_samples(
        paste0(
          failed.samples,
          " (",
          failure.messages[failed.samples],
          ")"
        ),
        sep = "; "
      ),
      ".",
      call. = FALSE
    )
  }

  skipped.h <- which(!h.applied & !failed)
  if (length(skipped.h)) {
    warning(
      length(skipped.h),
      " sample(s) calibrated without H (U plus upper-M): ",
      fixed_list_samples(skipped.h, sep = ", "),
      ".",
      call. = FALSE
    )
  }

  calibrated
}
