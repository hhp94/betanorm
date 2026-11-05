# bmiq_calibration_fixed() is bmiq_calibration() at its defaults with a
# built-in gold fit, except that every value is fitted where the general
# function fits a random subset of nfit. The stored gold numbers are
# bmiq_gold_fit() with nfit at the gold's length, pinned to the live fit bit
# for bit. Each sample is fitted on every value, in order, with no draw; the
# general function fits a sample the same way with nfit above its length
# and its draw replaced by the identity (local_fit_all_values()). So run,
# against the stored gold and over input the general function accepts, the
# output is pinned to bmiq_calibration() within 1e-14: the gold side of its
# quantile map is a precomputed table (FIXED_QTABLE, written to R/sysdata.rda
# by data-raw/bmiq_gold_qtable.R), not stats::qbeta(), so the two are not
# bit-identical. The table is pinned to the gold shapes and to stats::qbeta()
# directly.

horvath_vector <- function() {
  betanorm::horvath_goldstandard[["goldstandard2"]]
}

# The fit FIXED_GOLD stores (data-raw/bmiq_gold_fixed.R).
horvath_fit <- function() {
  gold <- horvath_vector()
  bmiq_gold_fit(gold, nL = 3L, nfit = length(gold), verbose = FALSE)
}

# With nfit at least the length, draw_fit_indices() returns every index; the
# identity puts them in order.
local_fit_all_values <- function(env = parent.frame()) {
  local_mocked_bindings(
    draw_fit_indices = function(n, size, seed) {
      stopifnot(size == n)
      seq_len(n)
    },
    .package = "betanorm",
    .env = env
  )
}

# bmiq_calibration() against the stored gold's fit, fitting each sample as
# the fixed path does. The gold is fitted before the draw is replaced.
general_all_values <- function(datM, ...) {
  gold <- horvath_fit()
  local_fit_all_values()
  bmiq_calibration(
    datM,
    gold,
    nfit = .Machine$integer.max,
    verbose = FALSE,
    ...
  )
}

# GPL21145_sample without the 1062 gold-standard probes the EPIC array does
# not measure (zero in every sample).
epic_measured <- function() {
  datM <- betanorm::GPL21145_sample
  datM[, colSums(datM == 0) < nrow(datM)]
}

# The numerics contract of the fixed function, over input the general function
# accepts: same NA pattern and dimnames, every other element within 1e-14
# (absolute).
expect_within_pin <- function(fixed, general) {
  expect_identical(dim(fixed), dim(general))
  expect_identical(dimnames(fixed), dimnames(general))
  expect_identical(is.na(fixed), is.na(general))
  expect_lte(max(abs(fixed - general), na.rm = TRUE), 1e-14)
}

# Samples in rows, named, over 1200 probes. No value is 0 or 1.
synthetic_samples <- function(n_samples) {
  n_probes <- 1200L
  datM <- t(vapply(
    seq_len(n_samples),
    function(i) simulate_beta_mixture(n_probes, nL = 3L, seed = 500L + i),
    numeric(n_probes)
  ))
  dimnames(datM) <- list(
    paste0("s", seq_len(n_samples)),
    paste0("cg", seq_len(n_probes))
  )
  datM
}

# The scanner refuses each of these. The last two are not caught by the
# +/-1e-12 tolerance on calibrated values, which is for output, not input.
bad_values <- list(
  na = NA_real_,
  nan = NaN,
  inf = Inf,
  minus.inf = -Inf,
  below = -1e-12,
  above = 1 + 1e-12
)

failed_samples_warning <- paste(
  "Some samples could not be calibrated, or hold a missing, infinite or",
  "out-of-range value, and are NA."
)

# Each bad value in turn at one probe of one sample: one warning, that row all
# NA (not NaN), dimnames kept, the other rows identical() to the clean run's.
expect_sample_fails_alone <- function(datM, row, probe) {
  clean <- bmiq_calibration_fixed(datM)
  expect_false(anyNA(clean))
  for (name in names(bad_values)) {
    bad <- datM
    bad[row, probe] <- bad_values[[name]]
    warnings <- capture_warnings(out <- bmiq_calibration_fixed(bad))
    expect_identical(warnings, failed_samples_warning, info = name)
    expect_identical(dimnames(out), dimnames(datM), info = name)
    expect_identical(unname(out[row, ]), rep(NA_real_, ncol(datM)), info = name)
    expect_identical(
      out[-row, , drop = FALSE],
      clean[-row, , drop = FALSE],
      info = name
    )
  }
}

# Dense p per (component, tail), after dev/qtable-exp.R without its random
# points: the grid nodes, three series of 20000 points spaced about 0.41 h
# apart from different starting offsets (about seven points per interval),
# their mirrors above 1/2, and the boundary cases. p is built first and
# the log taken from it, as in production.
sweep_p <- function(half) {
  eps <- .Machine$double.eps
  p.min <- 1e-30
  u0 <- half$u0
  u1 <- half$u1
  h <- half$h
  nodes <- u1 - ((half$n - 1):0) * h
  interior <- unlist(lapply(c(0.137, 0.5, 0.863), function(offset) {
    seq(u0 + offset * h, u1 - (1 - offset) * h, length.out = 20000L)
  }))
  low <- c(nodes, interior)
  sort(unique(c(
    exp(low),
    -expm1(low),
    0, 1,
    0.5, 0.5 * (1 - eps / 2), 0.5 * (1 + eps),
    p.min * (1 - eps), p.min, p.min * (1 + eps), 1e-300, .Machine$double.xmin,
    1 - eps / 2, 1 - eps, 1 - 1e-12, 1 - 1e-8
  )))
}

test_that("the stored Horvath gold fit is the live bmiq_gold_fit() of every value", {
  live <- horvath_fit()
  stored <- betanorm:::FIXED_GOLD[["horvath"]]

  expect_identical(names(betanorm:::FIXED_GOLD), "horvath")
  expect_identical(stored$a, live$a)
  expect_identical(stored$b, live$b)
  expect_identical(stored$thresholds, live$thresholds)
  expect_identical(stored$unmethylated.mode, live$unmethylated.mode)
  expect_identical(stored$methylated.mode, live$methylated.mode)
})

test_that("the stored table is built for every gold set on the documented grid", {
  expect_identical(names(FIXED_QTABLE), names(FIXED_GOLD))
  half.names <- c("a", "b", "n", "u0", "u1", "h", paste0("c", 0:5))
  for (gold in names(FIXED_GOLD)) {
    qtable <- FIXED_QTABLE[[gold]]
    # The shapes check every call makes.
    expect_identical(fixed_gold_qtable(gold), qtable)
    expect_identical(names(qtable), c("U", "M"))
    for (component in c("U", "M")) {
      expect_identical(names(qtable[[component]]), c("a", "b", "ab", "ba"))
      for (name in c("ab", "ba")) {
        half <- qtable[[component]][[name]]
        expect_identical(names(half), half.names)
        expect_identical(half$n, 8192L)
        expect_identical(half$u1, log(0.5))
        expect_identical(half$u0, half$u1 - (half$n - 1) * half$h)
        expect_equal(exp(half$u0), 1e-30)
        for (coefficient in paste0("c", 0:5)) {
          expect_true(is.double(half[[coefficient]]))
          expect_length(half[[coefficient]], half$n - 1L)
          expect_true(all(is.finite(half[[coefficient]])))
        }
      }
    }
  }
})

test_that("each gold's table file in data-raw/ is its FIXED_QTABLE entry", {
  # data-raw/ is build-ignored, so R CMD check runs without it.
  data.raw <- test_path("..", "..", "data-raw")
  skip_if_not(dir.exists(data.raw), "data-raw/ is absent")

  expect_setequal(
    list.files(data.raw, pattern = "^fixed_qtable_.*[.]rds$"),
    paste0("fixed_qtable_", names(FIXED_QTABLE), ".rds")
  )
  for (gold in names(FIXED_QTABLE)) {
    path <- file.path(data.raw, paste0("fixed_qtable_", gold, ".rds"))
    expect_true(identical(readRDS(path), FIXED_QTABLE[[gold]]))
  }
})

test_that("the gold quantile table reproduces stats::qbeta for every gold set", {
  for (gold in names(FIXED_GOLD)) {
    qtable <- FIXED_QTABLE[[gold]]

    for (component in c("U", "M")) {
      ctab <- qtable[[component]]
      for (lower.tail in c(TRUE, FALSE)) {
        p <- sweep_p(ctab$ab)
        reference <- stats::qbeta(p, ctab$a, ctab$b, lower.tail = lower.tail)
        table <- fixed_qbeta_table(p, ctab, lower.tail)
        expect_lte(max(abs(table - reference)), 1e-14)
        # p = 0 and 1 are exact ends; p at or below p.min is Rmath's own.
        ends <- p == 0 | p == 1 | log(p) <= ctab$ab$u0
        expect_identical(table[ends], reference[ends])
        # The ends of the sweep are the right ends.
        expect_identical(table[p == 0], if (lower.tail) 0 else 1)
        expect_identical(table[p == 1], if (lower.tail) 1 else 0)
      }
    }
  }
})

test_that("the table evaluated at its nodes returns the node values", {
  # Build and evaluation must place the nodes with the same arithmetic;
  # nodes from seq() reconstructed as u0 + i * h put a 1e-12 floor under
  # the interpolant (dev/qtable-exp.R).
  qtable <- FIXED_QTABLE[["horvath"]]
  for (component in c("U", "M")) {
    for (name in c("ab", "ba")) {
      half <- qtable[[component]][[name]]
      u <- half$u1 - ((half$n - 1):0) * half$h
      x <- stats::qbeta(u, half$a, half$b, log.p = TRUE)
      interior <- u > half$u0
      expect_lte(
        max(abs(fixed_eval_half(u[interior], half) - x[interior])),
        1e-14
      )
    }
  }
})

test_that("an element the table cannot place falls back to stats::qbeta", {
  qtable <- FIXED_QTABLE[["horvath"]]
  half <- qtable$U$ab
  p <- exp(half$u1 - (half$n - 4000.5) * half$h) # inside interval 4000
  # p <= 1/2 reads the (a, b) half on the lower tail and the (b, a) half on
  # the upper tail.
  for (lower.tail in c(TRUE, FALSE)) {
    ctab <- qtable$U
    name <- if (lower.tail) "ab" else "ba"
    ctab[[name]]$c0[4000L] <- Inf
    expect_identical(
      fixed_qbeta_table(p, ctab, lower.tail = lower.tail),
      stats::qbeta(p, ctab$a, ctab$b, lower.tail = lower.tail)
    )
  }
})

test_that("bmiq_calibration_fixed reproduces bmiq_calibration on real data", {
  datM <- epic_measured()

  general <- general_all_values(datM)
  fixed <- bmiq_calibration_fixed(datM)

  expect_true(all(general$success))
  expect_true(all(general$h.applied))
  expect_within_pin(fixed, general$calibrated)
})

test_that("the table reproduces stats::qbeta on the real gold-side inputs", {
  # Every gold-side quantile the six EPIC samples ask for, rebuilt from the
  # general function's diagnostics and compared with the table directly.
  datM <- epic_measured()
  gold <- betanorm:::FIXED_GOLD[["horvath"]]
  qtable <- FIXED_QTABLE[["horvath"]]
  result <- general_all_values(datM, debug = TRUE)

  evaluations <- 0L
  for (i in seq_len(nrow(datM))) {
    d <- result$diagnostics$samples[[i]]
    x <- unname(datM[i, ])
    class <- fixed_class_by_thresh(x, d$thresholds)
    mu <- d$component_means
    selU <- which(class == 1L)
    selM <- which(class == 3L)
    sites <- list(
      list(idx = selU[x[selU] <= mu[1L]], k = 1L, lower.tail = TRUE),
      list(idx = selU[x[selU] > mu[1L]], k = 1L, lower.tail = FALSE),
      list(idx = selM[x[selM] >= mu[3L]], k = 3L, lower.tail = FALSE)
    )
    for (site in sites) {
      k <- site$k
      p <- stats::pbeta(
        x[site$idx],
        d$component_a[k],
        d$component_b[k],
        lower.tail = site$lower.tail
      )
      reference <- stats::qbeta(
        p,
        gold$a[k],
        gold$b[k],
        lower.tail = site$lower.tail
      )
      ctab <- qtable[[if (k == 1L) "U" else "M"]]
      expect_lte(
        max(abs(fixed_qbeta_table(p, ctab, site$lower.tail) - reference)),
        1e-14
      )
      evaluations <- evaluations + length(p)
    }
  }
  expect_gt(evaluations, 50000L)
})

test_that("bmiq_calibration_fixed reproduces bmiq_calibration across blocks", {
  # sample.block.size is 8, so 11 samples exercises a full block plus a
  # partial one.
  datM <- synthetic_samples(11L)

  general <- general_all_values(datM)
  fixed <- bmiq_calibration_fixed(datM)

  expect_true(all(general$success))
  expect_within_pin(fixed, general$calibrated)
})

test_that("a stored table not built for the gold stops the call", {
  # FIXED_QTABLE travels apart from the R file when vendored, so a table
  # from another commit must not reach the quantile map.
  inputs <- make_calibration_inputs(n_probes = 600L, seed = 3L)
  stored <- FIXED_QTABLE
  refused <- function(table) {
    local_mocked_bindings(FIXED_QTABLE = table, .package = "betanorm")
    expect_error(
      bmiq_calibration_fixed(inputs$datM),
      regexp = "FIXED_QTABLE does not match FIXED_GOLD"
    )
  }

  bad <- stored
  bad$horvath$U <- stored$horvath$M
  bad$horvath$M <- stored$horvath$U
  refused(bad)

  bad <- stored
  bad$horvath$M$ab <- stored$horvath$M$ba
  bad$horvath$M$ba <- stored$horvath$M$ab
  refused(bad)

  bad <- stored
  bad$horvath$U$a <- stored$horvath$U$a * (1 + 1e-15)
  refused(bad)
  bad <- stored
  bad$horvath$M$ba$b <- stored$horvath$M$ba$b * (1 + 1e-15)
  refused(bad)

  refused(list(other = stored$horvath))
})

test_that("a skipped H step matches bmiq_calibration's optional policy", {
  # Push the methylated-component mean above every observation in both
  # fitters so the upper-M tail is empty and H fails (see
  # test-ordering-h-policy.R).
  real_fixed <- fixed_fit_mixture
  local_mocked_bindings(
    fit_mixture = mock_h_failure(fit_mixture),
    fixed_fit_mixture = function(beta, thresholds) {
      fit <- real_fixed(beta, thresholds)
      fit$mu <- c(0.15, 0.5, 1.5)
      fit
    },
    .package = "betanorm"
  )

  inputs <- make_calibration_inputs(n_probes = 1500L, seed = 11L)

  general <- suppressWarnings(general_all_values(inputs$datM))
  expect_false(general$h.applied[1L])

  expect_warning(
    fixed <- bmiq_calibration_fixed(inputs$datM),
    regexp = "^Some samples were calibrated in part, their intermediate values left as they were\\.$"
  )
  expect_within_pin(fixed, general$calibrated)
})

test_that("bmiq_calibration_fixed validates its inputs", {
  inputs <- make_calibration_inputs(n_probes = 600L, seed = 3L)
  datM <- inputs$datM

  expect_error(bmiq_calibration_fixed(datM, gold = "other"), regexp = "gold")
  expect_error(bmiq_calibration_fixed(datM, gold = 1L), regexp = "gold")
  expect_error(
    bmiq_calibration_fixed(matrix(1L, nrow = 1L, ncol = 10L)),
    regexp = "double"
  )
  expect_error(
    bmiq_calibration_fixed(as.numeric(datM)),
    regexp = "matrix"
  )
  expect_error(
    bmiq_calibration_fixed(datM[, 1:5, drop = FALSE]),
    regexp = "Too few probes"
  )

  expect_silent(bmiq_calibration_fixed(datM))
})

test_that("only a run of more than 200 samples reports progress, on one line", {
  # Every sample fails at its first step (no value below 0.4), so this is
  # quick; failed samples count toward progress.
  datM <- matrix(0.5, nrow = 201L, ncol = 6L)
  expect_warning(
    messages <- capture_messages(out <- bmiq_calibration_fixed(datM)),
    regexp = "NA"
  )
  expect_identical(
    paste(messages, collapse = ""),
    "Calibrating 201 samples: 20%...40%...60%...80%...100%\n"
  )
  expect_true(all(is.na(out)))

  messages <- capture_messages(
    suppressWarnings(bmiq_calibration_fixed(datM[1:200, ]))
  )
  expect_length(messages, 0L)
})

test_that("a failed sample becomes an NA row and the others proceed", {
  inputs <- make_calibration_inputs(n_probes = 600L, seed = 3L)
  good <- inputs$datM[1L, ]
  # No value below the lower mode window: the sample mode cannot be
  # estimated.
  bad <- pmax(good, 0.45)
  datM <- rbind(s1 = good, s2 = bad, s3 = good)

  warnings <- capture_warnings(fixed <- bmiq_calibration_fixed(datM))
  expect_identical(warnings, failed_samples_warning)
  expect_identical(dimnames(fixed), dimnames(datM))
  expect_true(all(is.na(fixed[2L, ])))
  expect_false(anyNA(fixed[c(1L, 3L), ]))

  general <- suppressWarnings(general_all_values(
    datM,
    on.sample.error = "continue"
  ))
  expect_identical(general$success, c(TRUE, FALSE, TRUE))
  expect_within_pin(fixed, general$calibrated)
})

test_that("a sample holding a bad value fails alone", {
  # At the last probe, so that the whole row has to be read.
  expect_sample_fails_alone(synthetic_samples(3L), row = 2L, probe = 1200L)
})

test_that("a sample holding a bad value fails alone across the block boundary", {
  # sample.block.size is 8, so sample 10 of 11 is in the second block, and the
  # samples before and after it have to land in their own rows.
  expect_sample_fails_alone(synthetic_samples(11L), row = 10L, probe = 1L)
})

test_that("a bad sample is failed before any of its work", {
  # Estimating the two sample modes is the first work a sample does. -1e-12 is
  # the case that would otherwise run to the end: nothing else refuses it.
  datM <- synthetic_samples(3L)
  datM[2L, 1L] <- -1e-12

  real_mode <- fixed_estimate_mode
  calls <- 0L
  local_mocked_bindings(
    fixed_estimate_mode = function(x) {
      calls <<- calls + 1L
      real_mode(x)
    },
    .package = "betanorm"
  )
  out <- suppressWarnings(bmiq_calibration_fixed(datM))

  expect_true(all(is.na(out[2L, ])))
  expect_identical(calls, 4L)
})

test_that("a sample holding an exact 0 and an exact 1 is calibrated", {
  # Probe 1 is in the unmethylated class and the last in the methylated one
  # (simulate_beta_mixture), so they run through both tails of the map.
  clean <- bmiq_calibration_fixed(synthetic_samples(3L))
  datM <- synthetic_samples(3L)
  datM[2L, c(1L, 1200L)] <- c(0, 1)

  expect_silent(out <- bmiq_calibration_fixed(datM))
  expect_false(anyNA(out))
  expect_identical(unname(out[2L, c(1L, 1200L)]), c(0, 1))
  expect_identical(out[-2L, ], clean[-2L, ])
  # The general function accepts this input, so the contract holds on it.
  expect_within_pin(out, general_all_values(datM)$calibrated)
})

test_that("every sample bad gives all NA, one warning and no error", {
  # Two blocks, every bad value in use.
  datM <- synthetic_samples(9L)
  for (i in seq_len(nrow(datM))) {
    datM[i, i] <- bad_values[[(i - 1L) %% length(bad_values) + 1L]]
  }
  warnings <- capture_warnings(out <- bmiq_calibration_fixed(datM))
  expect_identical(warnings, failed_samples_warning)
  expect_identical(
    out,
    matrix(
      NA_real_,
      nrow = nrow(datM),
      ncol = ncol(datM),
      dimnames = dimnames(datM)
    )
  )

  # One probe missing in every sample, as a filtered probe would be.
  datM <- synthetic_samples(3L)
  datM[, 100L] <- NA_real_
  warnings <- capture_warnings(out <- bmiq_calibration_fixed(datM))
  expect_identical(warnings, failed_samples_warning)
  expect_true(all(is.na(out)))

  # An all-missing numeric matrix is the same case. An integer one is coerced
  # by the gather.
  for (all.na in list(NA_real_, NA_integer_)) {
    datM <- matrix(all.na, nrow = 2L, ncol = 6L)
    warnings <- capture_warnings(out <- bmiq_calibration_fixed(datM))
    expect_identical(warnings, failed_samples_warning, info = typeof(all.na))
    expect_identical(
      out,
      matrix(NA_real_, nrow = 2L, ncol = 6L),
      info = typeof(all.na)
    )
  }
})

test_that("a matrix not stored as numeric stops the call at the front", {
  # checkmate skips its mode check when every value is missing, so these pass
  # it. Past this stop a character or list matrix reaches Rcpp's conversion
  # error at the gather, which ends the process in a -O0 build, under
  # devtools::test() too.
  for (datM in list(
    matrix(NA, nrow = 2L, ncol = 6L),
    matrix(NA_character_, nrow = 2L, ncol = 6L),
    matrix(NA_complex_, nrow = 2L, ncol = 6L),
    matrix(list(NULL), nrow = 2L, ncol = 6L)
  )) {
    expect_error(
      bmiq_calibration_fixed(datM),
      regexp = "^datM must be numeric\\.$",
      info = typeof(datM)
    )
  }

  # A matrix holding any value is refused by checkmate, as before.
  expect_error(
    bmiq_calibration_fixed(matrix("a", nrow = 2L, ncol = 6L)),
    regexp = "double"
  )
  expect_error(
    bmiq_calibration_fixed(matrix(TRUE, nrow = 2L, ncol = 6L)),
    regexp = "double"
  )
})

test_that("too few probes still stops when the matrix also holds a bad value", {
  # About the call, not a sample: no row could be calibrated at this width.
  datM <- matrix(0.5, nrow = 2L, ncol = 5L)
  for (name in names(bad_values)) {
    bad <- datM
    bad[1L, 1L] <- bad_values[[name]]
    expect_error(
      bmiq_calibration_fixed(bad),
      regexp = "^Too few probes\\.$",
      info = name
    )
  }
})

test_that("bmiq_calibration stops on the values the fixed function makes NA rows of", {
  # Why the contract with the general function ends at input it accepts.
  gold <- horvath_fit()
  datM <- synthetic_samples(3L)
  for (name in names(bad_values)) {
    bad <- datM
    bad[2L, 1200L] <- bad_values[[name]]
    expect_error(
      bmiq_calibration(
        bad,
        gold,
        on.sample.error = "continue",
        verbose = FALSE
      ),
      regexp = "datM",
      info = name
    )
  }
})

test_that("neither warning of bmiq_calibration_fixed carries a call", {
  # A downstream package's users see these warnings as they stand, and the
  # internal call means nothing to them. One failed sample, and one whose
  # upper-M tail is empty so that H fails, give both warnings in one run.
  real_fixed <- fixed_fit_mixture
  local_mocked_bindings(
    fixed_fit_mixture = function(beta, thresholds) {
      fit <- real_fixed(beta, thresholds)
      fit$mu <- c(0.15, 0.5, 1.5)
      fit
    },
    .package = "betanorm"
  )
  good <- make_calibration_inputs(n_probes = 1500L, seed = 11L)$datM[1L, ]
  datM <- rbind(s1 = good, s2 = pmax(good, 0.45))

  seen <- list()
  withCallingHandlers(
    bmiq_calibration_fixed(datM),
    warning = function(w) {
      seen[[length(seen) + 1L]] <<- w
      invokeRestart("muffleWarning")
    }
  )
  expect_identical(lapply(seen, conditionCall), list(NULL, NULL))
})

test_that("bmiq_calibration_fixed draws no random number", {
  # A draw would create a .Random.seed in a session with none, or move the
  # one that was set.
  inputs <- make_calibration_inputs(n_probes = 500L, seed = 3L)

  expect_random_state_kept(function() {
    bmiq_calibration_fixed(inputs$datM)
  })
})
