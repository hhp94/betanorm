# bmiq_calibration_fixed() is bmiq_calibration() at its defaults with a
# built-in gold fit, and its output is pinned to it within 1e-14: the gold
# side of its quantile map is a precomputed table (bmiq_gold_qtable()), not
# stats::qbeta(), so the two are not bit-identical. The stored gold numbers
# are pinned to the live bmiq_gold_fit() bit for bit, and the table is
# pinned to stats::qbeta() directly.

horvath_vector <- function() {
  betanorm::horvath_goldstandard[["goldstandard2"]]
}

# GPL21145_sample without the 1062 gold-standard probes the EPIC array does
# not measure (zero in every sample).
epic_measured <- function() {
  datM <- betanorm::GPL21145_sample
  datM[, colSums(datM == 0) < nrow(datM)]
}

# The numerics contract of the fixed function: same NA pattern and dimnames
# as the general function, every other element within 1e-14 (absolute).
expect_within_pin <- function(fixed, general) {
  expect_identical(dim(fixed), dim(general))
  expect_identical(dimnames(fixed), dimnames(general))
  expect_identical(is.na(fixed), is.na(general))
  expect_lte(max(abs(fixed - general), na.rm = TRUE), 1e-14)
}

# Dense p per (component, tail), after dev/qtable-exp.R without its random
# points: the grid nodes, three series of 20000 points spaced about 0.41 h
# apart from different starting offsets (about seven points per interval),
# their mirrors above 1/2, and the boundary cases. p is built first and
# the log taken from it, as in production.
sweep_p <- function(half) {
  eps <- .Machine$double.eps
  p.min <- FIXED_QTABLE_P_MIN
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

test_that("the stored Horvath gold fit is the live bmiq_gold_fit()", {
  live <- bmiq_gold_fit(horvath_vector(), nL = 3L, verbose = FALSE)
  stored <- betanorm:::FIXED_GOLD[["horvath"]]

  expect_identical(names(betanorm:::FIXED_GOLD), "horvath")
  expect_identical(stored$a, live$a)
  expect_identical(stored$b, live$b)
  expect_identical(stored$thresholds, live$thresholds)
  expect_identical(stored$unmethylated.mode, live$unmethylated.mode)
  expect_identical(stored$methylated.mode, live$methylated.mode)
})

test_that("the gold quantile table reproduces stats::qbeta for every gold set", {
  for (gold in names(FIXED_GOLD)) {
    qtable <- bmiq_gold_qtable(gold)
    expect_identical(names(qtable), c("version", "gold", "U", "M"))
    expect_identical(qtable$version, FIXED_QTABLE_VERSION)
    expect_identical(qtable$gold, gold)

    for (component in c("U", "M")) {
      ctab <- qtable[[component]]
      for (lower.tail in c(TRUE, FALSE)) {
        p <- sweep_p(ctab$ab)
        reference <- stats::qbeta(p, ctab$a, ctab$b, lower.tail = lower.tail)
        table <- fixed_qbeta_table(p, ctab, lower.tail)
        expect_lte(max(abs(table - reference)), 1e-14)
        # p = 0 and 1 are exact ends; p at or below p.min is Rmath's own.
        ends <- p == 0 | p == 1 | p <= FIXED_QTABLE_P_MIN
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
  qtable <- bmiq_gold_qtable()
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
  qtable <- bmiq_gold_qtable()
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

  general <- bmiq_calibration(datM, horvath_vector(), verbose = FALSE)
  fixed <- bmiq_calibration_fixed(datM, verbose = FALSE)

  expect_true(all(general$success))
  expect_true(all(general$h.applied))
  expect_within_pin(fixed, general$calibrated)
})

test_that("the table reproduces stats::qbeta on the real gold-side inputs", {
  # Every gold-side quantile the six EPIC samples ask for, rebuilt from the
  # general function's diagnostics and compared with the table directly.
  datM <- epic_measured()
  gold <- betanorm:::FIXED_GOLD[["horvath"]]
  qtable <- bmiq_gold_qtable()
  result <- bmiq_calibration(
    datM,
    horvath_vector(),
    debug = TRUE,
    verbose = FALSE
  )

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
  n_probes <- 1200L
  n_samples <- 11L
  datM <- t(vapply(
    seq_len(n_samples),
    function(i) {
      clip01(simulate_beta_mixture(n_probes, nL = 3L, seed = 500L + i))
    },
    numeric(n_probes)
  ))
  rownames(datM) <- paste0("s", seq_len(n_samples))

  general <- bmiq_calibration(datM, horvath_vector(), verbose = FALSE)
  fixed <- bmiq_calibration_fixed(datM, verbose = FALSE)

  expect_true(all(general$success))
  expect_within_pin(fixed, general$calibrated)
})

test_that("a passed table equals a built one and gives identical output", {
  inputs <- make_calibration_inputs(n_probes = 600L, seed = 3L)
  qtable <- bmiq_gold_qtable()

  expect_identical(qtable, bmiq_gold_qtable("horvath"))
  expect_identical(
    bmiq_calibration_fixed(inputs$datM, verbose = FALSE, qtable = qtable),
    bmiq_calibration_fixed(inputs$datM, verbose = FALSE)
  )
  # Positional calls written before qtable existed still work.
  expect_identical(
    bmiq_calibration_fixed(inputs$datM, "horvath", FALSE),
    bmiq_calibration_fixed(inputs$datM, verbose = FALSE, qtable = qtable)
  )
})

test_that("a malformed table is refused", {
  inputs <- make_calibration_inputs(n_probes = 600L, seed = 3L)
  datM <- inputs$datM
  qtable <- bmiq_gold_qtable()
  refused <- function(bad, regexp) {
    expect_error(
      bmiq_calibration_fixed(datM, verbose = FALSE, qtable = bad),
      regexp = paste0("qtable ", regexp, ".*bmiq_gold_qtable\\(\"horvath\"\\)")
    )
  }

  expect_error(bmiq_gold_qtable("other"), regexp = "gold")
  expect_error(bmiq_gold_qtable(1L), regexp = "gold")

  refused(list(1), "is not a list")
  refused(unclass(as.matrix(1)), "is not a list")

  bad <- qtable
  bad$gold <- "other"
  refused(bad, "was built for gold \"other\", not \"horvath\"")
  bad$gold <- sum
  refused(bad, "was not built for gold \"horvath\"")

  bad <- qtable
  bad$version <- 0L
  refused(bad, "was built by another version")
  bad$version <- 1
  refused(bad, "was built by another version")

  bad <- qtable
  tmp <- bad$U
  bad$U <- bad$M
  bad$M <- tmp
  refused(bad, "\\$U was built for other gold shapes")

  bad <- qtable
  bad$M$ab$a <- bad$M$ab$a * (1 + 1e-15)
  refused(bad, "\\$M\\$ab was built for other shapes")

  bad <- qtable
  tmp <- bad$M$ab
  bad$M$ab <- bad$M$ba
  bad$M$ba <- tmp
  refused(bad, "\\$M\\$ab was built for other shapes")

  bad <- qtable
  bad$U$ab$h <- bad$U$ab$h * (1 + 1e-15)
  refused(bad, "\\$U\\$ab was built on another grid")
  bad <- qtable
  bad$U$ba$n <- 4096L
  refused(bad, "\\$U\\$ba was built on another grid")

  bad <- qtable
  bad$M$ba$c3 <- bad$M$ba$c3[-1L]
  refused(bad, "\\$M\\$ba\\$c3 is not a finite double vector")
  bad <- qtable
  bad$U$ab$c5[10L] <- NaN
  refused(bad, "\\$U\\$ab\\$c5 is not a finite double vector")
  bad <- qtable
  bad$U$ab$c0 <- as.integer(bad$U$ab$c0)
  refused(bad, "\\$U\\$ab\\$c0 is not a finite double vector")

  bad <- qtable
  bad$M$ba$c5 <- NULL
  refused(bad, "\\$M\\$ba is not a half-table")
  bad <- qtable
  bad$U$extra <- 1
  refused(bad, "\\$U is not a component table")

  # Well-formed but wrong values: the structural checks pass, the spot
  # check against stats::qbeta does not.
  bad <- qtable
  bad$U$ab$c0 <- bad$U$ab$c0 + 1e-9
  refused(bad, "\\$U\\$ab does not reproduce stats::qbeta\\(\\)")
  bad <- qtable
  bad$M$ba$c3 <- 2 * bad$M$ba$c3
  refused(bad, "\\$M\\$ba does not reproduce stats::qbeta\\(\\)")
  bad <- qtable
  foreign <- fixed_make_half(3, 30)
  for (name in paste0("c", 0:5)) bad$M$ab[[name]] <- foreign[[name]]
  refused(bad, "\\$M\\$ab does not reproduce stats::qbeta\\(\\)")
})

test_that("a table saved to disk and read back is accepted", {
  inputs <- make_calibration_inputs(n_probes = 600L, seed = 3L)
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path))
  saveRDS(bmiq_gold_qtable(), path)
  expect_identical(
    bmiq_calibration_fixed(inputs$datM, verbose = FALSE, qtable = readRDS(path)),
    bmiq_calibration_fixed(inputs$datM, verbose = FALSE)
  )
})

test_that("the end-of-run warnings list at most ten samples", {
  expect_identical(fixed_list_samples(1:10, sep = ", "), paste(1:10, collapse = ", "))
  expect_identical(
    fixed_list_samples(1:12, sep = ", "),
    paste(c(1:10, "and 2 more"), collapse = ", ")
  )
  expect_identical(fixed_list_samples(character(), sep = "; "), "")
})

test_that("a skipped H step matches bmiq_calibration's optional policy", {
  # Push the methylated-component mean above every observation in both
  # fitters so the upper-M tail is empty and H fails (see
  # test-ordering-h-policy.R).
  real_fixed <- fixed_fit_mixture
  local_mocked_bindings(
    fit_mixture = mock_h_failure(fit_mixture),
    fixed_fit_mixture = function(beta, thresholds, fit.idx, context) {
      fit <- real_fixed(beta, thresholds, fit.idx, context)
      fit$mu <- c(0.15, 0.5, 1.5)
      fit
    },
    .package = "betanorm"
  )

  inputs <- make_calibration_inputs(n_probes = 1500L, seed = 11L)

  general <- suppressWarnings(bmiq_calibration(
    inputs$datM,
    horvath_vector(),
    verbose = FALSE
  ))
  expect_false(general$h.applied[1L])

  expect_warning(
    fixed <- bmiq_calibration_fixed(inputs$datM, verbose = FALSE),
    regexp = "without H.*: 1\\."
  )
  expect_within_pin(fixed, general$calibrated)

  messages <- capture_messages(
    suppressWarnings(bmiq_calibration_fixed(inputs$datM, verbose = TRUE))
  )
  expect_match(messages, "H skipped for sample 1", all = FALSE)
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
  bad <- datM
  bad[1L, 1L] <- NA_real_
  expect_error(bmiq_calibration_fixed(bad), regexp = "datM")
  bad[1L, 1L] <- 1.5
  expect_error(bmiq_calibration_fixed(bad), regexp = "datM")
  expect_error(
    bmiq_calibration_fixed(datM[, 1:5, drop = FALSE]),
    regexp = "Too few probes"
  )
  expect_error(bmiq_calibration_fixed(datM, verbose = NA), regexp = "verbose")

  expect_silent(bmiq_calibration_fixed(datM, verbose = FALSE))
  messages <- capture_messages(bmiq_calibration_fixed(datM))
  expect_match(messages, "Processing sample 1 of 1", all = FALSE)
})

test_that("a failed sample becomes an NA row and the others proceed", {
  inputs <- make_calibration_inputs(n_probes = 600L, seed = 3L)
  good <- inputs$datM[1L, ]
  # No value below the lower mode window: the sample mode cannot be
  # estimated.
  bad <- pmax(good, 0.45)
  datM <- rbind(s1 = good, s2 = bad, s3 = good)

  expect_warning(
    fixed <- bmiq_calibration_fixed(datM, verbose = FALSE),
    regexp = "1 sample\\(s\\) failed.*: 2 \\(Sample 2 values below 0.4 is empty"
  )
  expect_identical(dimnames(fixed), dimnames(datM))
  expect_true(all(is.na(fixed[2L, ])))
  expect_false(anyNA(fixed[c(1L, 3L), ]))

  general <- suppressWarnings(bmiq_calibration(
    datM,
    horvath_vector(),
    on.sample.error = "continue",
    verbose = FALSE
  ))
  expect_identical(general$success, c(TRUE, FALSE, TRUE))
  expect_within_pin(fixed, general$calibrated)

  messages <- capture_messages(
    suppressWarnings(bmiq_calibration_fixed(datM, verbose = TRUE))
  )
  expect_match(messages, "Sample 2 failed .*; replaced with NA", all = FALSE)

  # A failed sample is not also reported as calibrated without H.
  warnings <- capture_warnings(bmiq_calibration_fixed(datM, verbose = FALSE))
  expect_length(warnings, 1L)
})

test_that("bmiq_calibration_fixed leaves the random state as found", {
  inputs <- make_calibration_inputs(n_probes = 500L, seed = 3L)

  expect_random_state_kept(function() {
    bmiq_calibration_fixed(inputs$datM, verbose = FALSE)
  })
  expect_random_state_kept(bmiq_gold_qtable)
})
