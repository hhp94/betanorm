# bmiq_calibration_fixed() is bmiq_calibration() at its defaults with a
# built-in gold fit, and its output is pinned to it bit for bit. The stored
# gold numbers are pinned to the live bmiq_gold_fit() the same way.

horvath_vector <- function() {
  betanorm::horvath_goldstandard[["goldstandard2"]]
}

# GPL21145_sample without the 1062 gold-standard probes the EPIC array does
# not measure (zero in every sample).
epic_measured <- function() {
  datM <- betanorm::GPL21145_sample
  datM[, colSums(datM == 0) < nrow(datM)]
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

test_that("bmiq_calibration_fixed reproduces bmiq_calibration on real data", {
  datM <- epic_measured()

  general <- bmiq_calibration(datM, horvath_vector(), verbose = FALSE)
  fixed <- bmiq_calibration_fixed(datM, verbose = FALSE)

  expect_true(all(general$success))
  expect_true(all(general$h.applied))
  expect_identical(fixed, general$calibrated)
  expect_identical(dimnames(fixed), dimnames(datM))
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
  expect_identical(fixed, general$calibrated)
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
  expect_identical(fixed, general$calibrated)

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
  expect_identical(fixed, general$calibrated)

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
})
