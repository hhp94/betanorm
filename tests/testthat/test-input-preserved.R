# The calibration functions must leave datM as found. scatter_sample_block_cpp()
# writes into its destination in place, past R's copy-on-modify, so the output
# must be a freshly allocated matrix, never datM or an alias of it
# (`calibrated <- datM` would be one). The reference copy is serialized and
# read back: `saved <- datM` would be a second binding of the same vector,
# which an in-place write changes too.

test_that("bmiq_calibration and bmiq_calibration_fixed leave datM as found", {
  # 11 samples are a full block of 8 and a partial one. Sample 10 fails (no
  # value below the lower mode window), so the scatter writes at an offset
  # and the failed-sample paths run too.
  n_probes <- 1200L
  n_samples <- 11L
  datM <- t(vapply(
    seq_len(n_samples),
    function(i) {
      clip01(simulate_beta_mixture(n_probes, nL = 3L, seed = 500L + i))
    },
    numeric(n_probes)
  ))
  datM[10L, ] <- pmax(datM[10L, ], 0.45)
  dimnames(datM) <- list(
    paste0("s", seq_len(n_samples)),
    paste0("cg", seq_len(n_probes))
  )
  saved <- unserialize(serialize(datM, NULL))
  failed <- seq_len(n_samples) == 10L

  fixed <- suppressWarnings(bmiq_calibration_fixed(datM))
  expect_identical(unname(is.na(fixed[, 1L])), failed)
  expect_identical(datM, saved)

  gold <- bmiq_gold_fit(
    betanorm::horvath_goldstandard[["goldstandard2"]],
    verbose = FALSE
  )
  for (failed.sample in c("NA", "original")) {
    general <- suppressWarnings(bmiq_calibration(
      datM,
      gold,
      on.sample.error = "continue",
      failed.sample = failed.sample,
      verbose = FALSE
    ))
    expect_identical(general$success, !failed, info = failed.sample)
    expect_identical(datM, saved, info = failed.sample)
  }
})
