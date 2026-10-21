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

test_that("quantile_norm, quantile_norm_fixed and the kernel leave obj and target as found", {
  # The kernel reads obj in place, not through a transposed copy, so a write
  # through the wrong pointer would land in the caller's matrix; target is
  # sorted in a copy, and is unsorted here so that a sort in place would
  # show. Each path runs: equal lengths, unequal lengths, one column. The
  # rounding gives ties. The saved copies are serialized and read back, as
  # above.
  set.seed(31L)
  wide <- matrix(
    round(rnorm(7L * 50L), 1),
    nrow = 7L,
    dimnames = list(paste0("s", 1:7), paste0("cg", 1:50))
  )
  cases <- list(
    equal = list(obj = wide, target = rnorm(50L)),
    unequal = list(obj = wide, target = rnorm(17L)),
    one.column = list(obj = wide[, 1L, drop = FALSE], target = rnorm(5L))
  )
  fns <- list(
    quantile_norm = quantile_norm,
    quantile_norm_fixed = quantile_norm_fixed,
    qnorm_target_rows_cpp = qnorm_target_rows_cpp
  )

  for (case in names(cases)) {
    obj <- cases[[case]][["obj"]]
    target <- cases[[case]][["target"]]
    obj.saved <- unserialize(serialize(obj, NULL))
    target.saved <- unserialize(serialize(target, NULL))
    for (fn in names(fns)) {
      info <- paste(case, fn)
      out <- fns[[fn]](obj, target)
      expect_identical(obj, obj.saved, info = info)
      expect_identical(target, target.saved, info = info)
      # And the result is a new matrix, not obj handed back.
      expect_false(identical(out, obj.saved), info = info)
    }
  }
})
