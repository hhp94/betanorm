# h.policy, on.sample.error and failed.sample take one exact value: an
# abbreviation that match.arg() used to complete is refused.

test_that("abbreviated option values are refused", {
  inputs <- make_calibration_inputs(n_probes = 400L, seed = 1L)
  calibrate <- function(...) {
    bmiq_calibration(
      datM = inputs$datM,
      goldstandard.beta = inputs$gold,
      nfit = 400L,
      verbose = FALSE,
      ...
    )
  }

  # The same call with full values runs, so each error below comes from the
  # abbreviation alone.
  expect_no_error(calibrate(
    h.policy = "require",
    on.sample.error = "continue",
    failed.sample = "original"
  ))

  expect_error(calibrate(h.policy = "opt"))
  expect_error(calibrate(on.sample.error = "cont"))
  expect_error(calibrate(failed.sample = "orig"))
})
