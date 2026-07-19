test_that("bmiq_calibration rejects out-of-range beta values", {
  inputs <- make_calibration_inputs(n_probes = 200L, seed = 1L)
  bad <- inputs$datM
  bad[1L, 1L] <- 1.05

  expect_error(
    bmiq_calibration(
      datM = bad,
      goldstandard.beta = inputs$gold,
      nfit = 200L,
      verbose = FALSE
    ),
    regexp = "\\[0, 1\\]|must have all values"
  )
})

test_that("bmiq_calibration rejects non-finite beta values", {
  inputs <- make_calibration_inputs(n_probes = 200L, seed = 2L)
  bad <- inputs$datM
  bad[1L, 1L] <- Inf

  expect_error(
    bmiq_calibration(
      datM = bad,
      goldstandard.beta = inputs$gold,
      nfit = 200L,
      verbose = FALSE
    ),
    regexp = "finite|Inf|\\[0, 1\\]"
  )
})
