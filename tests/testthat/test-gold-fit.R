# The gold standard is consumed purely distributionally: bmiq_calibration()
# accepts either a raw gold vector (fitted internally) or a prefitted
# bmiq_gold_fit object, and the two paths must be interchangeable.

test_that("a prefitted gold standard reproduces the vector path exactly", {
  inputs <- make_calibration_inputs(n_probes = 1200L, nL = 3L, seed = 42L)

  from_vector <- bmiq_calibration(
    datM = inputs$datM,
    goldstandard.beta = inputs$gold,
    nfit = 1200L,
    verbose = FALSE
  )

  gold <- bmiq_gold_fit(inputs$gold, nfit = 1200L, verbose = FALSE)
  from_prefit <- bmiq_calibration(
    datM = inputs$datM,
    goldstandard.beta = gold,
    nfit = 1200L,
    verbose = FALSE
  )

  expect_s3_class(gold, "bmiq_gold_fit")
  expect_identical(as.matrix(from_prefit), as.matrix(from_vector))
})

test_that("the gold vector does not need to match ncol(datM)", {
  inputs <- make_calibration_inputs(n_probes = 1200L, nL = 3L, seed = 42L)
  gold_short <- simulate_beta_mixture(800L, nL = 3L, seed = 7L)

  result <- bmiq_calibration(
    datM = inputs$datM,
    goldstandard.beta = gold_short,
    nfit = 800L,
    verbose = FALSE
  )
  expect_true(all(result$success))
})

test_that("a prefit with a different nL is rejected", {
  inputs <- make_calibration_inputs(n_probes = 1200L, nL = 3L, seed = 42L)
  gold <- bmiq_gold_fit(inputs$gold, nL = 3L, nfit = 1200L, verbose = FALSE)

  expect_error(
    bmiq_calibration(
      datM = inputs$datM,
      goldstandard.beta = gold,
      nL = 2L,
      nfit = 1200L,
      verbose = FALSE
    ),
    regexp = "fitted with nL = 3"
  )
})

test_that("bmiq_gold_fit validates its inputs", {
  gold <- simulate_beta_mixture(600L, nL = 3L, seed = 3L)

  expect_error(bmiq_gold_fit(gold, nL = 4L, verbose = FALSE), regexp = "nL")
  expect_error(
    bmiq_gold_fit(c(gold, NA_real_), verbose = FALSE),
    regexp = "missing|finite"
  )
  expect_error(
    bmiq_gold_fit(gold, th1.v = c(0.9, 0.2), verbose = FALSE),
    regexp = "increasing"
  )
})

test_that("bmiq_gold_fit attaches diagnostics only when debug = TRUE", {
  gold <- simulate_beta_mixture(900L, nL = 3L, seed = 8L)

  plain <- bmiq_gold_fit(gold, nfit = 900L, verbose = FALSE)
  expect_null(plain$diagnostics)

  loud <- bmiq_gold_fit(gold, nfit = 900L, verbose = FALSE, debug = TRUE)
  expect_type(loud$diagnostics, "list")
  expect_length(loud$diagnostics$component_means, 3L)
})
