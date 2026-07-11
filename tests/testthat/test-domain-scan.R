test_that("BMIQcalibration rejects out-of-range beta values", {
  inputs <- make_calibration_inputs(n_probes = 200L, nL = 3L, seed = 1L)
  bad <- inputs$datM
  bad[1L, 1L] <- 1.05

  expect_error(
    BMIQcalibration(
      datM = bad,
      goldstandard.beta = inputs$gold,
      nL = 3L,
      nfit = 200L,
      verbose = FALSE
    ),
    regexp = "\\[0, 1\\]|must have all values"
  )
})

test_that("BMIQcalibration rejects non-finite beta values", {
  inputs <- make_calibration_inputs(n_probes = 200L, nL = 3L, seed = 2L)
  bad <- inputs$datM
  bad[1L, 1L] <- Inf

  expect_error(
    BMIQcalibration(
      datM = bad,
      goldstandard.beta = inputs$gold,
      nL = 3L,
      nfit = 200L,
      verbose = FALSE
    ),
    regexp = "finite|Inf|\\[0, 1\\]"
  )
})

test_that("Boost special functions are usable on open-interval betas", {
  # digamma has a pole at 0: Boost should surface a domain error via Rcpp.
  digamma_zero <- tryCatch(
    boost_digamma_cpp(0),
    error = function(e) e
  )
  expect_true(
    inherits(digamma_zero, "error") ||
      !is.finite(as.numeric(digamma_zero)[1L])
  )

  z <- .Machine$double.eps
  expect_equal(
    boost_digamma_cpp(z),
    digamma(z),
    tolerance = 1e-10
  )

  q_open <- c(1e-12, 1 - 1e-12)
  p_open <- boost_pbeta_cpp(q_open, 2, 5, lower_tail = TRUE)
  expect_true(all(is.finite(p_open)))
  expect_true(all(p_open >= 0 & p_open <= 1))
  expect_equal(
    p_open,
    pbeta(q_open, 2, 5, lower.tail = TRUE),
    tolerance = 1e-10
  )

  p_safe <- c(1e-12, 1 - 1e-12)
  q_safe <- boost_qbeta_cpp(p_safe, 2, 5, lower_tail = TRUE)
  expect_true(all(is.finite(q_safe)))
  expect_true(all(q_safe > 0 & q_safe < 1))
  expect_equal(
    q_safe,
    qbeta(p_safe, 2, 5, lower.tail = TRUE),
    tolerance = 1e-8
  )
})
