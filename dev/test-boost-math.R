test_that("Boost digamma and trigamma match base R on positive reals", {
  x <- c(0.1, 0.5, 1, 1.5, 2, 5, 10, 20, 50, 100)

  expect_equal(boost_digamma_cpp(x), digamma(x), tolerance = 1e-10)
  expect_equal(boost_trigamma_cpp(x), trigamma(x), tolerance = 1e-10)
})

test_that("Boost digamma/trigamma match R on Beta-shape grids", {
  # Domain used by the Newton Beta MLE (shapes above min_shape).
  shapes <- exp(seq(log(1e-4), log(1e3), length.out = 40))

  expect_equal(boost_digamma_cpp(shapes), digamma(shapes), tolerance = 1e-9)
  expect_equal(boost_trigamma_cpp(shapes), trigamma(shapes), tolerance = 1e-9)
})

test_that("Boost pbeta matches base R lower and upper tails", {
  q <- seq(0.01, 0.99, by = 0.01)
  shapes <- list(
    c(0.5, 0.5),
    c(1, 1),
    c(2, 5),
    c(5, 2),
    c(10, 10),
    c(0.8, 3.2)
  )

  for (ab in shapes) {
    a <- ab[1L]
    b <- ab[2L]
    expect_equal(
      boost_pbeta_cpp(q, a, b, lower_tail = TRUE),
      pbeta(q, a, b, lower.tail = TRUE),
      tolerance = 1e-10
    )
    expect_equal(
      boost_pbeta_cpp(q, a, b, lower_tail = FALSE),
      pbeta(q, a, b, lower.tail = FALSE),
      tolerance = 1e-10
    )
  }
})

test_that("Boost qbeta matches base R lower and upper tails", {
  p <- seq(0.01, 0.99, by = 0.01)
  shapes <- list(
    c(0.5, 0.5),
    c(1, 1),
    c(2, 5),
    c(5, 2),
    c(10, 10),
    c(0.8, 3.2)
  )

  for (ab in shapes) {
    a <- ab[1L]
    b <- ab[2L]
    expect_equal(
      boost_qbeta_cpp(p, a, b, lower_tail = TRUE),
      qbeta(p, a, b, lower.tail = TRUE),
      tolerance = 1e-8
    )
    expect_equal(
      boost_qbeta_cpp(p, a, b, lower_tail = FALSE),
      qbeta(p, a, b, lower.tail = FALSE),
      tolerance = 1e-8
    )
  }
})

test_that("Boost pbeta/qbeta round-trip agrees with R for BMIQ-like tails", {
  # Mimic U/M quantile normalization tails used in bmiq_calibration.
  set.seed(99)
  y <- pmin(0.999, pmax(0.001, rbeta(200, 2, 8)))
  a_sample <- 2.1
  b_sample <- 7.8
  a_gold <- 2.0
  b_gold <- 8.0

  p_lower <- pbeta(y, a_sample, b_sample, lower.tail = TRUE)
  p_upper <- pbeta(y, a_sample, b_sample, lower.tail = FALSE)

  boost_lower <- boost_qbeta_cpp(
    boost_pbeta_cpp(y, a_sample, b_sample, lower_tail = TRUE),
    a_gold,
    b_gold,
    lower_tail = TRUE
  )
  r_lower <- qbeta(p_lower, a_gold, b_gold, lower.tail = TRUE)

  boost_upper <- boost_qbeta_cpp(
    boost_pbeta_cpp(y, a_sample, b_sample, lower_tail = FALSE),
    a_gold,
    b_gold,
    lower_tail = FALSE
  )
  r_upper <- qbeta(p_upper, a_gold, b_gold, lower.tail = FALSE)

  expect_equal(boost_lower, r_lower, tolerance = 1e-8)
  expect_equal(boost_upper, r_upper, tolerance = 1e-8)
})
