test_that("density_thresholds handles three components", {
  a <- c(2, 3, 10)
  b <- c(10, 3, 2)
  eta <- c(0.4, 0.2, 0.4)
  means <- a / (a + b)
  thr <- density_thresholds(
    a = a,
    b = b,
    eta = eta,
    component.means = means,
    context = "test3"
  )
  expect_length(thr, 2L)
  expect_true(is.unsorted(thr, strictly = TRUE) == FALSE)
  expect_gt(thr[1L], means[1L])
  expect_lt(thr[1L], means[2L])
  expect_gt(thr[2L], means[2L])
  expect_lt(thr[2L], means[3L])

  # Weighted log-densities agree at each boundary.
  for (k in seq_along(thr)) {
    s_lo <- log(eta[k]) + dbeta(thr[k], a[k], b[k], log = TRUE)
    s_hi <- log(eta[k + 1L]) + dbeta(thr[k], a[k + 1L], b[k + 1L], log = TRUE)
    expect_equal(s_lo, s_hi, tolerance = 1e-6)
  }
})

test_that("density_thresholds solves two-component log-score ties", {
  a <- c(2, 8)
  b <- c(8, 2)
  eta <- c(0.5, 0.5)
  means <- a / (a + b)
  thr <- density_thresholds(
    a = a,
    b = b,
    eta = eta,
    component.means = means,
    context = "test2"
  )
  expect_length(thr, 1L)
  expect_gt(thr, means[1L])
  expect_lt(thr, means[2L])

  s1 <- log(eta[1L]) + dbeta(thr, a[1L], b[1L], log = TRUE)
  s2 <- log(eta[2L]) + dbeta(thr, a[2L], b[2L], log = TRUE)
  expect_equal(s1, s2, tolerance = 1e-6)
})

test_that("density_thresholds fails when means are not separated", {
  expect_error(
    density_thresholds(
      a = c(2, 2),
      b = c(2, 2),
      eta = c(0.5, 0.5),
      component.means = c(0.5, 0.5),
      context = "collapsed"
    ),
    regexp = "unordered|not separated|crossing"
  )
})

test_that("density_thresholds requires the lower-to-higher orientation", {
  # Two components that cross with the wrong orientation between their means
  # must be rejected rather than silently returning the wrong crossing.
  expect_error(
    density_thresholds(
      a = c(8, 2),
      b = c(2, 8),
      eta = c(0.5, 0.5),
      component.means = c(0.2, 0.8),
      context = "flipped"
    ),
    regexp = "lower-to-higher|crossing|unordered"
  )
})

test_that("density thresholds classify and H stays continuous", {
  inputs <- make_calibration_inputs(n_probes = 2500L, seed = 21L)

  result <- suppressWarnings(bmiq_calibration(
    datM = inputs$datM,
    goldstandard.beta = inputs$gold,
    nfit = 2500L,
    niter = 10L,
    verbose = FALSE,
    debug = TRUE
  ))

  expect_true(result$success[1L])
  diag <- result$diagnostics$samples[[1L]]
  thr <- diag$thresholds
  means <- diag$component_means
  a <- diag$component_a
  b <- diag$component_b
  eta <- diag$eta

  expect_length(thr, 2L)
  expect_length(means, 3L)

  # Thresholds are density crossings (recompute and compare).
  thr2 <- density_thresholds(
    a = a,
    b = b,
    eta = eta,
    component.means = means,
    context = "recompute"
  )
  expect_equal(thr, thr2, tolerance = 1e-8)

  # Weighted log-scores match at each boundary.
  for (k in seq_along(thr)) {
    s_lo <- log(eta[k]) + dbeta(thr[k], a[k], b[k], log = TRUE)
    s_hi <- log(eta[k + 1L]) + dbeta(thr[k], a[k + 1L], b[k + 1L], log = TRUE)
    expect_equal(s_lo, s_hi, tolerance = 1e-5)
  }

  if (isTRUE(result$h.applied[1L])) {
    anchors <- diag$H_output_anchors
    expect_length(anchors, 2L)
    expect_true(is.finite(anchors[1L]) && is.finite(anchors[2L]))
    expect_gt(anchors[2L], anchors[1L])
  }
})
