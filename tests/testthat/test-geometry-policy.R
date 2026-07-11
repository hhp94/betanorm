test_that("clipBetaForFit matches legacy blc adaptive endpoints", {
  y <- c(0, 0.1, 0.5, 0.9, 1)
  out <- clipBetaForFit(y)

  expect_equal(length(out), length(y))
  expect_true(all(out > 0 & out < 1))
  # Exact blc: pmax(Y, min(Y[Y>0])/2), pmin(Y, 1 - (1 - max(Y[Y<1]))/2)
  expect_equal(out[1L], 0.1 / 2)
  expect_equal(out[5L], 1 - (1 - 0.9) / 2)
  expect_equal(out[2:4], y[2:4])
  # Materially larger than machine-eps clamping.
  expect_gt(out[1L], 1e-10)
})

test_that("thresholdsFromDensityCrossings solve weighted Beta log-score ties", {
  # Two well-separated components: equal weight, same concentration scale.
  # Crossing of eta * dbeta should sit between the means.
  a <- c(2, 8)
  b <- c(8, 2)
  eta <- c(0.5, 0.5)
  means <- a / (a + b)
  thr <- thresholdsFromDensityCrossings(
    a = a,
    b = b,
    eta = eta,
    component.means = means,
    context = "test"
  )
  expect_length(thr, 1L)
  expect_gt(thr, means[1L])
  expect_lt(thr, means[2L])

  # At the root, weighted log-densities agree.
  s1 <- log(eta[1L]) + dbeta(thr, a[1L], b[1L], log = TRUE)
  s2 <- log(eta[2L]) + dbeta(thr, a[2L], b[2L], log = TRUE)
  expect_equal(s1, s2, tolerance = 1e-6)
})

test_that("thresholdsFromDensityCrossings handles three components", {
  a <- c(2, 3, 10)
  b <- c(10, 3, 2)
  eta <- c(0.4, 0.2, 0.4)
  means <- a / (a + b)
  thr <- thresholdsFromDensityCrossings(
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
})

test_that("thresholdsFromDensityCrossings fails when means are not separated", {
  expect_error(
    thresholdsFromDensityCrossings(
      a = c(2, 2),
      b = c(2, 2),
      eta = c(0.5, 0.5),
      component.means = c(0.5, 0.5),
      context = "collapsed"
    ),
    regexp = "not separated|crossing"
  )
})

test_that("nL = 3 density thresholds classify and H stays continuous", {
  inputs <- make_calibration_inputs(n_probes = 2500L, nL = 3L, seed = 21L)

  result <- BMIQcalibration(
    datM = inputs$datM,
    goldstandard.beta = inputs$gold,
    nL = 3L,
    nfit = 2500L,
    niter = 10L,
    verbose = FALSE,
    debug = TRUE
  )

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
  thr2 <- thresholdsFromDensityCrossings(
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
