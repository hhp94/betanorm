test_that("nL = 4 is rejected", {
  inputs <- make_calibration_inputs(n_probes = 400L, nL = 3L, seed = 1L)
  expect_error(
    BMIQcalibration(
      datM = inputs$datM,
      goldstandard.beta = inputs$gold,
      nL = 4L,
      nfit = 400L,
      verbose = FALSE
    ),
    regexp = "nL|2|3|upper"
  )
})

test_that("nL = 2 with dynamic defaults works and leaves h.applied as NA", {
  inputs <- make_calibration_inputs(n_probes = 1200L, nL = 2L, seed = 2L)

  result <- BMIQcalibration(
    datM = inputs$datM,
    goldstandard.beta = inputs$gold,
    nL = 2L,
    nfit = 1200L,
    niter = 8L,
    verbose = FALSE,
    debug = TRUE
  )

  expect_true(result$success[1L])
  expect_true(is.na(result$h.applied[1L]))
  expect_false(isTRUE(result$settings$doH))
  expect_equal(result$settings$nL, 2L)
  expect_equal(ncol(result$calibrated), ncol(inputs$datM))
  expect_true(all(is.finite(result$calibrated[1L, ])))
  expect_true(all(result$calibrated[1L, ] >= 0 & result$calibrated[1L, ] <= 1))
})

test_that("nL = 2, doH = TRUE fails clearly", {
  inputs <- make_calibration_inputs(n_probes = 400L, nL = 2L, seed = 3L)
  expect_error(
    BMIQcalibration(
      datM = inputs$datM,
      goldstandard.beta = inputs$gold,
      nL = 2L,
      doH = TRUE,
      nfit = 400L,
      verbose = FALSE
    ),
    regexp = "doH|nL = 2|no H"
  )
})

test_that("nL = 2 truncated map remaps both sides and is continuous at the cut", {
  inputs <- make_calibration_inputs(n_probes = 2000L, nL = 2L, seed = 4L)

  result <- BMIQcalibration(
    datM = inputs$datM,
    goldstandard.beta = inputs$gold,
    nL = 2L,
    nfit = 2000L,
    niter = 10L,
    verbose = FALSE,
    debug = TRUE
  )

  expect_true(result$success[1L])
  diag <- result$diagnostics$samples[[1L]]
  expect_true(!is.null(diag$tail_counts))
  expect_true(diag$tail_counts[["M_left"]] >= 0L)
  expect_true(diag$tail_counts[["M_right"]] >= 0L)

  original <- as.numeric(inputs$datM[1L, ])
  calibrated <- as.numeric(result$calibrated[1L, ])
  expect_true(any(abs(calibrated - original) > 1e-12))

  thr <- diag$thresholds
  expect_length(thr, 1L)
  t_s <- thr[1L]
  t_g <- diag$nl2_gold_threshold
  expect_true(is.finite(t_g))

  class <- ifelse(original <= t_s, 1L, 2L)
  u_idx <- which(class == 1L)
  m_idx <- which(class == 2L)
  expect_true(length(u_idx) > 0L && length(m_idx) > 0L)

  # Both sides of the mean are covered by the single truncated M map.
  mu_m <- diag$component_means[2L]
  m_left <- which(class == 2L & original < mu_m)
  m_right <- which(class == 2L & original > mu_m)
  if (length(m_left) > 0L) {
    expect_true(any(abs(calibrated[m_left] - original[m_left]) > 1e-12))
  }
  if (length(m_right) > 0L) {
    expect_true(any(abs(calibrated[m_right] - original[m_right]) > 1e-12))
  }

  # Continuity: no downward jump; edges meet at gold cut (within tol).
  gap <- min(calibrated[m_idx]) - max(calibrated[u_idx])
  expect_gte(gap, -1e-8)
  # Map of the cut from either side is the gold threshold.
  cut_u <- normalizeNL2Truncated(
    beta = t_s,
    class = 1L,
    sample.a = diag$component_a,
    sample.b = diag$component_b,
    gold.a = result$diagnostics$gold$component_a,
    gold.b = result$diagnostics$gold$component_b,
    sample.threshold = t_s,
    gold.threshold = t_g,
    context = "cut-U"
  )
  # Force class M at the cut to check the M-side formula joins at t_g.
  cut_m <- normalizeNL2Truncated(
    beta = t_s,
    class = 2L,
    sample.a = diag$component_a,
    sample.b = diag$component_b,
    gold.a = result$diagnostics$gold$component_a,
    gold.b = result$diagnostics$gold$component_b,
    sample.threshold = t_s,
    gold.threshold = t_g,
    context = "cut-M"
  )
  expect_equal(cut_u, t_g, tolerance = 1e-6)
  expect_equal(as.numeric(cut_m), t_g, tolerance = 1e-6)
})

test_that("normalizeNL2Truncated maps sample cut to gold cut from both sides", {
  # Symmetric Beta components; analytic check of the join.
  sample.a <- c(2, 5)
  sample.b <- c(5, 2)
  gold.a <- c(3, 6)
  gold.b <- c(6, 3)
  t_s <- 0.45
  t_g <- 0.55

  g_u <- normalizeNL2Truncated(
    beta = t_s,
    class = 1L,
    sample.a = sample.a,
    sample.b = sample.b,
    gold.a = gold.a,
    gold.b = gold.b,
    sample.threshold = t_s,
    gold.threshold = t_g
  )
  # M-side formula at the cut (class forced to 2) also yields t_g.
  g_m <- normalizeNL2Truncated(
    beta = t_s,
    class = 2L,
    sample.a = sample.a,
    sample.b = sample.b,
    gold.a = gold.a,
    gold.b = gold.b,
    sample.threshold = t_s,
    gold.threshold = t_g
  )
  expect_equal(g_u, t_g, tolerance = 1e-8)
  expect_equal(g_m, t_g, tolerance = 1e-8)

  # Interior U and M stay on the correct side of t_g.
  x_u <- 0.2
  x_m <- 0.8
  g <- normalizeNL2Truncated(
    beta = c(x_u, x_m),
    class = c(1L, 2L),
    sample.a = sample.a,
    sample.b = sample.b,
    gold.a = gold.a,
    gold.b = gold.b,
    sample.threshold = t_s,
    gold.threshold = t_g
  )
  expect_lte(g[1L], t_g + 1e-10)
  expect_gte(g[2L], t_g - 1e-10)
  expect_true(all(g > 0 & g < 1))
})

test_that("nL = 3 retains U/H/M behavior with dynamic defaults", {
  inputs <- make_calibration_inputs(n_probes = 2000L, nL = 3L, seed = 5L)

  result <- BMIQcalibration(
    datM = inputs$datM,
    goldstandard.beta = inputs$gold,
    nL = 3L,
    nfit = 2000L,
    niter = 8L,
    verbose = FALSE,
    debug = TRUE
  )

  expect_true(result$success[1L])
  expect_true(isTRUE(result$settings$doH))
  expect_true(result$h.applied[1L] %in% c(TRUE, FALSE))
  expect_true(all(is.finite(result$calibrated[1L, ])))

  gold_diag <- result$diagnostics$gold
  expect_length(gold_diag$component_means, 3L)
  expect_true(is.unsorted(gold_diag$component_means, strictly = TRUE) == FALSE)

  sample_diag <- result$diagnostics$samples[[1L]]
  expect_true(!is.null(sample_diag$parameter_criterion))
  expect_true(!is.null(sample_diag$loglik_criterion))
  expect_null(sample_diag$em_criterion)
})
