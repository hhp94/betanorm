make_em <- function(a, b, eta, mu, n = 6L) {
  nL <- length(mu)
  list(
    a = matrix(a, ncol = 1L),
    b = matrix(b, ncol = 1L),
    eta = eta,
    mu = matrix(mu, ncol = 1L),
    w = matrix(seq_len(n * nL), nrow = n, ncol = nL),
    fit_status = rep("converged", nL),
    fit_reason = paste0("reason", seq_len(nL))
  )
}

test_that("canonicalize_em_components sorts complete tuples by mean", {
  em <- make_em(
    a = c(8, 2, 5),
    b = c(2, 8, 5),
    eta = c(0.5, 0.3, 0.2),
    mu = c(0.8, 0.2, 0.5)
  )
  ord <- order(c(0.8, 0.2, 0.5)) # c(2, 3, 1)

  out <- canonicalize_em_components(em, "test")

  expect_equal(as.numeric(out$mu[, 1L]), c(0.2, 0.5, 0.8))
  expect_equal(as.numeric(out$a[, 1L]), c(2, 5, 8))
  expect_equal(as.numeric(out$b[, 1L]), c(8, 5, 2))
  expect_equal(out$eta, c(0.3, 0.2, 0.5))
  expect_equal(out$fit_reason, paste0("reason", ord))
  # Responsibility columns move with their component.
  expect_equal(out$w, em$w[, ord, drop = FALSE])
})

test_that("canonicalize_em_components rejects invalid mixture parameters", {
  bad <- make_em(
    a = c(-1, 2, 5),
    b = c(2, 8, 5),
    eta = c(0.5, 0.3, 0.2),
    mu = c(0.2, 0.5, 0.8)
  )
  expect_error(
    canonicalize_em_components(bad, "test"),
    regexp = "invalid mixture parameters"
  )
})

test_that("canonicalize_em_components rejects indistinguishable means", {
  tied <- make_em(
    a = c(2, 2),
    b = c(2, 2),
    eta = c(0.5, 0.5),
    mu = c(0.5, 0.5)
  )
  expect_error(
    canonicalize_em_components(tied, "test"),
    regexp = "indistinguishable"
  )
})

# Force H normalization to fail through a guard that survives canonicalization:
# push the methylated-component mean above every observation so there are no
# methylated probes above it (empty upper-M tail). Returns a mock function for
# fit_mixture; the caller must install it with local_mocked_bindings() in its
# own frame so the binding stays live for the test body.
mock_h_failure <- function(real) {
  # Force `real` now: if it stayed a lazy promise it would resolve to the
  # mocked binding (itself) once installed, causing infinite recursion.
  force(real)
  function(
    beta,
    thresholds,
    nL,
    nfit,
    niter,
    tol,
    beta.maxit,
    beta.score.tol,
    context,
    debug = FALSE,
    seed = 1L,
    fit.idx = NULL
  ) {
    fit <- real(
      beta = beta,
      thresholds = thresholds,
      nL = nL,
      nfit = nfit,
      niter = niter,
      tol = tol,
      beta.maxit = beta.maxit,
      beta.score.tol = beta.score.tol,
      context = context,
      debug = debug,
      seed = seed,
      fit.idx = fit.idx
    )
    if (grepl("^Sample", context) && nL == 3L) {
      fit$component_means <- c(0.15, 0.5, 1.5)
    }
    fit
  }
}

test_that("optional H skips gracefully when H normalization fails", {
  inputs <- make_calibration_inputs(n_probes = 1500L, seed = 11L)
  local_mocked_bindings(
    fit_mixture = mock_h_failure(fit_mixture),
    .package = "betanorm"
  )

  result <- suppressWarnings(bmiq_calibration(
    datM = inputs$datM,
    goldstandard.beta = inputs$gold,
    doH = TRUE,
    h.policy = "optional",
    nfit = 1500L,
    niter = 8L,
    verbose = FALSE
  ))

  expect_true(result$success[1L])
  expect_false(isTRUE(result$h.applied[1L]))
  expect_true(all(is.finite(result$calibrated[1L, ])))
})

test_that("required H fails for the same fit", {
  inputs <- make_calibration_inputs(n_probes = 1500L, seed = 11L)
  local_mocked_bindings(
    fit_mixture = mock_h_failure(fit_mixture),
    .package = "betanorm"
  )

  expect_error(
    suppressWarnings(bmiq_calibration(
      datM = inputs$datM,
      goldstandard.beta = inputs$gold,
      doH = TRUE,
      h.policy = "require",
      on.sample.error = "stop",
      nfit = 1500L,
      niter = 8L,
      verbose = FALSE
    )),
    regexp = "H|failed|methylated"
  )
})
