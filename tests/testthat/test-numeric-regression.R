# Numeric regression net for the Armadillo -> Rcpp port.
#
# No arithmetic changes in that port: the Newton step solves its 2x2 system by
# hand (explicit determinant, no arma::solve), and every special function comes
# from R's Rmath (R::digamma, R::trigamma) or libstdc++ (std::lgamma). Nothing
# routes through BLAS/LAPACK. The port should therefore be bit-for-bit
# identical, and these snapshots hold it to that.
#
# CAVEAT: no baseline snapshot was captured before the port, so the values in
# _snaps/ were generated on the Rcpp build and do NOT prove the port preserved
# numerics. To establish that retroactively, check out the last RcppArmadillo
# commit, run this file to write the snapshots, then return to HEAD and rerun.
# From here on they serve their normal purpose: guarding future changes.

test_that("EM fit is numerically stable for nL = 3", {
  y <- simulate_beta_mixture(2000L, nL = 3L, seed = 1234L)
  w0 <- hard_responsibility(y, c(0.25, 0.75), 3L)

  fit <- beta_mixture_em_cpp(
    y = y,
    initial_responsibility = w0,
    nL = 3L,
    maxiter = 25L,
    tol = 1e-6
  )

  expect_snapshot_value(
    list(
      a = as.numeric(fit$a[, 1L]),
      b = as.numeric(fit$b[, 1L]),
      eta = as.numeric(fit$eta),
      mu = as.numeric(fit$mu[, 1L]),
      llike = fit$llike,
      iterations = fit$iterations,
      converged = fit$converged,
      parameter_criterion = fit$parameter_criterion,
      loglik_criterion = fit$loglik_criterion,
      fit_status = as.character(fit$fit_status),
      # Full responsibility matrix is too large to snapshot; column sums are
      # sensitive to every element of it.
      w_colsums = as.numeric(colSums(fit$w))
    ),
    style = "serialize",
    tolerance = 1e-12
  )
})

test_that("EM fit is numerically stable for nL = 2", {
  y <- simulate_beta_mixture(2000L, nL = 2L, seed = 4321L)
  w0 <- hard_responsibility(y, 0.5, 2L)

  fit <- beta_mixture_em_cpp(
    y = y,
    initial_responsibility = w0,
    nL = 2L,
    maxiter = 40L,
    tol = 1e-8
  )

  expect_snapshot_value(
    list(
      a = as.numeric(fit$a[, 1L]),
      b = as.numeric(fit$b[, 1L]),
      eta = as.numeric(fit$eta),
      mu = as.numeric(fit$mu[, 1L]),
      llike = fit$llike,
      iterations = fit$iterations,
      converged = fit$converged,
      fit_status = as.character(fit$fit_status),
      w_colsums = as.numeric(colSums(fit$w))
    ),
    style = "serialize",
    tolerance = 1e-12
  )
})

test_that("EM log-likelihood is monotonically non-decreasing", {
  # Property-based backstop that does not depend on a stored snapshot: the
  # generalized-EM ascent guard in the M-step must keep the objective from
  # dropping. Cheap insurance against a port that reorders accumulation.
  for (seed in c(11L, 12L, 13L)) {
    y <- simulate_beta_mixture(900L, nL = 3L, seed = seed)
    w0 <- hard_responsibility(y, c(0.25, 0.75), 3L)

    likelihoods <- vapply(
      seq_len(12L),
      function(maxiter) {
        beta_mixture_em_cpp(
          y = y,
          initial_responsibility = w0,
          nL = 3L,
          maxiter = maxiter,
          tol = 0
        )$llike
      },
      numeric(1L)
    )

    expect_true(
      all(diff(likelihoods) >= -1e-8),
      info = paste(
        "seed",
        seed,
        ":",
        paste(signif(likelihoods, 10), collapse = " ")
      )
    )
  }
})

test_that("bmiq_calibration output is numerically stable end to end", {
  inputs <- make_calibration_inputs(n_probes = 1500L, nL = 3L, seed = 42L)

  result <- bmiq_calibration(
    datM = inputs$datM,
    goldstandard.beta = inputs$gold,
    nL = 3L,
    nfit = 1500L,
    niter = 5L,
    verbose = FALSE
  )

  expect_snapshot_value(
    list(
      calibrated = unname(as.matrix(result)),
      success = result$success,
      h_applied = result$h.applied
    ),
    style = "serialize",
    tolerance = 1e-12
  )
})

test_that("bmiq_calibration output is numerically stable for nL = 2", {
  inputs <- make_calibration_inputs(n_probes = 1500L, nL = 2L, seed = 43L)

  result <- bmiq_calibration(
    datM = inputs$datM,
    goldstandard.beta = inputs$gold,
    nL = 2L,
    nfit = 1500L,
    niter = 10L,
    verbose = FALSE
  )

  expect_snapshot_value(
    list(
      calibrated = unname(as.matrix(result)),
      success = result$success
    ),
    style = "serialize",
    tolerance = 1e-12
  )
})

test_that("bmiq_calibration handles multiple samples across gather/scatter blocks", {
  # sample.block.size is 8, so 11 samples exercises a full block plus a
  # partial one. This is the only path that touches the block gather/scatter
  # code with a non-trivial second block.
  n_probes <- 1200L
  n_samples <- 11L
  gold <- simulate_beta_mixture(n_probes, nL = 3L, seed = 500L)
  datM <- t(vapply(
    seq_len(n_samples),
    function(i) {
      clip01(simulate_beta_mixture(n_probes, nL = 3L, seed = 500L + i))
    },
    numeric(n_probes)
  ))
  rownames(datM) <- paste0("s", seq_len(n_samples))

  result <- bmiq_calibration(
    datM = datM,
    goldstandard.beta = gold,
    nL = 3L,
    nfit = n_probes,
    verbose = FALSE
  )

  calibrated <- as.matrix(result)
  expect_identical(dim(calibrated), c(n_samples, n_probes))
  expect_identical(rownames(calibrated), rownames(datM))
  expect_true(all(result$success))
  expect_true(all(is.finite(calibrated)))
  expect_true(all(calibrated >= 0 & calibrated <= 1))

  # Samples must be independent: calibrating one sample alone must reproduce
  # its row from the multi-sample run exactly. A gather/scatter offset bug
  # shows up here and nowhere else.
  for (i in c(1L, 8L, 9L, 11L)) {
    solo <- bmiq_calibration(
      datM = datM[i, , drop = FALSE],
      goldstandard.beta = gold,
      nL = 3L,
      nfit = n_probes,
      verbose = FALSE
    )
    expect_equal(
      as.numeric(as.matrix(solo)[1L, ]),
      as.numeric(calibrated[i, ]),
      tolerance = 0,
      info = paste("sample", i)
    )
  }
})

test_that("quantile_norm output is numerically stable", {
  # Already pure Rcpp, but Makevars changes affect its compilation too.
  set.seed(7L)
  obj <- matrix(rnorm(6L * 40L), nrow = 6L, ncol = 40L)
  target <- sort(rnorm(40L))

  expect_snapshot_value(
    unname(quantile_norm(obj, target)),
    style = "serialize",
    tolerance = 1e-12
  )

  # Unequal-length target path.
  expect_snapshot_value(
    unname(quantile_norm(obj, sort(rnorm(17L)))),
    style = "serialize",
    tolerance = 1e-12
  )
})
