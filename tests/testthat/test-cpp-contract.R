# Contract tests for the compiled layer.
#
# These pin the parts of the C++ interface that Armadillo used to provide
# implicitly: argument acceptance, deep-copy semantics, and the exact R shape
# of every element of the returned list. R code in bmiq_calibration.R indexes
# em$a[, 1L], em$a[ord, , drop = FALSE] and em$w[, ord, drop = FALSE], so the
# K x 1 / n x K shapes are load-bearing even though nothing else asserts them.
#
# test-ordering-h-policy.R exercises canonicalize_em_components() against a
# hand-built mock list, so it does NOT constrain what the compiled code
# actually returns. That is the gap these tests close.

em_fit <- function(n = 600L, nL = 3L, seed = 101L, ...) {
  y <- simulate_beta_mixture(n, nL = nL, seed = seed)
  thresholds <- if (nL == 2L) 0.5 else c(0.25, 0.75)
  beta_mixture_em_cpp(
    y = y,
    initial_responsibility = hard_responsibility(y, thresholds, nL),
    nL = nL,
    ...
  )
}

test_that("beta_mixture_em_cpp returns a, b, mu as nL x 1 double matrices", {
  for (nL in c(2L, 3L)) {
    fit <- em_fit(nL = nL, maxiter = 6L)

    for (field in c("a", "b", "mu")) {
      value <- fit[[field]]
      expect_true(is.matrix(value), info = paste(field, "nL =", nL))
      expect_identical(dim(value), c(nL, 1L), info = paste(field, "nL =", nL))
      expect_identical(typeof(value), "double", info = paste(field, "nL =", nL))
      expect_null(dimnames(value), info = paste(field, "nL =", nL))
    }
  }
})

test_that("beta_mixture_em_cpp returns eta as a plain numeric vector", {
  # RcppArmadillo wraps arma::vec without a dim attribute. A NumericMatrix
  # replacement would silently turn this into an nL x 1 matrix; downstream
  # `em$eta[ord]` would still work, so nothing else would catch the change.
  for (nL in c(2L, 3L)) {
    eta <- em_fit(nL = nL, maxiter = 6L)$eta
    expect_null(dim(eta), info = paste("nL =", nL))
    expect_identical(typeof(eta), "double", info = paste("nL =", nL))
    expect_length(eta, nL)
    expect_equal(sum(eta), 1, tolerance = 1e-10)
  }
})

test_that("beta_mixture_em_cpp returns w as an n x nL double matrix", {
  n <- 600L
  for (nL in c(2L, 3L)) {
    w <- em_fit(n = n, nL = nL, maxiter = 6L)$w
    expect_true(is.matrix(w), info = paste("nL =", nL))
    expect_identical(dim(w), c(n, nL), info = paste("nL =", nL))
    expect_identical(typeof(w), "double", info = paste("nL =", nL))
    expect_equal(rowSums(w), rep(1, n), tolerance = 1e-10)
  }
})

test_that("beta_mixture_em_cpp returns the documented scalar and status types", {
  fit <- em_fit(nL = 3L, maxiter = 6L)

  expect_type(fit$llike, "double")
  expect_length(fit$llike, 1L)
  expect_type(fit$iterations, "integer")
  expect_length(fit$iterations, 1L)
  expect_type(fit$converged, "logical")
  expect_length(fit$converged, 1L)
  expect_type(fit$parameter_criterion, "double")
  expect_type(fit$loglik_criterion, "double")
  expect_identical(fit$nL, 3L)

  expect_type(fit$fit_status, "character")
  expect_length(fit$fit_status, 3L)
  expect_type(fit$fit_reason, "character")
  expect_length(fit$fit_reason, 3L)

  # Element order matters: bmiq_calibration.R reorders by position elsewhere.
  expect_identical(
    names(fit),
    c(
      "a", "b", "eta", "mu", "w", "llike", "iterations", "converged",
      "parameter_criterion", "loglik_criterion", "fit_status", "fit_reason",
      "nL"
    )
  )
})

test_that("beta_mixture_em_cpp does not mutate initial_responsibility", {
  # arma::mat responsibility = initial_responsibility makes a deep copy. A
  # NumericMatrix assignment is a shallow SEXP copy, so the E-step would write
  # straight into the caller's matrix. This is the single highest-risk
  # difference in an Armadillo -> Rcpp port; keep this test.
  y <- simulate_beta_mixture(400L, nL = 3L, seed = 77L)
  w0 <- hard_responsibility(y, c(0.25, 0.75), 3L)
  w0_before <- w0
  storage.mode(w0_before) <- "double"

  fit <- beta_mixture_em_cpp(
    y = y,
    initial_responsibility = w0,
    nL = 3L,
    maxiter = 5L
  )

  expect_identical(w0, w0_before)
  # And the returned matrix must not alias the input either.
  expect_false(isTRUE(all.equal(unname(fit$w), unname(w0_before))))
})

test_that("beta_mixture_em_cpp validates nL and responsibility dimensions", {
  y <- simulate_beta_mixture(300L, nL = 3L, seed = 5L)
  w0 <- hard_responsibility(y, c(0.25, 0.75), 3L)

  expect_error(
    beta_mixture_em_cpp(y = y, initial_responsibility = w0, nL = 4L),
    regexp = "nL must be 2 or 3"
  )
  expect_error(
    beta_mixture_em_cpp(y = y, initial_responsibility = w0, nL = 1L),
    regexp = "nL must be 2 or 3"
  )
  # nL disagrees with ncol(w0)
  expect_error(
    beta_mixture_em_cpp(y = y, initial_responsibility = w0, nL = 2L),
    regexp = "initial_responsibility"
  )
  # nrow(w0) disagrees with length(y)
  expect_error(
    beta_mixture_em_cpp(
      y = y[1:100],
      initial_responsibility = w0,
      nL = 3L
    ),
    regexp = "initial_responsibility"
  )
})

test_that("beta_mixture_em_cpp accepts an integer responsibility matrix", {
  # Both arma::mat and NumericMatrix coerce INTSXP, but via different paths;
  # pin it so a hand-rolled REAL() port cannot regress into reading garbage.
  y <- simulate_beta_mixture(300L, nL = 2L, seed = 6L)
  w0 <- hard_responsibility(y, 0.5, 2L)
  storage.mode(w0) <- "integer"

  expect_no_error(
    beta_mixture_em_cpp(
      y = y,
      initial_responsibility = w0,
      nL = 2L,
      maxiter = 3L
    )
  )
})

test_that("beta_mixture_em_cpp is deterministic", {
  a <- em_fit(nL = 3L, seed = 31L, maxiter = 8L)
  b <- em_fit(nL = 3L, seed = 31L, maxiter = 8L)
  expect_identical(a, b)
})

test_that("debug traces are attached only when requested", {
  quiet <- em_fit(nL = 3L, maxiter = 4L, debug = FALSE)
  expect_null(quiet$parameter_criterion_trace)
  expect_null(quiet$loglik_criterion_trace)

  loud <- em_fit(nL = 3L, maxiter = 4L, debug = TRUE)
  expect_type(loud$parameter_criterion_trace, "double")
  expect_null(dim(loud$parameter_criterion_trace))
  expect_length(loud$parameter_criterion_trace, loud$iterations)
  expect_length(loud$loglik_criterion_trace, loud$iterations)
})

test_that("scan_finite_unit_interval_cpp accepts the shapes callers pass", {
  # bmiq_calibration() passes datM (a matrix) and goldstandard.beta (a bare
  # vector) directly. The NumericVector signature accepts both; the old
  # arma::mat one forced a matrix(..., ncol = 1L) wrapper around the gold
  # vector, which the port deleted. Keep all three shapes covered so a future
  # narrowing to NumericMatrix cannot pass silently.
  m <- matrix(c(0.1, 0.5, 0.9, 0.2), nrow = 2L)
  expect_no_error(scan_finite_unit_interval_cpp(m, name = "m"))
  expect_no_error(
    scan_finite_unit_interval_cpp(matrix(c(0.1, 0.9), ncol = 1L), name = "v")
  )
  expect_no_error(scan_finite_unit_interval_cpp(c(0.1, 0.9), name = "bare"))
})

test_that("scan_finite_unit_interval_cpp enforces closed and open intervals", {
  boundary <- matrix(c(0, 0.5, 1), ncol = 1L)

  expect_no_error(
    scan_finite_unit_interval_cpp(boundary, name = "x", require_open = FALSE)
  )
  expect_error(
    scan_finite_unit_interval_cpp(boundary, name = "x", require_open = TRUE),
    regexp = "strictly inside"
  )
  expect_error(
    scan_finite_unit_interval_cpp(matrix(1.5, ncol = 1L), name = "x"),
    regexp = "\\[0, 1\\]"
  )
  expect_error(
    scan_finite_unit_interval_cpp(matrix(-1e-12, ncol = 1L), name = "x"),
    regexp = "\\[0, 1\\]"
  )
  for (bad in c(NA_real_, NaN, Inf, -Inf)) {
    expect_error(
      scan_finite_unit_interval_cpp(matrix(bad, ncol = 1L), name = "zz"),
      regexp = "finite"
    )
  }
})

test_that("scan_finite_unit_interval_cpp reports the supplied name", {
  expect_error(
    scan_finite_unit_interval_cpp(matrix(2, ncol = 1L), name = "my_argument"),
    regexp = "my_argument"
  )
})

test_that("gather/scatter round-trip preserves a samples x probes matrix", {
  set.seed(3L)
  x <- matrix(runif(9L * 5L), nrow = 9L, ncol = 5L)

  destination <- matrix(NA_real_, nrow = 9L, ncol = 5L)
  for (start in seq.int(1L, 9L, by = 4L)) {
    count <- min(4L, 9L - start + 1L)
    block <- gather_sample_block_cpp(x, start, count)
    expect_identical(dim(block), c(5L, count))
    expect_identical(block, t(x[start:(start + count - 1L), , drop = FALSE]))
    scatter_sample_block_cpp(destination, block, start)
  }
  expect_identical(destination, x)
})

test_that("gather_sample_block_cpp does not alias its input", {
  x <- matrix(runif(12L), nrow = 4L, ncol = 3L)
  x_before <- x
  block <- gather_sample_block_cpp(x, 1L, 4L)
  block[] <- -1
  expect_identical(x, x_before)
})

test_that("gather/scatter reject out-of-range blocks", {
  x <- matrix(runif(12L), nrow = 4L, ncol = 3L)
  expect_error(gather_sample_block_cpp(x, 0L, 2L), regexp = "Invalid")
  expect_error(gather_sample_block_cpp(x, 4L, 2L), regexp = "Invalid")
  expect_error(gather_sample_block_cpp(x, 1L, 0L), regexp = "Invalid")

  destination <- matrix(0, nrow = 4L, ncol = 3L)
  good <- gather_sample_block_cpp(x, 1L, 2L)
  expect_error(
    scatter_sample_block_cpp(destination, good, 4L),
    regexp = "Invalid"
  )
  expect_error(
    scatter_sample_block_cpp(destination, matrix(0, 2L, 2L), 1L),
    regexp = "Invalid"
  )
})
