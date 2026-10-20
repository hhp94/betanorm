test_that("EM always returns dual criteria and never converges after one iteration", {
  y <- simulate_beta_mixture(800L, nL = 3L, seed = 11L)
  w0 <- hard_responsibility(y, c(0.25, 0.75), 3L)

  fit1 <- beta_mixture_em_cpp(
    y = y,
    initial_responsibility = w0,
    nL = 3L,
    maxiter = 1L,
    tol = 0.1,
    debug = FALSE
  )

  expect_false(isTRUE(fit1$converged))
  expect_equal(fit1$iterations, 1L)
  expect_true(
    is.finite(fit1$parameter_criterion) ||
      is.infinite(fit1$parameter_criterion)
  )
  expect_true(
    is.finite(fit1$loglik_criterion) ||
      is.infinite(fit1$loglik_criterion)
  )
  expect_null(fit1$parameter_criterion_trace)
  expect_null(fit1$loglik_criterion_trace)
})

test_that("EM convergence requires both parameter and relative loglik criteria", {
  y <- simulate_beta_mixture(1200L, nL = 3L, seed = 22L)
  w0 <- hard_responsibility(y, c(0.25, 0.75), 3L)

  fit <- beta_mixture_em_cpp(
    y = y,
    initial_responsibility = w0,
    nL = 3L,
    maxiter = 25L,
    tol = 1e-4,
    debug = TRUE
  )

  expect_true(is.numeric(fit$parameter_criterion))
  expect_true(is.numeric(fit$loglik_criterion))
  expect_equal(length(fit$parameter_criterion_trace), fit$iterations)
  expect_equal(length(fit$loglik_criterion_trace), fit$iterations)

  if (isTRUE(fit$converged)) {
    expect_gte(fit$iterations, 2L)
    expect_lt(fit$parameter_criterion, 1e-4)
    expect_lt(fit$loglik_criterion, 1e-4)
    expect_true(all(fit$fit_status == "converged"))
  } else {
    # Exhausted iterations without dual-criterion success.
    expect_equal(fit$iterations, 25L)
    expect_true(
      fit$parameter_criterion >= 1e-4 || fit$loglik_criterion >= 1e-4
    )
  }
})

test_that("a component that keeps its previous shapes keeps their status", {
  # The nL = 2 snapshot case: at its last iteration the ascent guard keeps
  # both components' previous shapes. Rerunning one iteration short gives
  # the fit those shapes came from.
  y <- simulate_beta_mixture(2000L, nL = 2L, seed = 4321L)
  w0 <- hard_responsibility(y, 0.5, 2L)
  run <- function(maxiter) {
    beta_mixture_em_cpp(
      y = y,
      initial_responsibility = w0,
      nL = 2L,
      maxiter = maxiter,
      tol = 1e-8
    )
  }

  fit <- run(40L)
  previous <- run(fit$iterations - 1L)

  expect_true(all(startsWith(fit$fit_reason, "retained previous shapes")))
  expect_identical(fit$a, previous$a)
  expect_identical(fit$b, previous$b)
  # The status is that of the fit that produced the shapes, not a blanket
  # "converged".
  expect_identical(fit$fit_status, previous$fit_status)
})
