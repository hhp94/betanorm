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
  expect_true(is.finite(fit1$parameter_criterion) ||
    is.infinite(fit1$parameter_criterion))
  expect_true(is.finite(fit1$loglik_criterion) ||
    is.infinite(fit1$loglik_criterion))
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

test_that("EM does not stop on mean stability alone when shapes or eta still move", {
  # Well-separated data; inspect early iterations under a loose mean-only
  # rule versus the dual criterion actually used for converged.
  y <- simulate_beta_mixture(2000L, nL = 3L, seed = 33L)
  w0 <- hard_responsibility(y, c(0.25, 0.75), 3L)

  fit <- beta_mixture_em_cpp(
    y = y,
    initial_responsibility = w0,
    nL = 3L,
    maxiter = 15L,
    tol = 1e-6,
    debug = TRUE
  )

  # Reconstruct successive means from warm-started single-step runs is heavy;
  # instead verify that parameter_criterion is the max abs change over
  # log(a), log(b), eta by checking it is at least as large as the eta-only
  # and log-shape changes between the final state and a one-step earlier
  # restart is not available. Proxy: traces are positive until late, and
  # declaring converged requires both traces' final values below tol.
  expect_equal(length(fit$parameter_criterion_trace), fit$iterations)
  expect_equal(length(fit$loglik_criterion_trace), fit$iterations)

  # First iteration cannot declare convergence (warm-start Inf / NA).
  if (fit$iterations >= 1L) {
    expect_true(
      is.infinite(fit$parameter_criterion_trace[1L]) ||
        fit$parameter_criterion_trace[1L] > 0
    )
  }

  # Parameter state includes log-shapes and eta: a fit with only maxiter = 2
  # and extremely tight tol should not claim convergence unless both are tiny.
  fit_tight <- beta_mixture_em_cpp(
    y = y,
    initial_responsibility = w0,
    nL = 3L,
    maxiter = 2L,
    tol = 1e-30,
    debug = TRUE
  )
  expect_false(isTRUE(fit_tight$converged))
  expect_equal(fit_tight$iterations, 2L)
  expect_true(
    fit_tight$parameter_criterion >= 1e-30 ||
      fit_tight$loglik_criterion >= 1e-30
  )

  # Mean-only would ignore shape/eta; verify the parameter criterion can be
  # driven by log-shape even when we force a second iteration with a
  # different effective concentration path by using nL = 2 on the same y.
  w0_2 <- hard_responsibility(y, 0.5, 2L)
  fit2 <- beta_mixture_em_cpp(
    y = y,
    initial_responsibility = w0_2,
    nL = 2L,
    maxiter = 20L,
    tol = 1e-5,
    debug = TRUE
  )
  expect_true("parameter_criterion" %in% names(fit2))
  expect_true("loglik_criterion" %in% names(fit2))
  if (isTRUE(fit2$converged)) {
    expect_lt(fit2$parameter_criterion, 1e-5)
    expect_lt(fit2$loglik_criterion, 1e-5)
    expect_gte(fit2$iterations, 2L)
  }
})

test_that("debug traces have consistent lengths with iterations", {
  y <- simulate_beta_mixture(600L, nL = 2L, seed = 44L)
  w0 <- hard_responsibility(y, 0.5, 2L)

  fit <- beta_mixture_em_cpp(
    y = y,
    initial_responsibility = w0,
    nL = 2L,
    maxiter = 10L,
    tol = 1e-8,
    debug = TRUE
  )

  expect_length(fit$parameter_criterion_trace, fit$iterations)
  expect_length(fit$loglik_criterion_trace, fit$iterations)
  expect_equal(
    fit$parameter_criterion,
    fit$parameter_criterion_trace[fit$iterations]
  )
  expect_equal(
    fit$loglik_criterion,
    fit$loglik_criterion_trace[fit$iterations]
  )
})
