# Helpers shared by the BMIQ test files.

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

# Run f() from a session with no .Random.seed (none may be left behind),
# then from a set seed (it must come back identical). The caller's own
# state is restored on exit.
expect_random_state_kept <- function(f) {
  global <- globalenv()
  had.seed <- exists(".Random.seed", envir = global, inherits = FALSE)
  saved.seed <- if (had.seed) get(".Random.seed", envir = global)
  on.exit(
    if (had.seed) {
      assign(".Random.seed", saved.seed, envir = global)
    } else if (exists(".Random.seed", envir = global, inherits = FALSE)) {
      rm(".Random.seed", envir = global)
    },
    add = TRUE
  )

  if (exists(".Random.seed", envir = global, inherits = FALSE)) {
    rm(".Random.seed", envir = global)
  }
  f()
  expect_false(exists(".Random.seed", envir = global, inherits = FALSE))

  set.seed(20261020L)
  before <- get(".Random.seed", envir = global)
  f()
  expect_identical(
    get0(".Random.seed", envir = global, inherits = FALSE),
    before
  )
}
