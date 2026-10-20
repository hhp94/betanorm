# The exported functions must leave the caller's random state as found.
# draw_fit_indices() restores it around its own draw, but Rcpp's default
# export wrapper opens an RNGScope, whose GetRNGstate() creates a
# time-seeded .Random.seed in a session that has none. No kernel draws a
# random number, so every export is marked rng = false.

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

test_that("bmiq_calibration leaves the random state as found (raw gold)", {
  inputs <- make_calibration_inputs(n_probes = 500L, seed = 3L)

  expect_random_state_kept(function() {
    bmiq_calibration(
      datM = inputs$datM,
      goldstandard.beta = inputs$gold,
      nfit = 500L,
      verbose = FALSE
    )
  })
})

test_that("bmiq_calibration leaves the random state as found (prefit)", {
  inputs <- make_calibration_inputs(n_probes = 500L, seed = 3L)
  gold.fit <- bmiq_gold_fit(inputs$gold, nfit = 500L, verbose = FALSE)

  expect_random_state_kept(function() {
    bmiq_calibration(
      datM = inputs$datM,
      goldstandard.beta = gold.fit,
      nfit = 500L,
      verbose = FALSE
    )
  })
})

test_that("bmiq_gold_fit leaves the random state as found", {
  inputs <- make_calibration_inputs(n_probes = 500L, seed = 3L)

  expect_random_state_kept(function() {
    bmiq_gold_fit(inputs$gold, nfit = 500L, verbose = FALSE)
  })
})

test_that("quantile_norm leaves the random state as found", {
  obj <- matrix(seq_len(24L) / 25, nrow = 4L)
  target <- seq(0, 1, length.out = 6L)

  expect_random_state_kept(function() quantile_norm(obj, target))
})
