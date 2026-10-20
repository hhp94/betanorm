# The exported functions must leave the caller's random state as found.
# draw_fit_indices() restores it around its own draw, but Rcpp's default
# export wrapper opens an RNGScope, whose GetRNGstate() creates a
# time-seeded .Random.seed in a session that has none. No kernel draws a
# random number, so every export is marked rng = false.

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

test_that("the fit subset does not follow the caller's RNG kind", {
  global <- globalenv()
  saved.kind <- RNGkind()
  saved.seed <- get0(".Random.seed", envir = global, inherits = FALSE)
  on.exit(
    {
      suppressWarnings(RNGkind(
        kind = saved.kind[[1L]],
        normal.kind = saved.kind[[2L]],
        sample.kind = saved.kind[[3L]]
      ))
      if (is.null(saved.seed)) {
        rm(".Random.seed", envir = global)
      } else {
        assign(".Random.seed", saved.seed, envir = global)
      }
    },
    add = TRUE
  )

  inputs <- make_calibration_inputs(n_probes = 3000L, seed = 5L)
  run <- function() {
    suppressWarnings(bmiq_calibration(
      datM = inputs$datM,
      goldstandard.beta = inputs$gold,
      nfit = 2000L,
      verbose = FALSE
    ))$calibrated
  }

  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  reference <- run()

  RNGkind("L'Ecuyer-CMRG")
  expect_identical(run(), reference)
  expect_identical(RNGkind()[[1L]], "L'Ecuyer-CMRG")

  RNGkind("Mersenne-Twister")
  suppressWarnings(RNGkind(sample.kind = "Rounding"))
  expect_identical(run(), reference)
  expect_identical(RNGkind()[[3L]], "Rounding")
})
