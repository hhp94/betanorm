test_that("ordering helpers enforce strict vs anchor policies", {
  expect_error(
    requireOrderedComponentMeans(c(0.1, 0.9, 0.5), "disordered"),
    regexp = "not strictly increasing"
  )
  # U=0.1 < M=0.5 is enough for anchors even if H is disordered.
  expect_silent(
    requireOrderedAnchors(c(0.1, 0.9, 0.5), "disordered-H")
  )
  expect_error(
    requireOrderedAnchors(c(0.8, 0.5, 0.4), "collapsed"),
    regexp = "not separated"
  )
  expect_silent(
    requireOrderedComponentMeans(c(0.1, 0.4, 0.8), "ordered")
  )
})

test_that("gold-standard fits still require complete component ordering", {
  inputs <- make_calibration_inputs(n_probes = 800L, seed = 10L)
  real_pipeline <- fitMixturePipeline

  local_mocked_bindings(
    fitMixturePipeline = function(beta,
                                  thresholds,
                                  nL,
                                  nfit,
                                  niter,
                                  tol,
                                  beta.maxit,
                                  beta.score.tol,
                                  fit.policy,
                                  context,
                                  debug = FALSE,
                                  seed = 1L,
                                  mean.order = c("strict", "anchors")) {
      mean.order <- match.arg(mean.order)
      fit <- real_pipeline(
        beta = beta,
        thresholds = thresholds,
        nL = nL,
        nfit = nfit,
        niter = niter,
        tol = tol,
        beta.maxit = beta.maxit,
        beta.score.tol = beta.score.tol,
        fit.policy = fit.policy,
        context = context,
        debug = debug,
        seed = seed,
        mean.order = mean.order
      )
      if (identical(context, "Gold-standard")) {
        fit$component_means <- c(0.2, 0.9, 0.5)
        requireOrderedComponentMeans(
          fit$component_means,
          "Gold-standard mixture"
        )
      }
      fit
    },
    .package = "bmiqpp"
  )

  expect_error(
    bmiq_calibration(
      datM = inputs$datM,
      goldstandard.beta = inputs$gold,
      nfit = 800L,
      verbose = FALSE
    ),
    regexp = "not strictly increasing|Gold-standard"
  )
})

test_that("optional H accepts separated U/M anchors despite disordered H", {
  inputs <- make_calibration_inputs(n_probes = 1500L, seed = 11L)
  real_pipeline <- fitMixturePipeline

  local_mocked_bindings(
    fitMixturePipeline = function(beta,
                                  thresholds,
                                  nL,
                                  nfit,
                                  niter,
                                  tol,
                                  beta.maxit,
                                  beta.score.tol,
                                  fit.policy,
                                  context,
                                  debug = FALSE,
                                  seed = 1L,
                                  mean.order = c("strict", "anchors")) {
      mean.order <- match.arg(mean.order)
      fit <- real_pipeline(
        beta = beta,
        thresholds = thresholds,
        nL = nL,
        nfit = nfit,
        niter = niter,
        tol = tol,
        beta.maxit = beta.maxit,
        beta.score.tol = beta.score.tol,
        fit.policy = fit.policy,
        context = context,
        debug = debug,
        seed = seed,
        mean.order = mean.order
      )
      if (grepl("^Sample", context) && mean.order == "anchors" && nL == 3L) {
        fit$component_means <- c(0.15, 0.85, 0.70)
      }
      fit
    },
    .package = "bmiqpp"
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

test_that("required H fails for the same disordered fit", {
  inputs <- make_calibration_inputs(n_probes = 1500L, seed = 11L)
  real_pipeline <- fitMixturePipeline

  local_mocked_bindings(
    fitMixturePipeline = function(beta,
                                  thresholds,
                                  nL,
                                  nfit,
                                  niter,
                                  tol,
                                  beta.maxit,
                                  beta.score.tol,
                                  fit.policy,
                                  context,
                                  debug = FALSE,
                                  seed = 1L,
                                  mean.order = c("strict", "anchors")) {
      mean.order <- match.arg(mean.order)
      fit <- real_pipeline(
        beta = beta,
        thresholds = thresholds,
        nL = nL,
        nfit = nfit,
        niter = niter,
        tol = tol,
        beta.maxit = beta.maxit,
        beta.score.tol = beta.score.tol,
        fit.policy = fit.policy,
        context = context,
        debug = debug,
        seed = seed,
        mean.order = mean.order
      )
      if (grepl("^Sample", context) && mean.order == "anchors" && nL == 3L) {
        fit$component_means <- c(0.15, 0.85, 0.70)
      }
      fit
    },
    .package = "bmiqpp"
  )

  expect_error(
    bmiq_calibration(
      datM = inputs$datM,
      goldstandard.beta = inputs$gold,
      doH = TRUE,
      h.policy = "require",
      on.sample.error = "stop",
      nfit = 1500L,
      niter = 8L,
      verbose = FALSE
    ),
    regexp = "not strictly increasing|H|failed"
  )
})