# Shared synthetic Beta-mixture data for BMIQ tests.

clip01 <- function(x, eps = 1e-4) {
  pmin(1 - eps, pmax(eps, x))
}

#' Draw a simple U/M or U/H/M beta mixture on (0, 1).
simulate_beta_mixture <- function(
  n,
  nL = 3L,
  seed = 1L
) {
  set.seed(seed)
  if (nL == 2L) {
    n_u <- floor(n / 2)
    n_m <- n - n_u
    y <- c(rbeta(n_u, 2, 8), rbeta(n_m, 8, 2))
  } else {
    n_u <- floor(0.4 * n)
    n_m <- floor(0.4 * n)
    n_h <- n - n_u - n_m
    y <- c(
      rbeta(n_u, 2, 10),
      rbeta(n_h, 3, 3),
      rbeta(n_m, 10, 2)
    )
  }
  clip01(y)
}

#' Hard one-hot responsibilities from ordered thresholds.
hard_responsibility <- function(y, thresholds, nL) {
  class <- rep.int(1L, length(y))
  for (boundary in seq_along(thresholds)) {
    class[y > thresholds[boundary]] <- boundary + 1L
  }
  w <- matrix(0, nrow = length(y), ncol = nL)
  w[cbind(seq_along(y), class)] <- 1
  w
}

#' Build a one-row datM matrix and matching gold vector for BMIQcalibration.
make_calibration_inputs <- function(
  n_probes = 1500L,
  nL = 3L,
  seed = 42L
) {
  gold <- simulate_beta_mixture(n_probes, nL = nL, seed = seed)
  # Mild sample shift so calibration has something to do.
  sample <- simulate_beta_mixture(n_probes, nL = nL, seed = seed + 7L)
  # Slight rightward shift of the sample U mode relative to gold.
  sample <- clip01(sample + 0.02)

  datM <- matrix(sample, nrow = 1L, dimnames = list("s1", NULL))
  list(datM = datM, gold = gold)
}
