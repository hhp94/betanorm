## Builds FIXED_QTABLE, the gold quantile table bmiq_calibration_fixed()
## reads, into R/sysdata.rda, which holds nothing else. Run from the package
## root after data-raw/bmiq_gold_fixed.R. Rationale: dev/DECISIONS.md.
devtools::load_all(".", quiet = TRUE)

QTABLE_N <- 8192L
QTABLE_P_MIN <- 1e-30

# Node i is u1 - (n - i) * h; fixed_eval_half() rebuilds it the same way.
qtable_grid <- function() {
  n <- QTABLE_N
  u1 <- log(0.5)
  h <- (u1 - log(QTABLE_P_MIN)) / (n - 1)
  list(n = n, u0 = u1 - (n - 1) * h, u1 = u1, h = h)
}

make_half <- function(a, b) {
  grid <- qtable_grid()
  n <- grid[["n"]]
  h <- grid[["h"]]
  u <- grid[["u1"]] - ((n - 1):0) * h
  x <- stats::qbeta(u, a, b, log.p = TRUE)
  v <- log(x)
  # dv/du = p / (x f(x)), and its derivative in u, both analytic.
  m <- exp(u - v - stats::dbeta(x, a, b, log = TRUE))
  k <- m * (1 - m * (a - (b - 1) * x / (1 - x)))

  # Per interval, the quintic in t whose value, first and second derivative
  # match both ends (derivatives with respect to t are h m and h^2 k).
  i0 <- seq_len(n - 1L)
  i1 <- i0 + 1L
  v0 <- v[i0]
  v1 <- v[i1]
  d0 <- h * m[i0]
  d1 <- h * m[i1]
  e0 <- h^2 * k[i0]
  e1 <- h^2 * k[i1]
  A <- v1 - v0 - d0 - e0 / 2
  B <- d1 - d0 - e0
  C <- e1 - e0
  list(
    a = a,
    b = b,
    n = n,
    u0 = grid[["u0"]],
    u1 = grid[["u1"]],
    h = h,
    c0 = v0,
    c1 = d0,
    c2 = e0 / 2,
    c3 = 10 * A - 4 * B + C / 2,
    c4 = -15 * A + 7 * B - C,
    c5 = 6 * A - 3 * B + C / 2
  )
}

make_component <- function(a, b) {
  list(a = a, b = b, ab = make_half(a, b), ba = make_half(b, a))
}

# The U and M components of each gold set; H is never quantile-mapped.
FIXED_QTABLE <- lapply(FIXED_GOLD, function(shapes) {
  list(
    U = make_component(shapes[["a"]][1L], shapes[["b"]][1L]),
    M = make_component(shapes[["a"]][3L], shapes[["b"]][3L])
  )
})

save(FIXED_QTABLE, file = file.path("R", "sysdata.rda"), compress = "xz")

# Round trip: the saved table must load back to the same bits.
loaded <- new.env()
load(file.path("R", "sysdata.rda"), envir = loaded)
stopifnot(
  identical(ls(loaded), "FIXED_QTABLE"),
  identical(loaded[["FIXED_QTABLE"]], FIXED_QTABLE)
)
