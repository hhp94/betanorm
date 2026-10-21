## Builds the gold quantile tables bmiq_calibration_fixed() reads: one file
## per gold set, data-raw/fixed_qtable_<gold>.rds, then FIXED_QTABLE in
## R/sysdata.rda (which holds nothing else), the named list of those files.
## Run from the package root after data-raw/bmiq_gold_fixed.R. A file is
## written only when it is missing or holds another value, so adding a gold
## set leaves the files of the others byte for byte. Rationale:
## dev/DECISIONS.md.
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
tables <- lapply(FIXED_GOLD, function(shapes) {
  list(
    U = make_component(shapes[["a"]][1L], shapes[["b"]][1L]),
    M = make_component(shapes[["a"]][3L], shapes[["b"]][3L])
  )
})

qtable_path <- function(gold) {
  file.path("data-raw", paste0("fixed_qtable_", gold, ".rds"))
}

for (gold in names(tables)) {
  path <- qtable_path(gold)
  if (!file.exists(path) || !identical(readRDS(path), tables[[gold]])) {
    saveRDS(tables[[gold]], path, compress = "xz")
  }
}

# FIXED_QTABLE is read back from the files, so it holds what they hold.
FIXED_QTABLE <- lapply(
  stats::setNames(nm = names(tables)),
  function(gold) readRDS(qtable_path(gold))
)
stopifnot(identical(FIXED_QTABLE, tables))

# The same rule for R/sysdata.rda: rewritten only when its table changes.
sysdata <- file.path("R", "sysdata.rda")
stored <- new.env()
if (file.exists(sysdata)) {
  load(sysdata, envir = stored)
}
if (
  !identical(ls(stored), "FIXED_QTABLE") ||
    !identical(stored[["FIXED_QTABLE"]], FIXED_QTABLE)
) {
  save(FIXED_QTABLE, file = sysdata, compress = "xz")
}

# Round trip: the saved table must load back to the same bits.
loaded <- new.env()
load(sysdata, envir = loaded)
stopifnot(
  identical(ls(loaded), "FIXED_QTABLE"),
  identical(loaded[["FIXED_QTABLE"]], FIXED_QTABLE)
)
