## Phase 1 of dev/bmiq_calibration_fixed.md: interpolant and grid experiment.
## Gold-side qbeta of bmiq_calibration_fixed() replaced by a table; which
## interpolant (cubic + Newton vs quintic alone) and which grid size n.
## Reference everywhere: the production path qbeta(pbeta(x)) in probability
## space. Run: Rscript dev/qtable-exp.R (package loaded with load_all).

suppressPackageStartupMessages({ devtools::load_all(".", quiet = TRUE); library(bench) })
options(warn = 1, digits = 4, width = 160)
cat("betanorm", as.character(packageVersion("betanorm")), " R", R.version.string, "\n")
EPS <- .Machine$double.eps

class_by_thresh <- function(beta, thresholds) {
  class <- rep.int(1L, length(beta))
  for (k in seq_along(thresholds)) class[beta > thresholds[k]] <- k + 1L
  class
}

## ------------------------------------------------------------- real inputs
## Every gold-side qbeta input of the six EPIC samples at the three call
## sites, reconstructed from the general function's diagnostics and checked
## bit for bit against bmiq_calibration() (the Rmath path; from c40eb20
## bmiq_calibration_fixed() itself reads the table measured here).
data(horvath_goldstandard, package = "betanorm")
data(GPL21145_sample, package = "betanorm")
measured <- colSums(GPL21145_sample == 0) < nrow(GPL21145_sample)
datM <- GPL21145_sample[, measured]
gold <- betanorm:::FIXED_GOLD[["horvath"]]
gold3 <- bmiq_gold_fit(horvath_goldstandard$goldstandard2, nL = 3L, verbose = FALSE)
stopifnot(identical(gold3$a, gold$a), identical(gold3$b, gold$b))
res3 <- bmiq_calibration(datM, gold3, nL = 3L, debug = TRUE, verbose = FALSE)
general <- res3$calibrated

sites <- list()
for (i in seq_len(nrow(datM))) {
  d <- res3$diagnostics$samples[[i]]
  x <- unname(datM[i, ])
  cls <- class_by_thresh(x, d$thresholds)
  mu <- d$component_means
  as <- d$component_a
  bs <- d$component_b
  selU <- which(cls == 1L)
  selUL <- selU[x[selU] <= mu[1]]
  selUR <- selU[x[selU] > mu[1]]
  selM <- which(cls == 3L)
  selMR <- selM[x[selM] >= mu[3]]
  add <- function(name, idx, k, lower) {
    p <- pbeta(x[idx], as[k], bs[k], lower.tail = lower)
    ref <- qbeta(p, gold$a[k], gold$b[k], lower.tail = lower)
    stopifnot(identical(ref, unname(general[i, idx])))
    sites[[length(sites) + 1L]] <<- list(
      site = name, sample = i, comp = if (k == 1L) "U" else "M", lower = lower,
      x = x[idx], a.s = as[k], b.s = bs[k], p = p, u = log(p), ref = ref
    )
  }
  add("UL", selUL, 1L, TRUE)
  add("UR", selUR, 1L, FALSE)
  add("MR", selMR, 3L, FALSE)
}
n.real <- sum(vapply(sites, function(s) length(s$u), 1L))
cat("real inputs:", n.real, "gold-side evaluations over", length(sites),
    "call sites; production path reproduced bit for bit\n")
cat("p = 0 or 1 among them:", sum(vapply(sites, function(s) sum(s$u == -Inf | s$u >= 0), 1L)), "\n\n")

## ------------------------------------------------------------- the tables
P.MIN <- 1e-30

## One half: lower-tail quantile of Beta(a, b), p in [P.MIN, 1/2], grid
## uniform in u = log p. Stores per-interval Horner coefficients for both
## interpolants (v = log x as a polynomial in t = (u - u_i) / h).
## polish: Newton steps on log F(x) = u applied to the nodes at build time
## (Rmath's qbeta with log.p = TRUE is only accurate to ~1e-12).
## Grid anchored at u1 = log(1/2): node i is u1 - (n - i) * h, computed with
## the same arithmetic in eval_half(). (seq() nodes reconstructed from
## u0 + i * h drift by up to 2.5e-11 at n = 4096, a 2.6e-12 floor in x.)
make_half <- function(a, b, n, polish = 0L) {
  u1 <- log(0.5)
  h <- (u1 - log(P.MIN)) / (n - 1)
  u <- u1 - ((n - 1):0) * h
  x <- qbeta(u, a, b, log.p = TRUE)
  lb <- lbeta(a, b)
  for (z in seq_len(polish)) {
    lF <- pbeta(x, a, b, log.p = TRUE)
    lf <- (a - 1) * log(x) + (b - 1) * log1p(-x) - lb
    x <- x - (lF - u) * exp(lF - lf)
  }
  v <- log(x)
  m <- exp(u - v - dbeta(x, a, b, log = TRUE))            # dv/du
  k <- m * (1 - m * (a - (b - 1) * x / (1 - x)))          # d2v/du2
  i0 <- seq_len(n - 1L)
  i1 <- i0 + 1L
  v0 <- v[i0]; v1 <- v[i1]
  d0 <- h * m[i0]; d1 <- h * m[i1]
  e0 <- h^2 * k[i0]; e1 <- h^2 * k[i1]
  cubic <- list(
    c0 = v0, c1 = d0,
    c2 = 3 * (v1 - v0) - 2 * d0 - d1,
    c3 = 2 * (v0 - v1) + d0 + d1
  )
  A <- v1 - v0 - d0 - e0 / 2
  B <- d1 - d0 - e0
  C <- e1 - e0
  quintic <- list(
    c0 = v0, c1 = d0, c2 = e0 / 2,
    c3 = 10 * A - 4 * B + C / 2,
    c4 = -15 * A + 7 * B - C,
    c5 = 6 * A - 3 * B + C / 2
  )
  list(a = a, b = b, lbeta = lbeta(a, b), u0 = u[1L], u1 = u[n], h = h, n = n,
       cubic = cubic, quintic = quintic)
}
make_qtable <- function(a, b, n, polish = 0L) {
  list(a = a, b = b, ab = make_half(a, b, n, polish), ba = make_half(b, a, n, polish))
}

## u inside (u0, u1]: Horner, exp, optional one Newton step on log F(x) = u.
eval_half <- function(u, hf, variant, newton) {
  j <- pmin(floor((hf$u1 - u) / hf$h), hf$n - 2)   # intervals below u1
  i <- hf$n - 1 - j                                 # left node, 1-based
  t <- (u - (hf$u1 - (hf$n - i) * hf$h)) / hf$h     # node as built, bit for bit
  cf <- hf[[variant]]
  v <- if (variant == "cubic") {
    cf$c0[i] + t * (cf$c1[i] + t * (cf$c2[i] + t * cf$c3[i]))
  } else {
    cf$c0[i] + t * (cf$c1[i] + t * (cf$c2[i] + t * (cf$c3[i] + t * (cf$c4[i] + t * cf$c5[i]))))
  }
  x <- exp(v)
  if (newton) {
    lF <- pbeta(x, hf$a, hf$b, log.p = TRUE)
    lf <- (hf$a - 1) * log(x) + (hf$b - 1) * log1p(-x) - hf$lbeta
    x <- x - (lF - u) * exp(lF - lf)
  }
  x
}

## qbeta(p, a, b, lower.tail) with p given as u = log p (valid inputs only).
qbeta_table <- function(u, tab, lower.tail, variant, newton) {
  out <- numeric(length(u))
  zero <- u == -Inf
  one <- u >= 0
  small <- u > tab$ab$u0 & u <= tab$ab$u1                  # p in (P.MIN, 1/2]
  ubig <- log(-expm1(u))                                   # log(1 - p)
  big <- u > tab$ab$u1 & u < 0 & ubig > tab$ab$u0          # p in (1/2, 1 - P.MIN)
  fb <- !(zero | one | small | big)                        # p <= P.MIN: Rmath
  if (any(fb)) out[fb] <- qbeta(u[fb], tab$a, tab$b, lower.tail = lower.tail, log.p = TRUE)
  if (lower.tail) {
    out[small] <- eval_half(u[small], tab$ab, variant, newton)
    out[big] <- 1 - eval_half(ubig[big], tab$ba, variant, newton)
    out[zero] <- 0
    out[one] <- 1
  } else {
    out[small] <- 1 - eval_half(u[small], tab$ba, variant, newton)
    out[big] <- eval_half(ubig[big], tab$ab, variant, newton)
    out[zero] <- 1
    out[one] <- 0
  }
  out
}

## ------------------------------------------------------------- the sweep
## Dense p per (component, tail): grid nodes, interval interiors, random
## interiors, their mirrors above 1/2, and the boundary cases. p is built
## first and u = log(p) taken from it, as in production.
make_sweep_p <- function(hf) {
  nodes <- seq(hf$u0, hf$u1, length.out = hf$n)
  interior <- seq(hf$u0 + hf$h / 2, hf$u1 - hf$h / 2, length.out = 20000L)
  set.seed(1L)
  random <- runif(20000L, hf$u0, hf$u1)
  low <- c(nodes, interior, random)
  p <- c(
    exp(low),                                   # [P.MIN, 1/2]
    -expm1(low),                                # [1/2, 1 - P.MIN]
    0, 1,
    0.5, 0.5 * (1 - EPS / 2), 0.5 * (1 + EPS),  # 1/2 and its neighbours
    P.MIN * (1 - EPS), P.MIN, P.MIN * (1 + EPS), # the fallback boundary
    1 - EPS / 2, 1 - EPS, 1 - 1e-12, 1 - 1e-8    # next to 1
  )
  sort(unique(p))
}

## Monotonicity violations beyond ~1 ulp on sorted p (direction by tail).
mono_violations <- function(x, lower) {
  d <- diff(x) * if (lower) 1 else -1
  tol <- EPS * pmax(abs(x[-1L]), abs(x[-length(x)]))
  sum(d < -tol)
}

## ------------------------------------------------------------- the runs
GRIDS <- c(1024L, 2048L, 4096L, 8192L, 16384L)
VARIANTS <- list(
  cubic0 = list(variant = "cubic", newton = FALSE),
  cubic1 = list(variant = "cubic", newton = TRUE),
  quintic0 = list(variant = "quintic", newton = FALSE),
  quintic1 = list(variant = "quintic", newton = TRUE)
)
comps <- list(U = c(gold$a[1L], gold$b[1L]), M = c(gold$a[3L], gold$b[3L]))
s1 <- Filter(function(s) s$sample == 1L, sites)
n.s1 <- sum(vapply(s1, function(s) length(s$u), 1L))

## Rmath baseline.
rm.us <- bench::mark(
  for (s in s1) qbeta(s$p, gold$a[if (s$comp == "U") 1L else 3L], gold$b[if (s$comp == "U") 1L else 3L], lower.tail = s$lower),
  min_iterations = 20L, check = FALSE
)
rm.map <- bench::mark(
  for (s in s1) qbeta(pbeta(s$x, s$a.s, s$b.s, lower.tail = s$lower),
                      gold$a[if (s$comp == "U") 1L else 3L], gold$b[if (s$comp == "U") 1L else 3L], lower.tail = s$lower),
  min_iterations = 20L, check = FALSE
)
rows <- list(data.frame(
  variant = "rmath", n = NA_integer_, real_max = 0, real_bitid = 100, real_1ulp = 100,
  sweep_max = 0, sweep_worst = "", mono = NA_integer_, fallback = 0L,
  us_elem = as.numeric(rm.us$median) * 1e6 / n.s1, map_ms = as.numeric(rm.map$median) * 1e3,
  build_ms = 0, size_kb = 0, stringsAsFactors = FALSE
))

## Rmath's own monotonicity count on the sweep, for the record.
rm.mono <- 0L
for (cn in names(comps)) {
  hf0 <- make_half(comps[[cn]][1L], comps[[cn]][2L], 4096L)
  p <- make_sweep_p(hf0)
  for (lower in c(TRUE, FALSE)) {
    rm.mono <- rm.mono + mono_violations(qbeta(p, comps[[cn]][1L], comps[[cn]][2L], lower.tail = lower), lower)
  }
}
rows[[1L]]$mono <- rm.mono
rows[[1L]]$polish <- NA_integer_

## Node diagnostic at n = 4096: how well Rmath's qbeta solves log F(x) = u at
## the nodes (log.p = TRUE vs probability space), and after polishing.
cat("node diagnostic (n = 4096), max |pbeta(x, log.p = TRUE) - u| over the nodes:\n")
for (cn in names(comps)) for (swap in c(FALSE, TRUE)) {
  a <- comps[[cn]][if (swap) 2L else 1L]; b <- comps[[cn]][if (swap) 1L else 2L]
  u <- seq(log(P.MIN), log(0.5), length.out = 4096L)
  x.log <- qbeta(u, a, b, log.p = TRUE)
  x.p <- qbeta(exp(u), a, b)
  resid <- function(x) max(abs(pbeta(x, a, b, log.p = TRUE) - u))
  lb <- lbeta(a, b)
  x.pol <- x.log
  for (z in 1:2) {
    lF <- pbeta(x.pol, a, b, log.p = TRUE)
    lf <- (a - 1) * log(x.pol) + (b - 1) * log1p(-x.pol) - lb
    x.pol <- x.pol - (lF - u) * exp(lF - lf)
  }
  cat(sprintf("  %s %-10s log.p nodes %.1e  prob-space nodes %.1e  polished x2 %.1e  max|x.log - x.p| %.1e  max rel|x.log - x.pol| %.1e\n",
              cn, if (swap) "Beta(b,a)" else "Beta(a,b)", resid(x.log), resid(x.p), resid(x.pol),
              max(abs(x.log - x.p)), max(abs(x.log - x.pol) / x.pol)))
}
cat("\n")

for (polish in c(0L, 2L)) for (n in GRIDS) {
  build <- system.time(for (z in 1:5) tabs <- lapply(comps, function(ab) make_qtable(ab[1L], ab[2L], n, polish)))[["elapsed"]] / 5 * 1e3
  sweeps <- lapply(names(comps), function(cn) {
    p <- make_sweep_p(tabs[[cn]]$ab)
    lapply(c(lower = TRUE, upper = FALSE), function(lower) list(
      comp = cn, lower = lower, p = p, u = log(p),
      ref = qbeta(p, comps[[cn]][1L], comps[[cn]][2L], lower.tail = lower)
    ))
  })
  sweeps <- unlist(sweeps, recursive = FALSE)

  for (vn in names(VARIANTS)) {
    vr <- VARIANTS[[vn]]
    ## real inputs
    diffs <- unlist(lapply(sites, function(s) {
      qbeta_table(s$u, tabs[[s$comp]], s$lower, vr$variant, vr$newton) - s$ref
    }))
    refs <- unlist(lapply(sites, `[[`, "ref"))
    real_max <- max(abs(diffs))
    real_bitid <- 100 * mean(diffs == 0)
    real_1ulp <- 100 * mean(abs(diffs) <= EPS * abs(refs))
    fallback <- sum(unlist(lapply(sites, function(s) s$u > -Inf & s$u <= log(P.MIN))))
    ## sweep
    sw <- vapply(sweeps, function(w) {
      x <- qbeta_table(w$u, tabs[[w$comp]], w$lower, vr$variant, vr$newton)
      c(max = max(abs(x - w$ref)), mono = mono_violations(x, w$lower))
    }, numeric(2L))
    worst <- which.max(sw["max", ])
    ## speed: gold side only (u given) and the full production map, sample 1
    us <- bench::mark(
      for (s in s1) qbeta_table(s$u, tabs[[s$comp]], s$lower, vr$variant, vr$newton),
      min_iterations = 20L, check = FALSE
    )
    mp <- bench::mark(
      for (s in s1) qbeta_table(log(pbeta(s$x, s$a.s, s$b.s, lower.tail = s$lower)), tabs[[s$comp]], s$lower, vr$variant, vr$newton),
      min_iterations = 20L, check = FALSE
    )
    n.arr <- if (vr$variant == "cubic") 4L else 6L
    rows[[length(rows) + 1L]] <- data.frame(
      variant = vn, n = n, polish = polish, real_max = real_max, real_bitid = real_bitid, real_1ulp = real_1ulp,
      sweep_max = max(sw["max", ]),
      sweep_worst = paste0(sweeps[[worst]]$comp, if (sweeps[[worst]]$lower) "-lower" else "-upper"),
      mono = as.integer(sum(sw["mono", ])), fallback = fallback,
      us_elem = as.numeric(us$median) * 1e6 / n.s1, map_ms = as.numeric(mp$median) * 1e3,
      build_ms = build, size_kb = 8 * (n - 1) * n.arr * 4 / 1024,
      stringsAsFactors = FALSE
    )
    cat(sprintf("polish=%d n=%5d %-8s real_max %.1e sweep_max %.1e (%s) mono %d us/elem %.2f map %.1f ms\n",
                polish, n, vn, real_max, max(sw["max", ]), rows[[length(rows)]]$sweep_worst,
                sum(sw["mono", ]), rows[[length(rows)]]$us_elem, rows[[length(rows)]]$map_ms))
  }
}

res <- do.call(rbind, rows)
cat("\nsweep points per (component, tail):", length(sweeps[[1L]]$p), "; sample 1 gold-side evaluations:", n.s1, "\n\n")
print(res, row.names = FALSE)

## Markdown for the hand-off file.
fmt <- function(v, f) ifelse(is.na(v), "", sprintf(f, v))
md <- with(res, sprintf("| %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s |",
  variant, fmt(n, "%d"), fmt(polish, "%d"), fmt(real_max, "%.1e"), fmt(real_bitid, "%.1f"), fmt(real_1ulp, "%.1f"),
  fmt(sweep_max, "%.1e"), sweep_worst, fmt(mono, "%d"), fmt(us_elem, "%.2f"), fmt(map_ms, "%.1f"),
  fmt(build_ms, "%.0f"), fmt(size_kb, "%.0f")))
cat("\n| variant | n | polish | real max | bit-id % | <=1 ulp % | sweep max | worst | mono | us/elem | map ms | build ms | size KB |\n",
    "|---|---|---|---|---|---|---|---|---|---|---|---|---|\n", paste(md, collapse = "\n"), "\n", sep = "")
