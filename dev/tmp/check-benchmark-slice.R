devtools::load_all(".", quiet = TRUE)
source("dev/bmiq-calibration-legacy.R", local = TRUE)
data("GPL21145_sample", envir = environment())
data("horvath_goldstandard", envir = environment())

cn <- colnames(GPL21145_sample)
nm <- as.character(horvath_goldstandard$Name)
idx <- match(cn, nm)
cat("setequal probes:", setequal(cn, nm), "\n")
cat("identical order:", identical(cn, nm), "\n")
cat(
  "match all found:",
  !any(is.na(idx)),
  " identity idx:",
  identical(idx, seq_along(cn)),
  "\n"
)

# Always name-align gold (length equality alone is not enough).
gold_full <- as.numeric(horvath_goldstandard$goldstandard2)[idx]

for (sample_i in 1:2) {
  cat(
    "\n======== sample_i =",
    sample_i,
    rownames(GPL21145_sample)[sample_i],
    "========\n"
  )
  dat1 <- GPL21145_sample[sample_i, , drop = FALSE]
  gold <- gold_full
  ok <- is.finite(dat1[1L, ]) &
    is.finite(gold) &
    dat1[1L, ] >= 0 &
    dat1[1L, ] <= 1 &
    gold >= 0 &
    gold <= 1
  cat("ok count:", sum(ok), "/", length(ok), "\n")

  dat1 <- dat1[, ok, drop = FALSE]
  gold_v <- as.numeric(gold[ok])
  before_v <- as.numeric(dat1[1L, ])
  storage.mode(dat1) <- "double"

  # Slicing checks
  direct <- as.numeric(GPL21145_sample[sample_i, ok])
  cat("before == direct row[ok]:", isTRUE(all.equal(before_v, direct)), "\n")
  cat("length before/gold:", length(before_v), length(gold_v), "\n")
  cat(
    "means before/gold:",
    mean(before_v),
    mean(gold_v),
    " cor:",
    cor(before_v, gold_v),
    "\n"
  )

  # Density peaks (what the plot shows)
  db <- density(before_v, from = 0, to = 1, n = 512)
  dg <- density(gold_v, from = 0, to = 1, n = 512)
  cat(sprintf(
    "density peak before: y=%.2f at x=%.3f | gold: y=%.2f at x=%.3f\n",
    max(db$y),
    db$x[which.max(db$y)],
    max(dg$y),
    dg$x[which.max(dg$y)]
  ))

  nfit <- 5000L
  niter <- 25L
  res2 <- bmiq_calibration(
    dat1,
    gold_v,
    nL = 2L,
    nfit = nfit,
    niter = niter,
    verbose = FALSE,
    debug = TRUE
  )
  res3 <- bmiq_calibration(
    dat1,
    gold_v,
    nL = 3L,
    nfit = nfit,
    niter = niter,
    verbose = FALSE,
    debug = TRUE
  )
  leg <- BMIQcalibration.old(
    datM = dat1,
    goldstandard.beta = gold_v,
    nL = 3,
    doH = TRUE,
    nfit = nfit,
    niter = niter,
    tol = 0.001
  )

  c2 <- as.numeric(res2$calibrated[1L, ])
  c3 <- as.numeric(res3$calibrated[1L, ])
  cat("legacy class/dim:", class(leg), paste(dim(leg), collapse = "x"), "\n")
  cl <- as.numeric(leg[1L, ])

  cat("lengths c2/c3/cl:", length(c2), length(c3), length(cl), "\n")
  cat("success nL2/nL3:", res2$success, res3$success, "\n")
  cat("anyNA c2/c3/cl:", anyNA(c2), anyNA(c3), anyNA(cl), "\n")
  cat(
    "means c2/c3/cl:",
    mean(c2, na.rm = TRUE),
    mean(c3, na.rm = TRUE),
    mean(cl, na.rm = TRUE),
    "\n"
  )

  # dens_df label alignment (the plot-weird bug if lengths differ)
  vals <- c(gold_v, before_v, c2, c3, cl)
  labs <- rep(
    c("gold", "before", "nL=2", "nL=3", "legacy"),
    each = length(before_v)
  )
  cat(
    "dens_df lens match:",
    length(vals) == length(labs),
    " vals=",
    length(vals),
    " labs=",
    length(labs),
    "\n"
  )

  # Per-type means via the same construction as the plot
  dens_df <- data.frame(value = vals, type = labs)
  print(aggregate(value ~ type, dens_df, function(z) {
    c(n = length(z), mean = mean(z), maxd = max(density(z, from = 0, to = 1)$y))
  }))

  if (!is.null(res3$diagnostics$samples[[1L]]$thresholds)) {
    cat("nL3 thresholds:", res3$diagnostics$samples[[1L]]$thresholds, "\n")
  }
  if (!is.null(res2$diagnostics$samples[[1L]]$thresholds)) {
    cat("nL2 thresholds:", res2$diagnostics$samples[[1L]]$thresholds, "\n")
  }
}
