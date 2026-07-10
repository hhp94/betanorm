load("data/GPL21145_sample.rda")
load("data/horvath_goldstandard.rda")

dat <- GPL21145_sample
gold <- as.numeric(horvath_goldstandard[["goldstandard2"]])
b <- as.numeric(dat[1, ])
th <- c(0.2, 0.75)

dir.create("tmp-plots", showWarnings = FALSE)

png("tmp-plots/beta-densities.png", width = 1100, height = 1000, res = 120)
par(mfrow = c(3, 1), mar = c(4, 4, 3, 1))
plot(
  density(gold, from = 0, to = 1, n = 512),
  main = "Gold standard (n=21368)",
  xlab = "beta",
  xlim = c(0, 1),
  lwd = 2
)
abline(v = th, col = "red", lty = 2)
legend(
  "topright",
  c("density", "default th1.v"),
  col = c("black", "red"),
  lty = c(1, 2),
  bty = "n"
)

plot(
  density(b, from = 0, to = 1, n = 512),
  main = paste0("Sample 1: ", rownames(dat)[1]),
  xlab = "beta",
  xlim = c(0, 1),
  lwd = 2
)
abline(v = th, col = "red", lty = 2)

hist(
  b,
  breaks = seq(0, 1, by = 0.01),
  main = "Sample 1 histogram (1% bins)",
  xlab = "beta",
  xlim = c(0, 1),
  col = "grey80",
  border = "grey60",
  freq = FALSE
)
lines(density(b, from = 0, to = 1, n = 512), lwd = 2)
abline(v = th, col = "red", lty = 2)
dev.off()

message("wrote tmp-plots/beta-densities.png")

d <- density(b, from = 0, to = 1, n = 512)
message(
  "sample1 density mode approx at ",
  round(d$x[which.max(d$y)], 4),
  " (max dens ",
  round(max(d$y), 3),
  ")"
)

bands <- list(c(0.2, 0.75), c(0.25, 0.6), c(0.3, 0.55), c(0.35, 0.5))
for (band in bands) {
  nmid <- sum(b > band[1] & b < band[2])
  message(sprintf(
    "count in (%.2f, %.2f): %d (%.1f%%)",
    band[1],
    band[2],
    nmid,
    100 * nmid / length(b)
  ))
}

set.seed(1)
km <- kmeans(b, centers = 2, nstart = 10)
cl <- km$cluster
u_lab <- which.min(km$centers)
m_lab <- which.max(km$centers)
gap_lo <- max(b[cl == u_lab])
gap_hi <- min(b[cl == m_lab])
message(
  "2-means centers: ",
  paste(round(sort(as.numeric(km$centers)), 3), collapse = ", ")
)
message(
  "between-cluster gap (max low-cl, min high-cl): ",
  round(gap_lo, 4),
  " .. ",
  round(gap_hi, 4)
)
message(
  "probes strictly between that gap: ",
  sum(b > gap_lo & b < gap_hi)
)

# Illustrate legacy empty-middle fallback formula
# Suppose posterior hard classes on a subsample have no class 2
# (all assigned 1 or 3). Then:
fake_u <- b[b < 0.3]
fake_m <- b[b > 0.7]
th_fallback <- c(
  0.5 * max(fake_u) + 0.5 * mean(fake_m),
  (1 / 3) * max(fake_u) + (2 / 3) * mean(fake_m)
)
message(
  "legacy-style fallback thresholds if middle empty: ",
  paste(round(th_fallback, 4), collapse = ", ")
)
cls <- rep(2L, length(b))
cls[b <= th_fallback[1]] <- 1L
cls[b >= th_fallback[2]] <- 3L
message(
  "full-sample class counts under that fallback: ",
  paste(tabulate(cls, 3), collapse = ", ")
)
