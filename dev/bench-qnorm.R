# Benchmark qnorm_target_rows_cpp vs preprocessCore on a large matrix.
#
# Usage (from package root):
#   Rscript dev/bench-qnorm.R
#
# Notes:
# - bmiqpp API is samples x variables; preprocessCore is features x samples.
# - Debug / devtools builds are much slower than -O3 release installs.
# - Set BENCH_ITERS / BENCH_N_SAMPLES / BENCH_N_VARS to override defaults.

devtools::load_all(".", quiet = TRUE)

if (!requireNamespace("preprocessCore", quietly = TRUE)) {
  stop("preprocessCore is required for this benchmark (Suggests).")
}
if (!requireNamespace("bench", quietly = TRUE)) {
  stop("bench is required: install.packages(\"bench\")")
}

library(preprocessCore)
library(bench)

n_samples <- as.integer(Sys.getenv("BENCH_N_SAMPLES", "500"))
n_vars <- as.integer(Sys.getenv("BENCH_N_VARS", "20000"))
iters <- as.integer(Sys.getenv("BENCH_ITERS", "10"))

set.seed(1L)
obj <- matrix(rnorm(n_samples * n_vars), n_samples, n_vars)
target <- sort(rnorm(n_vars))
X_feat <- t(obj)

cat(
  "dim obj (samples x vars): ",
  n_samples,
  " x ",
  n_vars,
  "\n",
  "dim X_feat (feat x samp): ",
  n_vars,
  " x ",
  n_samples,
  "\n",
  "length(target): ",
  length(target),
  "\n",
  "iterations: ",
  iters,
  "\n",
  sep = ""
)

# Spot-check equal-length path (ULP noise ok on unequal path).
o1 <- quantile_norm(obj, target)
o2 <- t(normalize.quantiles.use.target(X_feat, target, copy = TRUE))
cat("max|diff| vs preprocessCore: ", max(abs(o1 - o2)), "\n", sep = "")

bm <- bench::mark(
  bmiqpp = quantile_norm(obj, target),
  preprocessCore = {
    t(normalize.quantiles.use.target(X_feat, target, copy = TRUE))
  },
  preprocessCore_no_t = normalize.quantiles.use.target(
    X_feat,
    target,
    copy = TRUE
  ),
  iterations = iters,
  check = FALSE,
  memory = TRUE
)

print(bm)
invisible(bm)
