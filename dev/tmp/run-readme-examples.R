devtools::load_all(".")
data(horvath_goldstandard)
data(GPL21145_sample)

cat("=== dim ===\n")
print(dim(GPL21145_sample))

cat("=== gold head ===\n")
print(head(horvath_goldstandard$goldstandard2, 3))

fit <- bmiq_calibration(
  datM = GPL21145_sample,
  goldstandard.beta = horvath_goldstandard$goldstandard2,
  verbose = FALSE
)
cat("=== fit ===\n")
print(fit)

calibrated <- as.matrix(fit)
cat("=== dim calibrated ===\n")
print(dim(calibrated))

fit2 <- bmiq_calibration(
  datM = GPL21145_sample,
  goldstandard.beta = horvath_goldstandard$goldstandard2,
  nL = 2L,
  verbose = FALSE
)
cat("=== fit2 ===\n")
print(fit2)

set.seed(1)
X <- matrix(rnorm(5 * 100), nrow = 5)
target <- sort(rnorm(100))
Y <- quantile_norm(X, target)
cat("=== Y dim ===\n")
print(dim(Y))
