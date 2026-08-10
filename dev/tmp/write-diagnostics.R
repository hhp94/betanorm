dir.create("dev/tmp", recursive = TRUE, showWarnings = FALSE)
devtools::load_all(".", quiet = TRUE)
data("GPL21145_sample", envir = environment())
data("horvath_goldstandard", envir = environment())

sample_i <- 2L
idx <- match(
  colnames(GPL21145_sample),
  as.character(horvath_goldstandard$Name)
)
gold <- as.numeric(horvath_goldstandard$goldstandard2)[idx]
dat1 <- GPL21145_sample[sample_i, , drop = FALSE]
v <- as.numeric(dat1)
ok <- is.finite(v) & is.finite(gold) & v >= 0 & v <= 1 & gold >= 0 & gold <= 1
dat1 <- dat1[, ok, drop = FALSE]
gold_v <- gold[ok]
storage.mode(dat1) <- "double"

strip_random_indices <- function(diagnostics) {
  if (!is.null(diagnostics$gold)) {
    diagnostics$gold$random_indices <- NULL
  }
  if (
    length(diagnostics$samples) >= 1L &&
      !is.null(diagnostics$samples[[1L]])
  ) {
    diagnostics$samples[[1L]]$random_indices <- NULL
  }
  diagnostics
}

for (nL in c(2L, 3L)) {
  message("bmiq_calibration nL=", nL, " ...")
  res <- bmiq_calibration(
    datM = dat1,
    goldstandard.beta = gold_v,
    nL = nL,
    nfit = 20000L,
    niter = 100L,
    tol = 0.001,
    verbose = FALSE,
    debug = TRUE
  )
  d <- strip_random_indices(res$diagnostics)
  f <- file.path("dev/tmp", paste0("benchmark-diagnostics-nl", nL, ".txt"))
  writeLines(
    capture.output(str(d, vec.len = 10L, digits.d = 6L)),
    con = f
  )
  message("Wrote ", f)
}
