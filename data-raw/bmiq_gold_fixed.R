## Writes the gold-standard fits that bmiq_calibration_fixed() carries as
## numbers (FIXED_GOLD in R/bmiq_calibration_fixed.R). Run from the package
## root after any change to the gold data or the EM, and paste the output
## over the FIXED_GOLD block, then run data-raw/bmiq_gold_qtable.R to
## rebuild the gold quantile table from the pasted numbers.
## tests/testthat/test-bmiq-fixed.R checks that the pasted numbers are the
## live bmiq_gold_fit() result and that the table was built for them.
##
## Each value is written as a hexadecimal floating-point literal, which R
## parses exactly on every platform (its decimal parser is not correctly
## rounded on all of them), with the 17-digit decimal beside it to read.
devtools::load_all(".", quiet = TRUE)

gold.sets <- list(
  # Horvath (2013), probeAnnotation21kdatMethUsed.csv; see
  # data-raw/horvath_gold_mean.R.
  horvath = horvath_goldstandard[["goldstandard2"]]
)

literal <- function(x) {
  # A comma after every element but the last, then the decimal as a comment.
  sprintf(
    "%s%s # %.17g",
    sprintf("%a", x),
    c(rep(",", length(x) - 1L), ""),
    x
  )
}
vector_lines <- function(name, x, indent) {
  c(
    paste0(indent, name, " = c("),
    paste0(indent, "  ", literal(x)),
    paste0(indent, "),")
  )
}
scalar_line <- function(name, x, indent, last = FALSE) {
  paste0(indent, name, " = ", sprintf("%a", x), if (!last) ",", " # ",
    sprintf("%.17g", x))
}

lines <- "FIXED_GOLD <- list("
for (set in names(gold.sets)) {
  fit <- bmiq_gold_fit(gold.sets[[set]], nL = 3L, verbose = FALSE)
  stopifnot(
    identical(fit[["settings"]][["nfit"]], 20000L),
    identical(fit[["settings"]][["niter"]], 5L),
    identical(fit[["settings"]][["th1.v"]], c(0.2, 0.75))
  )
  lines <- c(
    lines,
    paste0("  ", set, " = list("),
    vector_lines("a", fit[["a"]], "    "),
    vector_lines("b", fit[["b"]], "    "),
    vector_lines("thresholds", fit[["thresholds"]], "    "),
    scalar_line("unmethylated.mode", fit[["unmethylated.mode"]], "    "),
    scalar_line("methylated.mode", fit[["methylated.mode"]], "    ",
      last = TRUE),
    "  ),"
  )
  # Round trip: the printed literals must parse back to the same bits.
  parsed <- eval(parse(text = paste(
    c("list(", vector_lines("a", fit[["a"]], ""),
      vector_lines("b", fit[["b"]], ""),
      vector_lines("thresholds", fit[["thresholds"]], ""),
      scalar_line("unmethylated.mode", fit[["unmethylated.mode"]], ""),
      scalar_line("methylated.mode", fit[["methylated.mode"]], "",
        last = TRUE),
      ")"),
    collapse = "\n"
  )))
  stopifnot(identical(parsed, fit[c("a", "b", "thresholds",
    "unmethylated.mode", "methylated.mode")]))
}
lines[length(lines)] <- "  )"
lines <- c(lines, ")")
cat(lines, sep = "\n")
