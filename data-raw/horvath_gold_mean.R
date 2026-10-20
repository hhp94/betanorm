## code to prepare `horvath_goldstandard` dataset goes here
library(readr)

if (!file.exists("data-raw/probeAnnotation21kdatMethUsed.csv.gz")) {
  R.utils::gzip("data-raw/probeAnnotation21kdatMethUsed.csv", remove = TRUE)
}
# probeAnnotation21kdatMethUsed.csv: published with the R tutorial in the
# additional files of Horvath (2013), doi:10.1186/gb-2013-14-10-r115.
# Stored as a plain data frame: a readr tibble carries a column spec and a
# `problems` pointer that is dead once the .rda is loaded.
horvath_goldstandard <- as.data.frame(read_csv(
  "data-raw/probeAnnotation21kdatMethUsed.csv.gz",
  show_col_types = FALSE
))
stopifnot(isTRUE(anyDuplicated(horvath_goldstandard$Name) == 0))
usethis::use_data(horvath_goldstandard, compress = "xz", overwrite = TRUE)
