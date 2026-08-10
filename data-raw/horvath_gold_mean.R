## code to prepare `horvath_goldstandard` dataset goes here
library(readr)

if (!file.exists("data-raw/probeAnnotation21kdatMethUsed.csv.gz")) {
  R.utils::gzip("data-raw/probeAnnotation21kdatMethUsed.csv", remove = TRUE)
}
horvath_goldstandard <- read_csv(
  "data-raw/probeAnnotation21kdatMethUsed.csv.gz"
)
stopifnot(isTRUE(anyDuplicated(horvath_goldstandard$Name) == 0))
usethis::use_data(horvath_goldstandard, compress = "xz", overwrite = TRUE)
