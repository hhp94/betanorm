## code to prepare `GPL21145_sample` dataset goes here
library(dplyr)
library(readr)

if (!file.exists("data-raw/server_matrix.csv.gz")) {
  R.utils::gzip("data-raw/server_matrix.csv", remove = TRUE)
}

GPL21145_sample_raw <- read_csv("data-raw/server_matrix.csv.gz")

stopifnot(isTRUE(anyDuplicated(horvath_goldstandard$Name) == 0))

set.seed(1234)

# pick 6 random columns (besides ProbeID) to keep
other_cols <- sample(setdiff(names(GPL21145_sample_raw), "ProbeID"), 6)

GPL21145_sample <- GPL21145_sample_raw |>
  select(ProbeID, all_of(other_cols)) |>
  filter(ProbeID %in% horvath_goldstandard$Name)

stopifnot(isTRUE(all.equal(GPL21145_sample$ProbeID, horvath_goldstandard$Name)))

# convert to matrix: drop ProbeID column, use it as rownames, then transpose
GPL21145_sample <- GPL21145_sample |>
  as.data.frame() |>
  tibble::column_to_rownames("ProbeID") |>
  as.matrix() |>
  t()

# absolute 1s
usethis::use_data(GPL21145_sample, overwrite = TRUE)
