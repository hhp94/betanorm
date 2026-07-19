#' Quantile Normalization to a Target Distribution
#'
#' Minimal, no-missing-data port of
#' [preprocessCore::normalize.quantiles.use.target()]. Each sample row is
#' ranked and mapped onto `target` (Bolstad equal-length indexing when
#' `length(target) == ncol(obj)`; linear quantile interpolation otherwise).
#'
#' @param obj Numeric matrix with samples in rows and variables in columns.
#'   Must be finite with no missing values.
#' @param target Numeric vector giving the reference distribution. Must be
#'   finite with no missing values; length may differ from `ncol(obj)`.
#'   Need not be sorted (a sorted copy is made internally).
#'
#' @return A numeric matrix the same dimension as `obj`, with `dimnames`
#'   preserved when present.
#'
#' @references
#' Bolstad BM, Irizarry RA, Åstrand M, Speed TP (2003).
#' A comparison of normalization methods for high density oligonucleotide
#' array data based on variance and bias.
#' *Bioinformatics* 19(2), 185–193.
#' \doi{10.1093/bioinformatics/19.2.185}
#'
#' @export
quantile_norm <- function(obj, target) {
  checkmate::assert_matrix(
    obj,
    mode = "numeric",
    any.missing = FALSE,
    min.rows = 1L,
    min.cols = 1L
  )
  checkmate::assert_numeric(
    target,
    any.missing = FALSE,
    min.len = 1L,
    finite = TRUE,
    .var.name = "target"
  )
  qnorm_target_rows_cpp(obj, target)
}
