#' Horvath Gold-Standard Probe Annotation
#'
#' Probe annotation and gold-standard mean methylation values for Horvath-style
#' BMIQ calibration (about 21k probes).
#'
#' @format A data frame with 21368 rows and 7 variables:
#' \describe{
#'   \item{`Name`}{Probe identifier.}
#'   \item{`Gene_ID`}{Associated gene identifier (may be missing).}
#'   \item{`GenomeBuild`}{Genome build.}
#'   \item{`Chr`}{Chromosome.}
#'   \item{`MapInfo`}{Genomic coordinate.}
#'   \item{`SourceVersion`}{Annotation version.}
#'   \item{`goldstandard2`}{Gold-standard mean beta value.}
#' }
#'
#' @seealso [BMIQcalibration()], [GPL21145_sample]
#'
#' @examples
#' data(horvath_goldstandard)
#' head(horvath_goldstandard)
#'
#' @keywords datasets
"horvath_goldstandard"


#' Example Methylation Beta Matrix (GPL21145)
#'
#' Example DNA methylation beta matrix for probes in
#' [horvath_goldstandard], for use with [BMIQcalibration()].
#'
#' @format A numeric matrix with 6 samples (rows) and 21368 CpGs (columns).
#'   Column names match `horvath_goldstandard$Name`. Values are in
#'   \eqn{[0, 1]}.
#'
#' @seealso [BMIQcalibration()], [horvath_goldstandard]
#'
#' @examples
#' data(GPL21145_sample)
#' dim(GPL21145_sample)
#'
#' @keywords datasets
"GPL21145_sample"
