#' Horvath Gold-Standard Probe Annotation
#'
#' Probe annotation and gold-standard mean methylation values for Horvath-style
#' BMIQ calibration (about 21k probes).
#'
#' @format A data frame with 21368 rows and 7 variables:
#' \describe{
#'   \item{`Name`}{Probe identifier.}
#'   \item{`Gene_ID`}{Associated gene identifier.}
#'   \item{`GenomeBuild`}{NCBI genome build (36 for every probe).}
#'   \item{`Chr`}{Chromosome.}
#'   \item{`MapInfo`}{Genomic coordinate on build 36, not on a current
#'     build.}
#'   \item{`SourceVersion`}{Annotation version.}
#'   \item{`goldstandard2`}{Gold-standard mean beta value.}
#' }
#'
#' @examples
#' data(horvath_goldstandard)
#' head(horvath_goldstandard)
#'
#' @source `probeAnnotation21kdatMethUsed.csv`, published with the R tutorial
#'   in the additional files of Horvath (2013) under CC BY 2.0
#'   (\doi{10.1186/gb-2013-14-10-r115}). The values are unchanged.
#'
#' @references
#' Horvath S (2013).
#' DNA methylation age of human tissues and cell types.
#' *Genome Biology* 14(10), R115.
#' \doi{10.1186/gb-2013-14-10-r115}
#'
#' @keywords datasets
"horvath_goldstandard"


#' Example Methylation Beta Matrix (GPL21145)
#'
#' Example DNA methylation beta matrix for probes in
#' [horvath_goldstandard], for use with [bmiq_calibration()].
#'
#' @format A numeric matrix with 6 samples (rows) and 21368 CpGs (columns).
#'   Column names match `horvath_goldstandard[["Name"]]`. Values are in
#'   \eqn{[0, 1]}. The 1062 probes of [horvath_goldstandard] that the EPIC
#'   array does not measure are filled with 0 in every sample. They are not
#'   measurements, so drop them before calibrating (see the examples); the
#'   gold vector need not match the columns.
#'
#' @examples
#' data(GPL21145_sample)
#' dim(GPL21145_sample)
#'
#' # Probes the EPIC array does not measure are 0 in every sample.
#' measured <- colSums(GPL21145_sample == 0) < nrow(GPL21145_sample)
#' sum(!measured)
#' dim(GPL21145_sample[, measured])
#'
#' @source NCBI GEO, platform GPL21145 (Illumina MethylationEPIC):
#'   <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL21145>.
#'
#' @keywords datasets
"GPL21145_sample"
