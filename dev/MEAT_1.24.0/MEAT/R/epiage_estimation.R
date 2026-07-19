#' Estimates age in skeletal muscle from calibrated DNA methylation profiles.
#'
#' \code{epiage_estimation} takes as input a
#' \code{\link[SummarizedExperiment]{SummarizedExperiment-class}} object whose
#' assays contain a beta-matrix called "beta". This beta-matrix should contain
#' DNA methylation profiles in skeletal muscle that have been cleaned with
#' \code{\link{clean_beta}} and calibrated with \code{\link{BMIQcalibration}}.
#' \code{epiage_estimation} will use the muscle clock to estimate epigenetic age
#' in each sample.
#'
#' \code{epiage_estimation} estimates epigenetic age for each sample in the
#' input \code{SE} based on DNA methylation profiles. \code{SE} needs to be a
#' \code{\link[SummarizedExperiment]{SummarizedExperiment-class}} object containing
#' a matrix of beta-values called "beta" in assays. Beta must have been
#' calibrated to the gold standard GSE50498 using \code{\link{BMIQcalibration}}
#' to obtain good estimates of epigenetic age.
#'
#' @param SE A \code{\link[SummarizedExperiment]{SummarizedExperiment-class}} object.
#' The "assays" component of \code{SE} should contain a beta-matrix of
#' DNA methylation beta-values called "beta" that has been cleaned with
#' \code{\link{clean_beta}} and calibrated with \code{\link{BMIQcalibration}}.
#' \code{SE} may optionally contain annotation information on the CpGs
#' stored in "rowData" and sample phenotypes stored in "colData".
#' @param version A character specifying which version of the epigenetic clock
#' you would like to use. Dy default, \code{version} is set to "MEAT2.0" for the
#' second version of the epigenetic clock. If you would like to use the original
#' version, set \code{version} to "MEAT".
#' @param age_col_name The name of the column in colData from \code{SE} that
#' contains age (in years).
#'
#' @return A \code{\link[SummarizedExperiment]{SummarizedExperiment-class}} object
#' identical to the input \code{SE}, with components added to colData. If no
#' phenotypes were provided in the colData of the input \code{SE},
#' \code{epiage_estimation} will put in colData a tibble containing a single
#' column called "DNAmage", corresponding to epigenetic age (in years) for each
#' sample. If phenotypes were provided in the colData of the input \code{SE},
#' \code{epiage_estimation} will add to the existing colData three columns:
#' \enumerate{
#'   \item \code{DNAmage} epigenetic age (in years)
#'   \item \code{AAdiff} the difference between predicted and actual age
#'   (in years).
#'   \item \code{AAresid} the residuals of a linear model
#'   (using \code{\link[stats]{lm}}) of DNAmage against actual age.
#'   \code{AAresid} is only returned if the number of samples is > 2, as
#'   \code{AAresid} cannot be calculated with < 2 samples.
#' }
#' @export
#' @import SummarizedExperiment
#' @import glmnet
#' @importFrom dplyr pull
#' @import tibble
#' @seealso \code{\link[wateRmelon]{BMIQ}} for the original BMIQ algorithm,
#' \url{https://genomebiology.biomedcentral.com/articles/10.1186/gb-2013-14-10-r115}
#' for the adapted version of the BMIQ algorithm, and
#' \url{https://onlinelibrary.wiley.com/doi/full/10.1002/jcsm.12556}
#' for the elastic net model of the muscle clock.
#' @examples
#' # Load matrix of beta-values of two individuals from dataset GSE121961
#' data("GSE121961", envir = environment())
#' # Load phenotypes of the two individuals from dataset GSE121961
#' data("GSE121961_pheno", envir = environment())
#'
#' # Create a SummarizedExperiment object to coordinate phenotypes and
#' # methylation into one object.
#' library(SummarizedExperiment)
#' GSE121961_SE <- SummarizedExperiment(assays=list(beta=GSE121961),
#' colData=GSE121961_pheno)
#'
#' # Run clean_beta() to clean the beta-matrix
#' GSE121961_SE_clean <- clean_beta(SE = GSE121961_SE,
#' version = "MEAT2.0")
#'
#' # Run BMIQcalibration() to calibrate the clean beta-matrix
#' GSE121961_SE_calibrated <- BMIQcalibration(SE = GSE121961_SE_clean,
#' version = "MEAT2.0")
#'
#' # Run epiage_estimation() to obtain DNAmage + optionally AAdiff and AAresid
#' GSE121961_SE_epiage <- epiage_estimation(SE = GSE121961_SE_calibrated,
#' version = "MEAT2.0",
#' age_col_name = "Age")
#' colData(GSE121961_SE_epiage)
epiage_estimation <- function(SE = NULL,
                              version = "MEAT2.0",
                              age_col_name = NULL) {

  # Check whether SE is a SummarizedExperiment object
  if (!is(SE, "SummarizedExperiment"))
    stop("Please make sure SE is a SummarizedExperiment object.")

  # Check that the beta-matrix is indeed called "beta"
  if (names(assays(SE))!="beta")
    stop("Please make sure that the beta-matrix stored in the assays component of SE is called beta.")

  # Check version is correct
  if (!version %in% c("MEAT","MEAT2.0")) {
    stop("Please provide a valid version for the muscle clock. Set version to either MEAT or MEAT2.0")
  }

  # Load the elastic net model
  elasticnet_model_MEAT <- NULL
  elasticnet_model_MEAT2.0 <- NULL
  if (version=="MEAT")
  {
    data("elasticnet_model_MEAT",envir = environment())
    elasticnet_model <- elasticnet_model_MEAT
    lambda.glmnet.Training <- 0.025
  }

  else
  {
    data("elasticnet_model_MEAT2.0",envir = environment())
    elasticnet_model <- elasticnet_model_MEAT2.0
    lambda.glmnet.Training <- 0.033
  }


  # Predict age based on calibrated DNA methylation profile
  DNAmage <- NULL
  print(dim)
  DNAmage <- anti.trafo(predict(elasticnet_model,
                                t(assays(SE)$beta),
                                type = "response",
                                s = lambda.glmnet.Training))[,1]

  if (ncol(colData(SE))==0) {
    AAdiff <- NULL
    AAresid <- NULL
    SE2 <- SummarizedExperiment(assays = assays(SE),
                               rowData = rowData(SE),
                               colData = as.data.frame(DNAmage))
    message("There are no phenotypes provided in colData, so only DNAmage will
            be returned.")
  } else {
    pheno <- colData(SE)

    # Check that the column name for age does exist in colData
    if (!age_col_name %in% colnames(pheno))
      stop(paste0("colData does not contain a column called ",age_col_name,
                     ". Please check the column names of colData."))

    pheno <- as_tibble(pheno)
    Age <- pull(pheno[, age_col_name])
    AAdiff <- DNAmage - Age
    pheno <- data.frame(as.data.frame(pheno),
                        DNAmage,
                        AAdiff)
    if (nrow(pheno)>2)
    {
      AAresid <- AAdiff
      AAresid_noNA <- resid(lm(DNAmage ~ Age))
      AAresid[!is.na(AAresid)] <- AAresid_noNA
      pheno <- data.frame(pheno,
                          AAresid)
    }
    else
    {
      message("You only have two samples in SE, so only DNAmage and AAdiff
              will be returned.")
    }
    SE2 <- SummarizedExperiment(assays = assays(SE),
                                  rowData = rowData(SE),
                                  colData = pheno)
  }
  return(SE2)
}

anti.trafo <- function(x, adult.age = 20) {
  ifelse(x < 0, (1 + adult.age) * exp(x) - 1, (1 + adult.age) * x + adult.age)
}
