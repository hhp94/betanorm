## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  fig.width = 8,
  collapse = TRUE,
  comment = "#>"
)

## ----style, echo = FALSE, results = 'asis'------------------------------------
BiocStyle::markdown()

## ----setup, message=F---------------------------------------------------------
library(knitr)

## ----Information on CpGs in MEAT, eval=FALSE----------------------------------
# data("CpGs_in_MEAT",envir = environment())

## ----Information on CpGs in MEAT2.0, eval=FALSE-------------------------------
# data("CpGs_in_MEAT2.0",envir = environment())

## ----MEAT package installation, eval=FALSE------------------------------------
# if (!requireNamespace("BiocManager", quietly=TRUE))
#     install.packages("BiocManager")
# BiocManager::install("MEAT")

## ----MEAT package installation loading, message=FALSE, warning=FALSE----------
library(MEAT)

## ----Methylation matrix presentation------------------------------------------
data("GSE121961", envir = environment())

## ----Methylation matrix presentation table, echo = FALSE----------------------
kable(
  head(GSE121961),
  caption = "Top rows of the GSE121961 matrix before cleaning and calibration."
)

## ----Phenotype table presentation---------------------------------------------
data("GSE121961_pheno", envir = environment())

## ----Phenotype table presentation table, echo = FALSE-------------------------
kable(GSE121961_pheno, caption = "Phenotypes corresponding to GSE121961.")

## ----Data formatting, message=FALSE, warning=FALSE----------------------------
library(SummarizedExperiment)
GSE121961_SE <- SummarizedExperiment(
  assays = list(beta = GSE121961),
  colData = GSE121961_pheno
)
GSE121961_SE

## ----Data cleaning, message=FALSE, warning=FALSE------------------------------
GSE121961_SE_clean <- clean_beta(SE = GSE121961_SE, version = "MEAT2.0")

## ----Data cleaning table, echo = FALSE----------------------------------------
kable(
  head(assays(GSE121961_SE_clean)$beta),
  caption = "Top rows of the GSE121961 beta matrix after cleaning."
)

## ----Data calibration, message=FALSE, warning=FALSE---------------------------
GSE121961_SE_calibrated <- BMIQcalibration(
  SE = GSE121961_SE_clean,
  version = "MEAT2.0"
)

## ----Data calibration table, echo=FALSE---------------------------------------
kable(
  head(assays(GSE121961_SE_calibrated)$beta),
  caption = "Top rows of the GSE121961 beta matrix after cleaning and calibration."
)

## ----DNA methylation profile distribution before and after calibration, message=FALSE, warning=FALSE----
data("gold.mean.MEAT2.0", envir = environment())
GSE121961_SE_clean_with_gold_mean <- cbind(
  assays(GSE121961_SE_clean)$beta,
  gold.mean.MEAT2.0$gold.mean
) # add the gold mean
GSE121961_SE_calibrated_with_gold_mean <- cbind(
  assays(GSE121961_SE_calibrated)$beta,
  gold.mean.MEAT2.0$gold.mean
) # add the gold mean
groups <- c(rep("GSE121961", ncol(GSE121961_SE_clean)), "Gold mean")

library(minfi)
par(mfrow = c(2, 1))
densityPlot(
  GSE121961_SE_clean_with_gold_mean,
  sampGroups = groups,
  main = "Before calibration",
  legend = FALSE
)
densityPlot(
  GSE121961_SE_calibrated_with_gold_mean,
  sampGroups = groups,
  main = "After calibration"
)

## ----Epigenetic age estimation with phenotypes, message=FALSE, warning=FALSE----
GSE121961_SE_epiage <- epiage_estimation(
  SE = GSE121961_SE_calibrated,
  age_col_name = "Age",
  version = "MEAT2.0"
)

## ----Epigenetic age estimation with phenotypes table, echo=FALSE--------------
kable(
  colData(GSE121961_SE_epiage),
  caption = "Phenotypes corresponding to GSE121961 with AAdiff for each sample."
)

## ----session info-------------------------------------------------------------
sessionInfo()
