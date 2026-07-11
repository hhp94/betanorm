#' Calibrate methylation data to a gold standard.
#'
#' \code{BMIQcalibration} uses an adapted version of the BMIQ algorithm to
#' calibrate the beta-matrix stored in the input SummarizedExperiment object
#' \code{SE} to the gold standard dataset used in the muscle clock (GSE50498).
#'
#' \code{BMIQcalibration} was created by Steve Horvath,
#' largely based on the \code{\link[wateRmelon]{BMIQ}} function from
#' Teschendorff (2013) to adjust for the type-2 bias in Illumina HM450
#' and HMEPIC arrays. BMIQ stands for beta mixture quantile normalization.
#' Horvath fixed minor errors in the v_1.2 version of the BMIQ algorithm
#' and changed the optimization algorithm to make the code more robust.
#' He used method = "Nelder-Mead" in \code{\link[stats]{optim}} since
#' the other optimization method sometimes gets stuck. Toward this end,
#' the function \code{\link[RPMM]{blc}} was replaced by \code{blc2}.
#' \code{SE} needs to be a SummarizedExperiment object containing a matrix of
#' beta-values that has been cleaned using \code{\link{clean_beta}}.
#' Each sample in \code{SE} is iteratively calibrated to the
#' gold standard values, so the time it takes to run
#' \code{BMIQcalibration} is directly proportional to the number
#' of samples in \code{SE}. This step is essential to estimate
#' epigenetic age with accuracy.
#' @param SE A \code{\link[SummarizedExperiment]{SummarizedExperiment-class}} object.
#' The "assays" component of \code{SE} should contain a beta-matrix of
#' DNA methylation beta-values called "beta" that has been cleaned with
#' \code{\link{clean_beta}}.
#' \code{SE} may optionally contain annotation information on the CpGs stored
#' in "rowData" and sample phenotypes stored in "colData".
#' @param version A character specifying which version of the epigenetic clock
#' you would like to use. Dy default, \code{version} is set to "MEAT2.0" for the
#' second version of the epigenetic clock. If you would like to use the original
#' version, set \code{version} to "MEAT".
#' @return A calibrated version of the input \code{SE} calibrated to the gold
#' standard dataset GSE50498.
#' @export
#' @import SummarizedExperiment
#' @importFrom stats density pbeta qbeta lm coef optim dbeta predict resid lm
#' @importFrom dynamicTreeCut printFlush
#' @importFrom wateRmelon BMIQ
#' @import RPMM
#' @import grDevices
#' @import graphics
#' @import minfi
#' @seealso \code{\link{clean_beta}} to get the DNA methylation matrix ready
#' for calibration,
#' \code{\link[wateRmelon]{BMIQ}} for the original BMIQ algorithm and
#' \url{https://genomebiology.biomedcentral.com/articles/10.1186/gb-2013-14-10-r115}
#' for the original paper describing Horvath's adapted BMIQ algorithm, and
#' \code{\link[SummarizedExperiment]{SummarizedExperiment-class}} for more
#' details on how to create and manipulate SummarizedExperiment objects.
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
#' GSE121961_SE_clean <- clean_beta(SE = GSE121961_SE, version = "MEAT2.0")
#'
#' # Run BMIQcalibration() to calibrate the clean beta-matrix
#' GSE121961_SE_calibrated <- BMIQcalibration(SE = GSE121961_SE_clean, version = "MEAT2.0")
#'
BMIQcalibration <- function(SE,
                            version = "MEAT2.0") {
  gold.mean <- NULL
  nL <- 3
  doH <- TRUE
  nfit <- 20000
  th1.v <- c(0.2, 0.75)
  th2.v <- NULL
  niter <- 5
  tol <- 0.001
  calibrateUnitInterval <- TRUE

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

  gold.mean.MEAT <- NULL
  gold.mean.MEAT2.0 <- NULL
  if (version=="MEAT")
  {
    data("gold.mean.MEAT",envir = environment())
    gold.mean <- gold.mean.MEAT
  }
  else
  {
    data("gold.mean.MEAT2.0",envir = environment())
    gold.mean <- gold.mean.MEAT2.0
  }

  goldstandard.beta <- gold.mean[,"gold.mean"] #extract gold mean
  datM <- assays(SE)$beta

  datM <- t(datM)
  if (length(goldstandard.beta) != dim(datM)[[2]]) {
    stop(paste0("The number of probes of the the beta-matrix store in SE does not match
    the number of probes in the gold standard dataset (",length(goldstandard.beta),")."," Consider transposing the beta-matrix."))
  }
  beta1.v <- goldstandard.beta

  if (calibrateUnitInterval) {
    datM <- CalibrateUnitInterval(datM)
  }

  ### estimate initial weight matrix from type1 distribution
  w0.m <- matrix(0, nrow = length(beta1.v), ncol = nL)
  w0.m[which(beta1.v <= th1.v[1]), 1] <- 1
  w0.m[intersect(which(beta1.v > th1.v[1]), which(beta1.v <= th1.v[2])), 2] <- 1
  w0.m[which(beta1.v > th1.v[2]), 3] <- 1


  ### fit type1
  message("Fitting EM beta mixture to goldstandard probes.")
  rand.idx <- sample(seq_len(length(beta1.v)), min(c(nfit, length(beta1.v))), replace = FALSE)
  em1.o <- blc(matrix(beta1.v[rand.idx], ncol = 1), w = w0.m[rand.idx, ], maxiter = niter,
               tol = tol)
  subsetclass1.v <- apply(em1.o$w, 1, which.max)
  subsetth1.v <- c(mean(max(beta1.v[rand.idx[subsetclass1.v == 1]]), min(beta1.v[rand.idx[subsetclass1.v ==
                                                                                            2]])), mean(max(beta1.v[rand.idx[subsetclass1.v == 2]]), min(beta1.v[rand.idx[subsetclass1.v ==
                                                                                                                                                                            3]])))
  class1.v <- rep(2, length(beta1.v))
  class1.v[which(beta1.v < subsetth1.v[1])] <- 1
  class1.v[which(beta1.v > subsetth1.v[2])] <- 3
  nth1.v <- subsetth1.v
  message("Done")

  ### Estimate Modes
  if (sum(class1.v == 1) == 1) {
    mod1U <- beta1.v[class1.v == 1]
  }
  if (sum(class1.v == 3) == 1) {
    mod1M <- beta1.v[class1.v == 3]
  }
  if (sum(class1.v == 1) > 1) {
    d1U.o <- density(beta1.v[class1.v == 1])
    mod1U <- d1U.o$x[which.max(d1U.o$y)]
  }
  if (sum(class1.v == 3) > 1) {
    d1M.o <- density(beta1.v[class1.v == 3])
    mod1M <- d1M.o$x[which.max(d1M.o$y)]
  }

  ### BETA 2
  for (ii in seq_len(dim(datM)[[1]])) {
    printFlush(paste("ii=", ii))
    sampleID <- ii
    beta2.v <- as.numeric(datM[ii, ])

    d2U.o <- density(beta2.v[which(beta2.v < 0.4)])
    d2M.o <- density(beta2.v[which(beta2.v > 0.6)])
    mod2U <- d2U.o$x[which.max(d2U.o$y)]
    mod2M <- d2M.o$x[which.max(d2M.o$y)]

    ### now deal with type2 fit
    th2.v <- vector()
    th2.v[1] <- nth1.v[1] + (mod2U - mod1U)
    th2.v[2] <- nth1.v[2] + (mod2M - mod1M)

    ### estimate initial weight matrix
    w0.m <- matrix(0, nrow = length(beta2.v), ncol = nL)
    w0.m[which(beta2.v <= th2.v[1]), 1] <- 1
    w0.m[intersect(which(beta2.v > th2.v[1]), which(beta2.v <= th2.v[2])), 2] <- 1
    w0.m[which(beta2.v > th2.v[2]), 3] <- 1

    message("Fitting EM beta mixture to input probes")
    # I fixed an error in the following line (replaced beta1 by beta2)
    rand.idx <- sample(seq_len(length(beta2.v)), min(c(nfit, length(beta2.v)),
                                                     na.rm = TRUE), replace = FALSE)
    em2.o <- blc2(Y = matrix(beta2.v[rand.idx], ncol = 1), w = w0.m[rand.idx,
                                                                    ], maxiter = niter, tol = tol, verbose = TRUE)
    message("Done")

    ### for type II probes assign to state (un-, hemi- or full methylation)
    subsetclass2.v <- apply(em2.o$w, 1, which.max)


    if (sum(subsetclass2.v == 2) > 0) {
      subsetth2.v <- c(mean(max(beta2.v[rand.idx[subsetclass2.v == 1]]), min(beta2.v[rand.idx[subsetclass2.v ==
                                                                                                2]])), mean(max(beta2.v[rand.idx[subsetclass2.v == 2]]), min(beta2.v[rand.idx[subsetclass2.v ==
                                                                                                                                                                                3]])))
    }
    if (sum(subsetclass2.v == 2) == 0) {
      subsetth2.v <- c(1/2 * max(beta2.v[rand.idx[subsetclass2.v == 1]]) +
                         1/2 * mean(beta2.v[rand.idx[subsetclass2.v == 3]]), 1/3 * max(beta2.v[rand.idx[subsetclass2.v ==
                                                                                                          1]]) + 2/3 * mean(beta2.v[rand.idx[subsetclass2.v == 3]]))
    }



    class2.v <- rep(2, length(beta2.v))
    class2.v[which(beta2.v <= subsetth2.v[1])] <- 1
    class2.v[which(beta2.v >= subsetth2.v[2])] <- 3


    classAV1.v <- vector()
    classAV2.v <- vector()
    for (l in seq_len(nL)) {
      classAV1.v[l] <- em1.o$mu[l, 1]
      classAV2.v[l] <- em2.o$mu[l, 1]
    }

    ### start normalising input probes
    message("Start normalising input probes")
    nbeta2.v <- beta2.v
    ### select U probes
    lt <- 1
    selU.idx <- which(class2.v == lt)
    selUR.idx <- selU.idx[which(beta2.v[selU.idx] > classAV2.v[lt])]
    selUL.idx <- selU.idx[which(beta2.v[selU.idx] < classAV2.v[lt])]
    ### find prob according to typeII distribution
    p.v <- pbeta(beta2.v[selUR.idx], em2.o$a[lt, 1], em2.o$b[lt, 1], lower.tail = FALSE)
    ### find corresponding quantile in type I distribution
    q.v <- qbeta(p.v, em1.o$a[lt, 1], em1.o$b[lt, 1], lower.tail = FALSE)
    nbeta2.v[selUR.idx] <- q.v
    p.v <- pbeta(beta2.v[selUL.idx], em2.o$a[lt, 1], em2.o$b[lt, 1], lower.tail = TRUE)
    ### find corresponding quantile in type I distribution
    q.v <- qbeta(p.v, em1.o$a[lt, 1], em1.o$b[lt, 1], lower.tail = TRUE)
    nbeta2.v[selUL.idx] <- q.v

    ### select M probes
    lt <- 3
    selM.idx <- which(class2.v == lt)
    selMR.idx <- selM.idx[which(beta2.v[selM.idx] > classAV2.v[lt])]
    selML.idx <- selM.idx[which(beta2.v[selM.idx] < classAV2.v[lt])]
    ### find prob according to typeII distribution
    p.v <- pbeta(beta2.v[selMR.idx], em2.o$a[lt, 1], em2.o$b[lt, 1], lower.tail = FALSE)
    ### find corresponding quantile in type I distribution
    q.v <- qbeta(p.v, em1.o$a[lt, 1], em1.o$b[lt, 1], lower.tail = FALSE)
    nbeta2.v[selMR.idx] <- q.v


    if (doH) {
      ### if TRUE also correct type2 hemimethylated probes, select H probes and include
      ### ML probes (left ML tail is not well described by a beta-distribution).
      lt <- 2
      selH.idx <- c(which(class2.v == lt), selML.idx)
      minH <- min(beta2.v[selH.idx], na.rm = TRUE)
      maxH <- max(beta2.v[selH.idx], na.rm = TRUE)
      deltaH <- maxH - minH
      #### need to do some patching
      deltaUH <- -max(beta2.v[selU.idx], na.rm = TRUE) + min(beta2.v[selH.idx],
                                                             na.rm = TRUE)
      deltaHM <- -max(beta2.v[selH.idx], na.rm = TRUE) + min(beta2.v[selMR.idx],
                                                             na.rm = TRUE)

      ## new maximum of H probes should be
      nmaxH <- min(nbeta2.v[selMR.idx], na.rm = TRUE) - deltaHM
      ## new minimum of H probes should be
      nminH <- max(nbeta2.v[selU.idx], na.rm = TRUE) + deltaUH
      ndeltaH <- nmaxH - nminH

      ### perform conformal transformation (shift+dilation) new_beta_H(i) = a +
      ### hf*(beta_H(i)-minH);
      hf <- ndeltaH/deltaH
      ### fix lower point first
      nbeta2.v[selH.idx] <- nminH + hf * (beta2.v[selH.idx] - minH)
    }

    datM[ii, ] <- nbeta2.v
  }  # end of for (ii=1 loop
  t(datM)
SE2 <- SE
assays(SE2)$beta <- t(datM)
return(SE2)
}  # end of function BMIQcalibration




CalibrateUnitInterval <- function(datM, onlyIfOutside = TRUE) {
  rangeBySample <- data.frame(lapply(data.frame(t(datM)), range, na.rm = TRUE))
  minBySample <- as.numeric(rangeBySample[1, ])
  maxBySample <- as.numeric(rangeBySample[2, ])
  if (onlyIfOutside) {
    indexSamples <- which((minBySample < 0 | maxBySample > 1) & !is.na(minBySample) &
                            !is.na(maxBySample))
  }
  if (!onlyIfOutside) {
    indexSamples <- seq_len(length(minBySample))
  }
  if (length(indexSamples) >= 1) {
    for (i in indexSamples) {
      y1 <- c(0.001, 0.999)
      x1 <- c(minBySample[i], maxBySample[i])
      lm1 <- lm(y1 ~ x1)
      intercept1 <- coef(lm1)[[1]]
      slope1 <- coef(lm1)[[2]]
      datM[i, ] <- intercept1 + slope1 * datM[i, ]
    }  # end of for loop
  }
  datM
}
# end of function for calibrating to [0,1]


betaEst2 <- function(y, w, weights) {
  yobs <- !is.na(y)
  if (sum(yobs) <= 1) {
    return(c(1, 1))
  }
  y <- y[yobs]
  w <- w[yobs]
  weights <- weights[yobs]
  N <- sum(weights * w)
  p <- sum(weights * w * y)/N
  v <- sum(weights * w * y * y)/N - p * p
  logab <- log(c(p, 1 - p)) + log(pmax(1e-06, p * (1 - p)/v - 1))
  if (sum(yobs) == 2) {
    return(exp(logab))
  }
  opt <- try(optim(logab, betaObjf, ydata = y, wdata = w, weights = weights, method = "Nelder-Mead",
                   control = list(maxit = 50)), silent = TRUE)
  if (inherits(opt, "try-error")) {
    return(c(1, 1))
  }
  exp(opt$par)
}  # end of function betaEst



blc2 <- function(Y, w, maxiter = 25, tol = 1e-06, weights = NULL, verbose = TRUE) {
  Ymn <- min(Y[Y > 0], na.rm = TRUE)
  Ymx <- max(Y[Y < 1], na.rm = TRUE)
  Y <- pmax(Y, Ymn/2)
  Y <- pmin(Y, 1 - (1 - Ymx)/2)
  Yobs <- !is.na(Y)
  J <- dim(Y)[2]
  K <- dim(w)[2]
  n <- dim(w)[1]
  if (n != dim(Y)[1]) {
    stop("Dimensions of w and Y do not agree")
  }
  if (is.null(weights)) {
    weights <- rep(1, n)
  }
  mu <- a <- b <- matrix(Inf, K, J)
  crit <- Inf
  for (i in seq_len(maxiter)) {
    warn0 <- options()$warn
    options(warn = -1)
    eta <- apply(weights * w, 2, sum)/sum(weights)
    mu0 <- mu
    for (k in seq_len(K)) {
      for (j in seq_len(J)) {
        ab <- betaEst2(Y[, j], w[, k], weights)
        a[k, j] <- ab[1]
        b[k, j] <- ab[2]
        mu[k, j] <- ab[1]/sum(ab)
      }
    }
    ww <- array(0, dim = c(n, J, K))
    for (k in seq_len(K)) {
      for (j in seq_len(J)) {
        ww[Yobs[, j], j, k] <- dbeta(Y[Yobs[, j], j], a[k, j], b[k, j], log = TRUE)
      }
    }
    options(warn = warn0)
    w <- apply(ww, c(1, 3), sum, na.rm = TRUE)
    wmax <- apply(w, 1, max)
    for (k in seq_len(K)) w[, k] <- w[, k] - wmax
    w <- t(eta * t(exp(w)))
    like <- apply(w, 1, sum)
    w <- (1/like) * w
    llike <- weights * (log(like) + wmax)
    crit <- max(abs(mu - mu0))
    if (verbose) {
      message(crit)
    }
    if (crit < tol) {
      break
    }
  }
  return(list(a = a, b = b, eta = eta, mu = mu, w = w, llike = sum(llike)))
}
