### ORIGINAL AUTHOR: Andrew Teschendorff
# The original BMIQ function from Teschendorff 2013 adjusts for the type-2 bias in
# Illumina Infinium 450k data.
# Later functions and edits were provided by yours truly, Steve Horvath.
# I changed the code so that one can calibrate methylation data to a gold standard.
# Specifically, I took version v_1.2 by Teschendorff  and fixed minor issues.
# Also I made the code more robust e.g. by changing the optimization algorithm.
# Toward this end, I used the method="Nelder-Mead" in optim()
### Later functions and edits by Steve Horvath
### # Steve Horvath took version v_1.2 by Teschendorff
# and fixed minor errors. Also he made the code more robust.
# Importantly, SH changed the optimization algorithm to make it #more robust.
# SH used method="Nelder-Mead" in optim() since the other #optimization method sometimes gets stuck.
# Toward this end, the function blc was replaced by blc2.
# Requires complete data (no NAs) in datM and goldstandard.beta.
betaEst2 <- function(y, w, weights) {
  if (length(y) <= 1) {
    return(c(1, 1))
  }
  N <- sum(weights * w)
  p <- sum(weights * w * y) / N
  v <- sum(weights * w * y * y) / N - p * p
  logab <- log(c(p, 1 - p)) + log(pmax(1e-06, p * (1 - p) / v - 1))
  if (length(y) == 2) {
    return(exp(logab))
  }
  opt <- try(
    optim(
      logab,
      betaObjf,
      ydata = y,
      wdata = w,
      weights = weights,
      method = "Nelder-Mead",
      control = list(maxit = 50)
    ),
    silent = TRUE
  )
  if (inherits(opt, "try-error")) {
    return(c(1, 1))
  }
  exp(opt$par)
} # end of function betaEst
blc2 <- function(
  Y,
  w,
  maxiter = 25,
  tol = 1e-06,
  weights = NULL,
  verbose = TRUE
) {
  Ymn <- min(Y[Y > 0])
  Ymx <- max(Y[Y < 1])
  Y <- pmax(Y, Ymn / 2)
  Y <- pmin(Y, 1 - (1 - Ymx) / 2)
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
  for (i in 1:maxiter) {
    warn0 <- options()$warn
    options(warn = -1)
    eta <- apply(weights * w, 2, sum) / sum(weights)
    mu0 <- mu
    for (k in 1:K) {
      for (j in 1:J) {
        ab <- betaEst2(Y[, j], w[, k], weights)
        a[k, j] <- ab[1]
        b[k, j] <- ab[2]
        mu[k, j] <- ab[1] / sum(ab)
      }
    }
    ww <- array(0, dim = c(n, J, K))
    for (k in 1:K) {
      for (j in 1:J) {
        ww[, j, k] <- dbeta(Y[, j], a[k, j], b[k, j], log = TRUE)
      }
    }
    options(warn = warn0)
    w <- apply(ww, c(1, 3), sum)
    wmax <- apply(w, 1, max)
    for (k in 1:K) {
      w[, k] <- w[, k] - wmax
    }
    w <- t(eta * t(exp(w)))
    like <- apply(w, 1, sum)
    w <- (1 / like) * w
    llike <- weights * (log(like) + wmax)
    crit <- max(abs(mu - mu0))
    if (verbose) {
      print(crit)
    }
    if (crit < tol) {
      break
    }
  }
  return(list(a = a, b = b, eta = eta, mu = mu, w = w, llike = sum(llike)))
}

# The function BMIQcalibration was created by Steve Horvath by heavily recycling code
# from A. Teschendorff's BMIQ function.
# BMIQ stands for beta mixture quantile normalization.
# Explanation: datM is a data frame with Illumina beta values (rows are samples, colums are CpGs.
# goldstandard is a numeric vector with beta values that is used as gold standard for calibrating the columns of datM.
# The length of goldstandard has to equal the number of columns of datM.
# Both datM and goldstandard.beta must be finite and in [0, 1] (no NA/NaN/Inf, no values outside).
# Example:
# gold.mean=as.numeric(apply(datMeth,2,mean))
# datMethCalibrated=BMIQcalibration(datM=datMeth,goldstandard.beta=gold.mean)
BMIQcalibration <- function(
  datM,
  goldstandard.beta,
  nL = 3,
  doH = TRUE,
  nfit = 20000,
  th1.v = c(0.2, 0.75),
  th2.v = NULL,
  niter = 5,
  tol = 0.001
) {
  if (length(goldstandard.beta) != dim(datM)[[2]]) {
    stop(
      "Error in function arguments length(goldstandard.beta) !=dim(datM)[[2]]. Consider transposing datM."
    )
  }
  assertFiniteUnitInterval(datM, "datM")
  assertFiniteUnitInterval(goldstandard.beta, "goldstandard.beta")
  beta1.v <- goldstandard.beta
  ### estimate initial weight matrix from type1 distribution
  w0.m <- matrix(0, nrow = length(beta1.v), ncol = nL)
  w0.m[which(beta1.v <= th1.v[1]), 1] <- 1
  w0.m[intersect(which(beta1.v > th1.v[1]), which(beta1.v <= th1.v[2])), 2] <- 1
  w0.m[which(beta1.v > th1.v[2]), 3] <- 1
  ### fit type1
  print("Fitting EM beta mixture to goldstandard probes")
  set.seed(1)
  rand.idx <- sample(
    1:length(beta1.v),
    min(c(nfit, length(beta1.v))),
    replace = FALSE
  )
  em1.o <- blc(
    matrix(beta1.v[rand.idx], ncol = 1),
    w = w0.m[rand.idx, ],
    maxiter = niter,
    tol = tol
  )
  subsetclass1.v <- apply(em1.o$w, 1, which.max)
  subsetth1.v <- c(
    mean(
      max(beta1.v[rand.idx[subsetclass1.v == 1]]),
      min(beta1.v[rand.idx[subsetclass1.v == 2]])
    ),
    mean(
      max(beta1.v[rand.idx[subsetclass1.v == 2]]),
      min(beta1.v[rand.idx[subsetclass1.v == 3]])
    )
  )
  class1.v <- rep(2, length(beta1.v))
  class1.v[which(beta1.v < subsetth1.v[1])] <- 1
  class1.v[which(beta1.v > subsetth1.v[2])] <- 3
  nth1.v <- subsetth1.v
  print("Done")
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
  for (ii in 1:dim(datM)[[1]]) {
    print(paste("ii=", ii, "of", dim(datM)[[1]]))
    print(date())
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
    w0.m[
      intersect(which(beta2.v > th2.v[1]), which(beta2.v <= th2.v[2])),
      2
    ] <- 1
    w0.m[which(beta2.v > th2.v[2]), 3] <- 1
    print("Fitting EM beta mixture to input probes")
    # I fixed an error in the following line (replaced beta1 by beta2)
    set.seed(1)
    rand.idx <- sample(
      1:length(beta2.v),
      min(c(nfit, length(beta2.v))),
      replace = FALSE
    )
    em2.o <- blc2(
      Y = matrix(beta2.v[rand.idx], ncol = 1),
      w = w0.m[rand.idx, ],
      maxiter = niter,
      tol = tol,
      verbose = TRUE
    )
    print("Done")
    ### for type II probes assign to state (unmethylated, hemi or full methylation)
    subsetclass2.v <- apply(em2.o$w, 1, which.max)
    if (sum(subsetclass2.v == 2) > 0) {
      subsetth2.v <- c(
        mean(
          max(beta2.v[rand.idx[subsetclass2.v == 1]]),
          min(beta2.v[rand.idx[subsetclass2.v == 2]])
        ),
        mean(
          max(beta2.v[rand.idx[subsetclass2.v == 2]]),
          min(beta2.v[rand.idx[subsetclass2.v == 3]])
        )
      )
    }
    if (sum(subsetclass2.v == 2) == 0) {
      subsetth2.v <- c(
        1 /
          2 *
          max(beta2.v[rand.idx[subsetclass2.v == 1]]) +
          1 / 2 * mean(beta2.v[rand.idx[subsetclass2.v == 3]]),
        1 /
          3 *
          max(beta2.v[rand.idx[subsetclass2.v == 1]]) +
          2 / 3 * mean(beta2.v[rand.idx[subsetclass2.v == 3]])
      )
    }
    class2.v <- rep(2, length(beta2.v))
    class2.v[which(beta2.v <= subsetth2.v[1])] <- 1
    class2.v[which(beta2.v >= subsetth2.v[2])] <- 3
    classAV1.v <- vector()
    classAV2.v <- vector()
    for (l in 1:nL) {
      classAV1.v[l] <- em1.o$mu[l, 1]
      classAV2.v[l] <- em2.o$mu[l, 1]
    }
    ### start normalising input probes
    print("Start normalising input probes")
    nbeta2.v <- beta2.v
    ### select U probes
    lt <- 1
    selU.idx <- which(class2.v == lt)
    selUR.idx <- selU.idx[which(beta2.v[selU.idx] > classAV2.v[lt])]
    selUL.idx <- selU.idx[which(beta2.v[selU.idx] < classAV2.v[lt])]
    ### find prob according to typeII distribution
    p.v <- pbeta(
      beta2.v[selUR.idx],
      em2.o$a[lt, 1],
      em2.o$b[lt, 1],
      lower.tail = FALSE
    )
    ### find corresponding quantile in type I distribution
    q.v <- qbeta(p.v, em1.o$a[lt, 1], em1.o$b[lt, 1], lower.tail = FALSE)
    nbeta2.v[selUR.idx] <- q.v
    p.v <- pbeta(
      beta2.v[selUL.idx],
      em2.o$a[lt, 1],
      em2.o$b[lt, 1],
      lower.tail = TRUE
    )
    ### find corresponding quantile in type I distribution
    q.v <- qbeta(p.v, em1.o$a[lt, 1], em1.o$b[lt, 1], lower.tail = TRUE)
    nbeta2.v[selUL.idx] <- q.v
    ### select M probes
    lt <- 3
    selM.idx <- which(class2.v == lt)
    selMR.idx <- selM.idx[which(beta2.v[selM.idx] > classAV2.v[lt])]
    selML.idx <- selM.idx[which(beta2.v[selM.idx] < classAV2.v[lt])]
    ### find prob according to typeII distribution
    p.v <- pbeta(
      beta2.v[selMR.idx],
      em2.o$a[lt, 1],
      em2.o$b[lt, 1],
      lower.tail = FALSE
    )
    ### find corresponding quantile in type I distribution
    q.v <- qbeta(p.v, em1.o$a[lt, 1], em1.o$b[lt, 1], lower.tail = FALSE)
    nbeta2.v[selMR.idx] <- q.v
    if (doH) {
      ### if TRUE also correct type2 hemimethylated probes
      ### select H probes and include ML probes (left ML tail is not well described by a beta-distribution).
      lt <- 2
      selH.idx <- c(which(class2.v == lt), selML.idx)
      minH <- min(beta2.v[selH.idx])
      maxH <- max(beta2.v[selH.idx])
      deltaH <- maxH - minH
      #### need to do some patching
      deltaUH <- -max(beta2.v[selU.idx]) + min(beta2.v[selH.idx])
      deltaHM <- -max(beta2.v[selH.idx]) + min(beta2.v[selMR.idx])
      ## new maximum of H probes should be
      nmaxH <- min(nbeta2.v[selMR.idx]) - deltaHM
      ## new minimum of H probes should be
      nminH <- max(nbeta2.v[selU.idx]) + deltaUH
      ndeltaH <- nmaxH - nminH
      ### perform conformal transformation (shift+dilation)
      ## new_beta_H(i) = a + hf*(beta_H(i)-minH);
      hf <- ndeltaH / deltaH
      ### fix lower point first
      nbeta2.v[selH.idx] <- nminH + hf * (beta2.v[selH.idx] - minH)
    }
    datM[ii, ] <- nbeta2.v
  } # end of for (ii=1 loop
  datM
} # end of function BMIQcalibration

# Vendored from RPMM (Houseman et al.); requires complete data (no NAs).
betaObjf <- function(logab, ydata, wdata, weights) {
  ab <- exp(logab)
  -sum(wdata * weights * dbeta(ydata, ab[1], ab[2], log = TRUE))
}

# Vendored from RPMM betaEst (used by blc); BFGS optim, complete data only.
betaEst <- function(y, w, weights) {
  if (length(y) <= 1) {
    return(c(1, 1))
  }
  N <- sum(weights * w)
  p <- sum(weights * w * y) / N
  v <- sum(weights * w * y * y) / N - p * p
  logab <- log(c(p, 1 - p)) + log(pmax(1e-06, p * (1 - p) / v - 1))
  if (length(y) == 2) {
    return(exp(logab))
  }
  opt <- try(
    optim(logab, betaObjf, ydata = y, wdata = w, weights = weights, method = "BFGS"),
    silent = TRUE
  )
  if (inherits(opt, "try-error")) {
    return(c(1, 1))
  }
  exp(opt$par)
}

# Vendored from RPMM blc; complete data only.
blc <- function(
  Y,
  w,
  maxiter = 25,
  tol = 1e-06,
  weights = NULL,
  verbose = TRUE
) {
  Ymn <- min(Y[Y > 0])
  Ymx <- max(Y[Y < 1])
  Y <- pmax(Y, Ymn / 2)
  Y <- pmin(Y, 1 - (1 - Ymx) / 2)
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
  for (i in 1:maxiter) {
    warn0 <- options()$warn
    options(warn = -1)
    eta <- apply(weights * w, 2, sum) / sum(weights)
    mu0 <- mu
    for (k in 1:K) {
      for (j in 1:J) {
        ab <- betaEst(Y[, j], w[, k], weights)
        a[k, j] <- ab[1]
        b[k, j] <- ab[2]
        mu[k, j] <- ab[1] / sum(ab)
      }
    }
    ww <- array(0, dim = c(n, J, K))
    for (k in 1:K) {
      for (j in 1:J) {
        ww[, j, k] <- dbeta(Y[, j], a[k, j], b[k, j], log = TRUE)
      }
    }
    options(warn = warn0)
    w <- apply(ww, c(1, 3), sum)
    wmax <- apply(w, 1, max)
    for (k in 1:K) {
      w[, k] <- w[, k] - wmax
    }
    w <- t(eta * t(exp(w)))
    like <- apply(w, 1, sum)
    w <- (1 / like) * w
    llike <- weights * (log(like) + wmax)
    crit <- max(abs(mu - mu0))
    if (verbose) {
      print(crit)
    }
    if (crit < tol) {
      break
    }
  }
  return(list(a = a, b = b, eta = eta, mu = mu, w = w, llike = sum(llike)))
}

# Abort unless all values are finite and in [0, 1] (covers NA, NaN, +/-Inf, out-of-range).
assertFiniteUnitInterval <- function(x, name) {
  if (any(!is.finite(x))) {
    stop(
      name,
      " must contain only finite values (no NA, NaN, or +/-Inf).",
      call. = FALSE
    )
  }
  if (any(x < 0 | x > 1)) {
    stop(
      name,
      " must have all values in [0, 1]; values outside that range are not allowed.",
      call. = FALSE
    )
  }
  invisible(NULL)
}
