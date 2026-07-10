# Below is a complete revision implementing the agreed policies:

# - real midpoint thresholds, intentionally correcting the legacy `mean()` trap;
# - one consistent equality rule: equality stays in the lower class;
# - hard failures for empty or unusable classes and normalization anchors;
# - strictly ordered fitted component means;
# - sample-level continuation through `on.sample.error = "continue"`;
# - failed rows become `NA` by default;
# - gold-standard failures always abort;
# - detailed diagnostics are collected only when `debug = TRUE`;
# - `nL > 3` remains supported but explicitly experimental;
# - no OpenMP;
# - no legacy invented thresholds or silent `Beta(1,1)` fallback;
# - structured result returned consistently.

# -----------------------------------------------------------------------------
# Clip raw beta values from [0, 1] into (0, 1) for Beta-distribution fitting.
#
# This follows the endpoint treatment in historical blc()/blc2(). The complete
# raw input is validated separately before fitting.
# -----------------------------------------------------------------------------
clipBetaForFit <- function(y, name = "y") {
  positive <- y[y > 0]
  below.one <- y[y < 1]

  if (!length(positive)) {
    stop(
      name,
      " contains no values greater than zero; ",
      "a Beta mixture cannot be fitted.",
      call. = FALSE
    )
  }

  if (!length(below.one)) {
    stop(
      name,
      " contains no values below one; ",
      "a Beta mixture cannot be fitted.",
      call. = FALSE
    )
  }

  ymin <- min(positive)
  ymax <- max(below.one)

  y <- pmax(y, ymin / 2)
  y <- pmin(y, 1 - (1 - ymax) / 2)

  if (any(!is.finite(y)) || any(y <= 0 | y >= 1)) {
    stop(
      name,
      " could not be clipped safely into (0, 1).",
      call. = FALSE
    )
  }

  y
}


# -----------------------------------------------------------------------------
# Validate nL - 1 ordered thresholds.
# -----------------------------------------------------------------------------
validateThresholds <- function(
  thresholds,
  nL,
  name,
  require.unit.interval = FALSE
) {
  if (length(thresholds) != nL - 1L) {
    stop(
      name,
      " must have length nL - 1 (expected ",
      nL - 1L,
      ", got ",
      length(thresholds),
      ").",
      call. = FALSE
    )
  }

  if (any(!is.finite(thresholds))) {
    stop(
      name,
      " must contain only finite values.",
      call. = FALSE
    )
  }

  if (is.unsorted(thresholds, strictly = TRUE)) {
    stop(
      name,
      " must be strictly increasing.",
      call. = FALSE
    )
  }

  if (
    require.unit.interval &&
      any(thresholds <= 0 | thresholds >= 1)
  ) {
    stop(
      name,
      " must lie strictly inside (0, 1).",
      call. = FALSE
    )
  }

  invisible(NULL)
}


# -----------------------------------------------------------------------------
# Assign values to ordered classes.
#
# Equality convention:
#   A value equal to a threshold stays in the lower class.
#
# For nL = 3:
#   class 1: x <= threshold[1]
#   class 2: threshold[1] < x <= threshold[2]
#   class 3: x > threshold[2]
#
# This intentionally removes inconsistent legacy equality conventions.
# -----------------------------------------------------------------------------
classifyByThresholds <- function(beta, thresholds, nL) {
  validateThresholds(
    thresholds,
    nL = nL,
    name = "thresholds",
    require.unit.interval = FALSE
  )

  class <- rep.int(1L, length(beta))

  for (boundary in seq_along(thresholds)) {
    class[beta > thresholds[boundary]] <- boundary + 1L
  }

  class
}


# -----------------------------------------------------------------------------
# Construct hard initial responsibilities.
# -----------------------------------------------------------------------------
initialResponsibilities <- function(beta, thresholds, nL) {
  class <- classifyByThresholds(beta, thresholds, nL)

  responsibility <- matrix(
    0,
    nrow = length(beta),
    ncol = nL
  )

  responsibility[
    cbind(seq_along(beta), class)
  ] <- 1

  responsibility
}


# -----------------------------------------------------------------------------
# Require all expected classes to have at least min.count observations.
#
# The finite-unit-interval scanner cannot guarantee this. A valid beta vector
# may still have an empty U, H, M, mode window, MAP class, or normalization
# anchor.
# -----------------------------------------------------------------------------
requireAllClasses <- function(
  class,
  nL,
  context,
  min.count = 1L
) {
  counts <- tabulate(class, nbins = nL)
  missing <- which(counts < min.count)

  if (length(missing)) {
    stop(
      context,
      " does not contain enough observations in class(es): ",
      paste(missing, collapse = ", "),
      ". Counts: ",
      paste(counts, collapse = ", "),
      ". Required count per class: ",
      min.count,
      ".",
      call. = FALSE
    )
  }

  counts
}


# -----------------------------------------------------------------------------
# Compute boundaries between adjacent MAP classes.
#
# Intentional legacy correction:
#
# Legacy code used:
#
#   mean(max(lower), min(upper))
#
# In R the second argument is matched to mean.default()'s `trim` argument.
# It therefore normally returned max(lower), not the midpoint.
#
# This implementation uses the intended midpoint:
#
#   mean(c(max(lower), min(upper)))
#
# This can change thresholds and class labels even when the EM fit is the same.
# -----------------------------------------------------------------------------
thresholdsFromPosteriorClasses <- function(
  beta,
  class,
  nL,
  context
) {
  requireAllClasses(
    class,
    nL = nL,
    context = context,
    min.count = 1L
  )

  thresholds <- numeric(nL - 1L)

  for (boundary in seq_len(nL - 1L)) {
    lower.values <- beta[class == boundary]
    upper.values <- beta[class == boundary + 1L]

    thresholds[boundary] <- mean(
      c(
        max(lower.values),
        min(upper.values)
      )
    )
  }

  if (any(!is.finite(thresholds))) {
    stop(
      context,
      " produced non-finite thresholds.",
      call. = FALSE
    )
  }

  if (is.unsorted(thresholds, strictly = TRUE)) {
    stop(
      context,
      " produced thresholds that are not strictly increasing: ",
      paste(signif(thresholds, 8), collapse = ", "),
      ". MAP component labels may not follow beta-value order.",
      call. = FALSE
    )
  }

  thresholds
}


# -----------------------------------------------------------------------------
# Estimate a mode using density(), with explicit treatment of singleton and
# constant regions.
# -----------------------------------------------------------------------------
estimateMode <- function(x, context) {
  if (!length(x)) {
    stop(context, " is empty.", call. = FALSE)
  }

  if (any(!is.finite(x))) {
    stop(
      context,
      " contains non-finite values.",
      call. = FALSE
    )
  }

  if (length(x) == 1L || all(x == x[1L])) {
    return(x[1L])
  }

  estimate <- density(x)

  if (
    !length(estimate$x) ||
      !length(estimate$y) ||
      any(!is.finite(estimate$x)) ||
      any(!is.finite(estimate$y))
  ) {
    stop(
      "Density estimation failed for ",
      context,
      ".",
      call. = FALSE
    )
  }

  estimate$x[which.max(estimate$y)]
}


# -----------------------------------------------------------------------------
# Ensure that fitted component means follow U -> H/interior -> M order.
#
# EM labels are initialized in this order but are not mathematically constrained
# to remain ordered. A nonempty MAP component is not sufficient if, for example,
# the nominal H component migrates to a second high-beta peak.
#
# Components are not silently reordered because responsibilities and
# normalization semantics depend on their original identities.
# -----------------------------------------------------------------------------
requireOrderedComponentMeans <- function(mu, context) {
  mu <- as.numeric(mu)

  if (!length(mu) || any(!is.finite(mu))) {
    stop(
      context,
      " has missing or non-finite component means.",
      call. = FALSE
    )
  }

  if (is.unsorted(mu, strictly = TRUE)) {
    stop(
      context,
      " component means are not strictly increasing: ",
      paste(signif(mu, 8), collapse = ", "),
      ". The fitted components cannot safely be interpreted as ",
      "ordered U/H/M states.",
      call. = FALSE
    )
  }

  invisible(mu)
}


# -----------------------------------------------------------------------------
# Apply the requested convergence policy.
#
# "usable":
#   Accept finite converged, max-iteration, or stalled component estimates.
#   Unusable components already cause a hard error in C++.
#
# "converged":
#   Require both the outer EM and every most-recent component optimization to
#   report convergence. This is intentionally strict and may reject otherwise
#   usable fits.
# -----------------------------------------------------------------------------
enforceFitPolicy <- function(em, fit.policy, context) {
  if (fit.policy == "usable") {
    return(invisible(NULL))
  }

  if (
    !isTRUE(em$converged) ||
      any(em$fit_status != "converged")
  ) {
    stop(
      context,
      " did not fully converge under fit.policy = \"converged\". ",
      "EM converged: ",
      isTRUE(em$converged),
      "; component statuses: ",
      paste(em$fit_status, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  invisible(NULL)
}


# -----------------------------------------------------------------------------
# Fit one Beta mixture through C++.
# -----------------------------------------------------------------------------
fitBetaMixture <- function(
  beta,
  initial.responsibility,
  nL = 3L,
  maxiter = 5L,
  tol = 0.001,
  beta.maxit = 50L,
  beta.score.tol = 1e-10,
  debug = FALSE
) {
  fit.beta <- clipBetaForFit(beta, "EM beta vector")

  beta_mixture_em_cpp(
    y = fit.beta,
    initial_responsibility = initial.responsibility,
    nL = nL,
    maxiter = maxiter,
    tol = tol,
    beta_maxit = beta.maxit,
    beta_score_tol = beta.score.tol,
    debug = debug
  )
}


# -----------------------------------------------------------------------------
# Internal sample-error condition.
#
# Sample-specific errors can be caught when on.sample.error = "continue".
# Gold-standard and global input errors are never converted to this condition.
# -----------------------------------------------------------------------------
newBMIQSampleError <- function(
  sample.index,
  sample.name,
  stage,
  error,
  diagnostics = NULL
) {
  original.message <- conditionMessage(error)

  structure(
    list(
      message = paste0(
        "Sample ",
        sample.index,
        if (nzchar(sample.name)) {
          paste0(" (", sample.name, ")")
        } else {
          ""
        },
        " failed during ",
        stage,
        ": ",
        original.message
      ),
      call = NULL,
      sample.index = sample.index,
      sample.name = sample.name,
      stage = stage,
      original.message = original.message,
      diagnostics = diagnostics,
      parent = error
    ),
    class = c(
      "bmiq_sample_error",
      "error",
      "condition"
    )
  )
}


emptyBMIQFailureTable <- function() {
  data.frame(
    sample_index = integer(),
    sample_name = character(),
    stage = character(),
    message = character(),
    stringsAsFactors = FALSE
  )
}


#' Calibrate methylation beta values against a gold standard
#'
#' This implementation uses C++ Beta-mixture EM fitting and deliberately
#' corrects several legacy behaviors:
#'
#' * `mean(max(lower), min(upper))` is replaced by the intended midpoint.
#' * Threshold equality consistently remains in the lower class.
#' * Empty MAP classes and unusable normalization anchors are hard failures.
#' * No invented H thresholds or silent Beta(1,1) fallback is used.
#' * Fitted component means must remain strictly ordered.
#' * Detailed diagnostics are retained only when `debug = TRUE`.
#'
#' Gold-standard failures always abort the complete operation. Sample failures
#' either abort or are recorded and skipped according to `on.sample.error`.
#'
#' `nL = 3` is the standard BMIQ U/H/M interpretation. For `nL > 3`, every
#' interior component is currently treated collectively as intermediate/H.
#' This is experimental and may not be biologically appropriate when additional
#' components represent multiple high-methylation peaks.
#'
#' @param datM Numeric matrix with samples in rows and CpGs in columns.
#' @param goldstandard.beta Numeric vector with one value per column of datM.
#' @param nL Number of ordered mixture components. Default is 3.
#' @param doH Whether to perform intermediate/H conformal normalization.
#' @param nfit Maximum number of probes used in each EM fit.
#' @param th1.v Initial gold-standard boundaries; must have length nL - 1.
#' @param niter Maximum number of outer EM iterations.
#' @param tol Outer EM component-mean convergence tolerance.
#' @param beta.maxit Maximum Newton iterations per Beta-component fit.
#' @param beta.score.tol Newton score convergence tolerance.
#' @param fit.policy `"usable"` accepts finite stalled/max-iteration fits;
#'   `"converged"` requires complete component and outer-EM convergence.
#' @param on.sample.error `"stop"` aborts on the first sample failure;
#'   `"continue"` records the failure and processes later samples.
#' @param failed.sample For continued failures, use an `"NA"` row or retain the
#'   `"original"` uncalibrated row. `"NA"` is safer and is the default.
#' @param debug If TRUE, retain detailed gold and per-sample diagnostics.
#' @param verbose Print progress messages.
#'
#' @return A `BMIQcalibration_result` list containing:
#'   `calibrated`, `success`, `failures`, and optionally `diagnostics`.
BMIQcalibration <- function(
  datM,
  goldstandard.beta,
  nL = 3L,
  doH = TRUE,
  nfit = 20000L,
  th1.v = c(0.2, 0.75),
  niter = 5L,
  tol = 0.001,
  beta.maxit = 50L,
  beta.score.tol = 1e-10,
  fit.policy = c("usable", "converged"),
  on.sample.error = c("stop", "continue"),
  failed.sample = c("NA", "original"),
  debug = FALSE,
  verbose = TRUE
) {
  call <- match.call()

  fit.policy <- match.arg(fit.policy)
  on.sample.error <- match.arg(on.sample.error)
  failed.sample <- match.arg(failed.sample)

  if (
    length(debug) != 1L ||
      is.na(debug) ||
      !is.logical(debug)
  ) {
    stop("debug must be TRUE or FALSE.", call. = FALSE)
  }

  if (
    length(verbose) != 1L ||
      is.na(verbose) ||
      !is.logical(verbose)
  ) {
    stop("verbose must be TRUE or FALSE.", call. = FALSE)
  }

  if (
    length(doH) != 1L ||
      is.na(doH) ||
      !is.logical(doH)
  ) {
    stop("doH must be TRUE or FALSE.", call. = FALSE)
  }

  if (is.data.frame(datM)) {
    datM <- data.matrix(datM)
  }

  if (!is.matrix(datM) || !is.numeric(datM)) {
    stop(
      "datM must be a numeric matrix or data frame.",
      call. = FALSE
    )
  }

  storage.mode(datM) <- "double"

  if (!nrow(datM) || !ncol(datM)) {
    stop("datM must have at least one row and column.", call. = FALSE)
  }

  if (!is.numeric(goldstandard.beta)) {
    stop(
      "goldstandard.beta must be numeric.",
      call. = FALSE
    )
  }

  goldstandard.beta <- as.numeric(goldstandard.beta)

  nL <- as.integer(nL)
  nfit <- as.integer(nfit)
  niter <- as.integer(niter)
  beta.maxit <- as.integer(beta.maxit)

  if (length(nL) != 1L || is.na(nL) || nL < 3L) {
    stop(
      "nL must be a single integer of at least 3.",
      call. = FALSE
    )
  }

  if (nL > 3L) {
    warning(
      "nL > 3 is experimental. All components between the first ",
      "and last are treated as intermediate/H components. This may ",
      "not be appropriate if extra components represent multiple ",
      "high-methylation peaks.",
      call. = FALSE
    )
  }

  if (
    length(nfit) != 1L ||
      is.na(nfit) ||
      nfit < 2L * nL
  ) {
    stop(
      "nfit must be at least 2 * nL.",
      call. = FALSE
    )
  }

  if (
    length(niter) != 1L ||
      is.na(niter) ||
      niter < 1L
  ) {
    stop(
      "niter must be a positive integer.",
      call. = FALSE
    )
  }

  if (
    length(beta.maxit) != 1L ||
      is.na(beta.maxit) ||
      beta.maxit < 1L
  ) {
    stop(
      "beta.maxit must be a positive integer.",
      call. = FALSE
    )
  }

  if (
    length(tol) != 1L ||
      !is.finite(tol) ||
      tol <= 0
  ) {
    stop(
      "tol must be a positive finite number.",
      call. = FALSE
    )
  }

  if (
    length(beta.score.tol) != 1L ||
      !is.finite(beta.score.tol) ||
      beta.score.tol <= 0
  ) {
    stop(
      "beta.score.tol must be a positive finite number.",
      call. = FALSE
    )
  }

  if (length(goldstandard.beta) != ncol(datM)) {
    stop(
      "length(goldstandard.beta) must equal ncol(datM). ",
      "Consider transposing datM.",
      call. = FALSE
    )
  }

  validateThresholds(
    th1.v,
    nL = nL,
    name = "th1.v",
    require.unit.interval = TRUE
  )

  scan_finite_unit_interval_cpp(
    datM,
    name = "datM",
    require_open = FALSE
  )

  # The C++ scanner accepts arma::mat, so the vector is presented as a
  # one-column matrix.
  scan_finite_unit_interval_cpp(
    matrix(goldstandard.beta, ncol = 1L),
    name = "goldstandard.beta",
    require_open = FALSE
  )

  number.of.samples <- nrow(datM)
  number.of.probes <- ncol(datM)

  if (number.of.probes < 2L * nL) {
    stop(
      "There are too few probes to fit ",
      nL,
      " mixture components.",
      call. = FALSE
    )
  }

  sample.names <- rownames(datM)

  if (is.null(sample.names)) {
    sample.names <- rep.int("", number.of.samples)
  }

  # Preserve the original matrix so failed rows can explicitly remain
  # uncalibrated when failed.sample = "original".
  original.datM <- datM
  calibrated <- datM

  success <- rep.int(FALSE, number.of.samples)
  failures <- list()

  sample.diagnostics <- if (debug) {
    vector("list", number.of.samples)
  } else {
    NULL
  }

  # ===========================================================================
  # Gold-standard fit
  #
  # Gold-standard failures are global failures and are never skipped.
  # ===========================================================================
  beta1.v <- goldstandard.beta

  gold.w0 <- initialResponsibilities(
    beta = beta1.v,
    thresholds = th1.v,
    nL = nL
  )

  if (verbose) {
    message("Fitting EM beta mixture to gold-standard probes")
  }

  set.seed(1)

  gold.rand.idx <- sample.int(
    number.of.probes,
    min(nfit, number.of.probes),
    replace = FALSE
  )

  gold.initial.class <- max.col(
    gold.w0[gold.rand.idx, , drop = FALSE],
    ties.method = "first"
  )

  gold.initial.counts <- requireAllClasses(
    class = gold.initial.class,
    nL = nL,
    context = "Gold-standard initial mixture",
    min.count = 2L
  )

  em1.o <- fitBetaMixture(
    beta = beta1.v[gold.rand.idx],
    initial.responsibility =
      gold.w0[gold.rand.idx, , drop = FALSE],
    nL = nL,
    maxiter = niter,
    tol = tol,
    beta.maxit = beta.maxit,
    beta.score.tol = beta.score.tol,
    debug = debug
  )

  enforceFitPolicy(
    em1.o,
    fit.policy = fit.policy,
    context = "Gold-standard mixture"
  )

  classAV1.v <- as.numeric(em1.o$mu[, 1L])

  requireOrderedComponentMeans(
    classAV1.v,
    "Gold-standard mixture"
  )

  gold.subset.class <- max.col(
    em1.o$w,
    ties.method = "first"
  )

  gold.subset.counts <- requireAllClasses(
    class = gold.subset.class,
    nL = nL,
    context = "Gold-standard posterior mixture",
    min.count = 1L
  )

  nth1.v <- thresholdsFromPosteriorClasses(
    beta = beta1.v[gold.rand.idx],
    class = gold.subset.class,
    nL = nL,
    context = "Gold-standard posterior mixture"
  )

  gold.full.class <- classifyByThresholds(
    beta = beta1.v,
    thresholds = nth1.v,
    nL = nL
  )

  gold.full.counts <- requireAllClasses(
    class = gold.full.class,
    nL = nL,
    context = "Complete gold-standard mixture",
    min.count = 1L
  )

  mod1U <- estimateMode(
    beta1.v[gold.full.class == 1L],
    "Gold-standard unmethylated class"
  )

  mod1M <- estimateMode(
    beta1.v[gold.full.class == nL],
    "Gold-standard methylated class"
  )

  gold.diagnostics <- if (debug) {
    list(
      random_indices = gold.rand.idx,
      initial_class_counts = gold.initial.counts,
      eta = em1.o$eta,
      component_means = classAV1.v,
      component_a = as.numeric(em1.o$a[, 1L]),
      component_b = as.numeric(em1.o$b[, 1L]),
      component_fit_status = em1.o$fit_status,
      component_fit_reason = em1.o$fit_reason,
      em_iterations = em1.o$iterations,
      em_converged = em1.o$converged,
      em_criterion = em1.o$criterion,
      log_likelihood = em1.o$llike,
      subset_map_counts = gold.subset.counts,
      thresholds = nth1.v,
      complete_class_counts = gold.full.counts,
      unmethylated_mode = mod1U,
      methylated_mode = mod1M
    )
  } else {
    NULL
  }

  if (verbose) {
    message("Gold-standard mixture fit complete")
  }

  # ===========================================================================
  # One-sample worker
  #
  # All changes are made to a local beta vector. The output matrix is updated
  # only after the complete sample succeeds.
  # ===========================================================================
  processOneSample <- function(ii) {
    beta2.v <- as.numeric(original.datM[ii, ])
    sample.name <- sample.names[ii]
    stage <- "initialization"

    diagnostic <- if (debug) {
      list(
        sample_index = ii,
        sample_name = sample.name,
        input_range = range(beta2.v)
      )
    } else {
      NULL
    }

    tryCatch(
      {
        # ---------------------------------------------------------------------
        # Sample modes
        # ---------------------------------------------------------------------
        stage <- "sample mode estimation"

        low.mode.values <- beta2.v[beta2.v < 0.4]
        high.mode.values <- beta2.v[beta2.v > 0.6]

        mod2U <- estimateMode(
          low.mode.values,
          paste0("Sample ", ii, " values below 0.4")
        )

        mod2M <- estimateMode(
          high.mode.values,
          paste0("Sample ", ii, " values above 0.6")
        )

        if (debug) {
          diagnostic$low_mode_window_count <-
            length(low.mode.values)
          diagnostic$high_mode_window_count <-
            length(high.mode.values)
          diagnostic$unmethylated_mode <- mod2U
          diagnostic$methylated_mode <- mod2M
        }

        # ---------------------------------------------------------------------
        # Initial type-2 thresholds
        # ---------------------------------------------------------------------
        stage <- "initial threshold construction"

        unmethylated.shift <- mod2U - mod1U
        methylated.shift <- mod2M - mod1M

        # For nL = 3, alpha = c(0, 1), reproducing the two historical
        # end-mode shifts. For nL > 3, shifts are linearly interpolated.
        # That nL > 3 extension is computational rather than a claim that all
        # platforms have biologically equivalent interior components.
        alpha <- seq(
          from = 0,
          to = 1,
          length.out = nL - 1L
        )

        th2.initial <- nth1.v +
          (1 - alpha) * unmethylated.shift +
          alpha * methylated.shift

        validateThresholds(
          th2.initial,
          nL = nL,
          name = paste0("Sample ", ii, " initial thresholds"),
          require.unit.interval = FALSE
        )

        if (debug) {
          diagnostic$unmethylated_shift <-
            unmethylated.shift
          diagnostic$methylated_shift <-
            methylated.shift
          diagnostic$initial_thresholds <-
            th2.initial
        }

        # ---------------------------------------------------------------------
        # Initial hard classes and EM fit
        # ---------------------------------------------------------------------
        stage <- "initial class construction"

        sample.w0 <- initialResponsibilities(
          beta = beta2.v,
          thresholds = th2.initial,
          nL = nL
        )

        set.seed(1)

        sample.rand.idx <- sample.int(
          length(beta2.v),
          min(nfit, length(beta2.v)),
          replace = FALSE
        )

        sample.initial.class <- max.col(
          sample.w0[
            sample.rand.idx, ,
            drop = FALSE
          ],
          ties.method = "first"
        )

        sample.initial.counts <- requireAllClasses(
          class = sample.initial.class,
          nL = nL,
          context = paste0(
            "Sample ",
            ii,
            " initial mixture"
          ),
          min.count = 2L
        )

        if (debug) {
          diagnostic$random_indices <-
            sample.rand.idx
          diagnostic$initial_class_counts <-
            sample.initial.counts
        }

        stage <- "sample EM fitting"

        em2.o <- fitBetaMixture(
          beta = beta2.v[sample.rand.idx],
          initial.responsibility =
            sample.w0[
              sample.rand.idx, ,
              drop = FALSE
            ],
          nL = nL,
          maxiter = niter,
          tol = tol,
          beta.maxit = beta.maxit,
          beta.score.tol = beta.score.tol,
          debug = debug
        )

        enforceFitPolicy(
          em2.o,
          fit.policy = fit.policy,
          context = paste0("Sample ", ii, " mixture")
        )

        classAV2.v <- as.numeric(em2.o$mu[, 1L])

        stage <- "fitted component ordering"

        requireOrderedComponentMeans(
          classAV2.v,
          paste0("Sample ", ii, " mixture")
        )

        if (debug) {
          diagnostic$eta <- em2.o$eta
          diagnostic$component_means <- classAV2.v
          diagnostic$component_a <-
            as.numeric(em2.o$a[, 1L])
          diagnostic$component_b <-
            as.numeric(em2.o$b[, 1L])
          diagnostic$component_fit_status <-
            em2.o$fit_status
          diagnostic$component_fit_reason <-
            em2.o$fit_reason
          diagnostic$em_iterations <-
            em2.o$iterations
          diagnostic$em_converged <-
            em2.o$converged
          diagnostic$em_criterion <-
            em2.o$criterion
          diagnostic$log_likelihood <-
            em2.o$llike
        }

        # ---------------------------------------------------------------------
        # MAP classes and corrected midpoint thresholds
        # ---------------------------------------------------------------------
        stage <- "sample posterior MAP classes"

        subsetclass2.v <- max.col(
          em2.o$w,
          ties.method = "first"
        )

        subset.map.counts <- requireAllClasses(
          class = subsetclass2.v,
          nL = nL,
          context = paste0(
            "Sample ",
            ii,
            " posterior mixture"
          ),
          min.count = 1L
        )

        stage <- "sample posterior threshold construction"

        subsetth2.v <- thresholdsFromPosteriorClasses(
          beta = beta2.v[sample.rand.idx],
          class = subsetclass2.v,
          nL = nL,
          context = paste0(
            "Sample ",
            ii,
            " posterior mixture"
          )
        )

        # Equality remains in the lower class. This is intentionally
        # consistent and differs from some legacy final-class assignments.
        class2.v <- classifyByThresholds(
          beta = beta2.v,
          thresholds = subsetth2.v,
          nL = nL
        )

        full.class.counts <- requireAllClasses(
          class = class2.v,
          nL = nL,
          context = paste0(
            "Sample ",
            ii,
            " complete mixture"
          ),
          min.count = 1L
        )

        if (debug) {
          diagnostic$subset_map_counts <-
            subset.map.counts
          diagnostic$posterior_thresholds <-
            subsetth2.v
          diagnostic$complete_class_counts <-
            full.class.counts
        }

        nbeta2.v <- beta2.v

        # ---------------------------------------------------------------------
        # Unmethylated component
        # ---------------------------------------------------------------------
        stage <- "unmethylated quantile normalization"

        U <- 1L
        selU.idx <- which(class2.v == U)

        if (!length(selU.idx)) {
          stop(
            "Unmethylated class is empty.",
            call. = FALSE
          )
        }

        selUL.idx <- selU.idx[
          beta2.v[selU.idx] < classAV2.v[U]
        ]

        selUR.idx <- selU.idx[
          beta2.v[selU.idx] > classAV2.v[U]
        ]

        if (length(selUL.idx)) {
          probability <- pbeta(
            beta2.v[selUL.idx],
            em2.o$a[U, 1L],
            em2.o$b[U, 1L],
            lower.tail = TRUE
          )

          nbeta2.v[selUL.idx] <- qbeta(
            probability,
            em1.o$a[U, 1L],
            em1.o$b[U, 1L],
            lower.tail = TRUE
          )
        }

        if (length(selUR.idx)) {
          probability <- pbeta(
            beta2.v[selUR.idx],
            em2.o$a[U, 1L],
            em2.o$b[U, 1L],
            lower.tail = FALSE
          )

          nbeta2.v[selUR.idx] <- qbeta(
            probability,
            em1.o$a[U, 1L],
            em1.o$b[U, 1L],
            lower.tail = FALSE
          )
        }

        # ---------------------------------------------------------------------
        # Methylated component
        # ---------------------------------------------------------------------
        stage <- "methylated quantile normalization"

        M <- nL
        selM.idx <- which(class2.v == M)

        if (!length(selM.idx)) {
          stop(
            "Methylated class is empty.",
            call. = FALSE
          )
        }

        selML.idx <- selM.idx[
          beta2.v[selM.idx] < classAV2.v[M]
        ]

        selMR.idx <- selM.idx[
          beta2.v[selM.idx] > classAV2.v[M]
        ]

        if (length(selMR.idx)) {
          probability <- pbeta(
            beta2.v[selMR.idx],
            em2.o$a[M, 1L],
            em2.o$b[M, 1L],
            lower.tail = FALSE
          )

          nbeta2.v[selMR.idx] <- qbeta(
            probability,
            em1.o$a[M, 1L],
            em1.o$b[M, 1L],
            lower.tail = FALSE
          )
        }

        if (debug) {
          diagnostic$tail_counts <- c(
            U_left = length(selUL.idx),
            U_right = length(selUR.idx),
            M_left = length(selML.idx),
            M_right = length(selMR.idx)
          )
        }

        # ---------------------------------------------------------------------
        # Intermediate/H conformal transformation
        # ---------------------------------------------------------------------
        if (doH) {
          stage <- "intermediate/H normalization"

          if (!length(selMR.idx)) {
            stop(
              "No methylated observations lie above the fitted ",
              "methylated-component mean; H normalization has no ",
              "right-side methylated anchor.",
              call. = FALSE
            )
          }

          interior.components <- seq.int(2L, nL - 1L)

          # For nL > 3, all interior components are currently treated
          # collectively as H/intermediate. This is experimental.
          selH.idx <- unique(
            c(
              which(class2.v %in% interior.components),
              selML.idx
            )
          )

          if (!length(selH.idx)) {
            stop(
              "Intermediate/H normalization set is empty.",
              call. = FALSE
            )
          }

          minH <- min(beta2.v[selH.idx])
          maxH <- max(beta2.v[selH.idx])
          deltaH <- maxH - minH

          if (!is.finite(deltaH) || deltaH <= 0) {
            stop(
              "Intermediate/H values have zero or invalid range.",
              call. = FALSE
            )
          }

          deltaUH <-
            min(beta2.v[selH.idx]) -
            max(beta2.v[selU.idx])

          deltaHM <-
            min(beta2.v[selMR.idx]) -
            max(beta2.v[selH.idx])

          nmaxH <-
            min(nbeta2.v[selMR.idx]) -
            deltaHM

          nminH <-
            max(nbeta2.v[selU.idx]) +
            deltaUH

          ndeltaH <- nmaxH - nminH

          if (!is.finite(ndeltaH) || ndeltaH <= 0) {
            stop(
              "Normalized H anchors are crossed or degenerate.",
              call. = FALSE
            )
          }

          hf <- ndeltaH / deltaH

          if (!is.finite(hf) || hf <= 0) {
            stop(
              "Invalid H shift/dilation factor.",
              call. = FALSE
            )
          }

          nbeta2.v[selH.idx] <-
            nminH +
            hf * (beta2.v[selH.idx] - minH)

          if (debug) {
            diagnostic$H_count <- length(selH.idx)
            diagnostic$H_input_range <- c(minH, maxH)
            diagnostic$H_output_anchors <- c(nminH, nmaxH)
            diagnostic$H_scale <- hf
          }
        }

        # ---------------------------------------------------------------------
        # Output validation
        # ---------------------------------------------------------------------
        stage <- "sample output validation"

        if (any(!is.finite(nbeta2.v))) {
          stop(
            "Normalization produced non-finite beta values.",
            call. = FALSE
          )
        }

        # Only numerical-roundoff-sized excursions are clamped. Material
        # excursions are hard failures.
        range.tolerance <- 1e-12

        if (
          any(
            nbeta2.v < -range.tolerance |
              nbeta2.v > 1 + range.tolerance
          )
        ) {
          bad.range <- range(nbeta2.v)

          stop(
            "Normalization produced beta values outside [0, 1]. ",
            "Observed range: [",
            signif(bad.range[1L], 8),
            ", ",
            signif(bad.range[2L], 8),
            "].",
            call. = FALSE
          )
        }

        nbeta2.v <- pmin(1, pmax(0, nbeta2.v))

        if (debug) {
          diagnostic$output_range <- range(nbeta2.v)
          diagnostic$success <- TRUE
        }

        list(
          beta = nbeta2.v,
          diagnostics = diagnostic
        )
      },
      error = function(error) {
        if (inherits(error, "bmiq_sample_error")) {
          stop(error)
        }

        if (debug) {
          diagnostic$success <- FALSE
          diagnostic$failure_stage <- stage
          diagnostic$failure_message <-
            conditionMessage(error)
        }

        stop(
          newBMIQSampleError(
            sample.index = ii,
            sample.name = sample.name,
            stage = stage,
            error = error,
            diagnostics = diagnostic
          )
        )
      }
    )
  }

  # ===========================================================================
  # Sample processing loop
  # ===========================================================================
  for (ii in seq_len(number.of.samples)) {
    if (verbose) {
      message(
        "Processing sample ",
        ii,
        " of ",
        number.of.samples,
        if (nzchar(sample.names[ii])) {
          paste0(" (", sample.names[ii], ")")
        } else {
          ""
        }
      )
    }

    attempt <- tryCatch(
      processOneSample(ii),
      bmiq_sample_error = function(error) error
    )

    if (inherits(attempt, "bmiq_sample_error")) {
      if (on.sample.error == "stop") {
        stop(attempt)
      }

      success[ii] <- FALSE

      if (failed.sample == "NA") {
        calibrated[ii, ] <- NA_real_
      } else {
        calibrated[ii, ] <- original.datM[ii, ]
      }

      failures[[length(failures) + 1L]] <- data.frame(
        sample_index = attempt$sample.index,
        sample_name = attempt$sample.name,
        stage = attempt$stage,
        message = attempt$original.message,
        stringsAsFactors = FALSE
      )

      if (debug) {
        sample.diagnostics[[ii]] <-
          attempt$diagnostics
      }

      next
    }

    calibrated[ii, ] <- attempt$beta
    success[ii] <- TRUE

    if (debug) {
      sample.diagnostics[[ii]] <-
        attempt$diagnostics
    }
  }

  if (length(failures)) {
    failures <- do.call(rbind, failures)
    rownames(failures) <- NULL

    warning(
      nrow(failures),
      " sample(s) failed BMIQ calibration and were ",
      if (failed.sample == "NA") {
        "replaced with NA rows"
      } else {
        "left uncalibrated"
      },
      ". See result$failures.",
      call. = FALSE
    )
  } else {
    failures <- emptyBMIQFailureTable()
  }

  # If failed rows are NA, the complete output cannot pass the finite scanner.
  # Validate only successful rows here. Each successful row has already been
  # validated locally before assignment.
  if (all(success)) {
    scan_finite_unit_interval_cpp(
      calibrated,
      name = "calibrated datM",
      require_open = FALSE
    )
  }

  diagnostics <- if (debug) {
    list(
      gold = gold.diagnostics,
      samples = sample.diagnostics
    )
  } else {
    NULL
  }

  result <- list(
    calibrated = calibrated,
    success = success,
    failures = failures,
    diagnostics = diagnostics,
    call = call,
    settings = list(
      nL = nL,
      doH = doH,
      nfit = nfit,
      niter = niter,
      tol = tol,
      beta.maxit = beta.maxit,
      beta.score.tol = beta.score.tol,
      fit.policy = fit.policy,
      on.sample.error = on.sample.error,
      failed.sample = failed.sample,
      debug = debug
    )
  )

  class(result) <- "BMIQcalibration_result"
  result
}


# -----------------------------------------------------------------------------
# Convenience methods for the structured result.
# -----------------------------------------------------------------------------
print.BMIQcalibration_result <- function(x, ...) {
  number.of.samples <- length(x$success)
  number.succeeded <- sum(x$success)
  number.failed <- number.of.samples - number.succeeded

  cat("BMIQ calibration result\n")
  cat("  Samples:   ", number.of.samples, "\n", sep = "")
  cat("  Succeeded: ", number.succeeded, "\n", sep = "")
  cat("  Failed:    ", number.failed, "\n", sep = "")

  if (number.failed) {
    cat("\nFailures:\n")
    print(x$failures, row.names = FALSE)
  }

  invisible(x)
}

as.matrix.BMIQcalibration_result <- function(x, ...) {
  x$calibrated
}
