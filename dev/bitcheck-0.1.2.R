# Bit-identity check for the 0.1.2 cleanup (see dev/0.1.2.md).
#
# Usage, from the package root (R_LIBS_USER must reach Rcpp and checkmate):
#   Rscript dev/bitcheck-0.1.2.R <outdir>                 # capture only
#   Rscript dev/bitcheck-0.1.2.R <outdir> <baseline-dir>  # capture + compare
#
# Copies the working tree (tracked plus untracked, not-ignored files, minus
# dev/) into <outdir>/src so no object file from an earlier build (load_all
# builds at -O0) is reused, installs it with standard R CMD INSTALL flags
# into <outdir>/lib, runs every case, and saves <outdir>/results.rds.
#
# Compare mode checks each case against the baseline with
# identical(num.eq = FALSE), i.e. bitwise on doubles, and reruns the
# calibration on the gold fits saved by the baseline (a prefit made with the
# older version must give the older numbers). The random-state footprint is
# reported apart: phase 1 changes it on purpose.

args <- commandArgs(trailingOnly = TRUE)
stopifnot(
  length(args) %in% 1:2,
  file.exists("DESCRIPTION"),
  read.dcf("DESCRIPTION", "Package")[[1L]] == "betanorm"
)
outdir <- normalizePath(args[[1L]], winslash = "/", mustWork = FALSE)
baseline.dir <- if (length(args) == 2L) {
  normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
}

# ---- install --------------------------------------------------------------

src <- file.path(outdir, "src")
lib <- file.path(outdir, "lib")
unlink(c(src, lib), recursive = TRUE)
dir.create(src, recursive = TRUE)
dir.create(lib, recursive = TRUE)

files <- system2(
  "git",
  c("ls-files", "--cached", "--others", "--exclude-standard"),
  stdout = TRUE
)
files <- unique(files[file.exists(files) & !startsWith(files, "dev/")])
for (d in unique(dirname(files))) {
  dir.create(file.path(src, d), recursive = TRUE, showWarnings = FALSE)
}
stopifnot(all(file.copy(files, file.path(src, files), copy.date = TRUE)))

status <- system2(
  file.path(R.home("bin"), "R"),
  c("CMD", "INSTALL", "--no-multiarch", paste0("--library=", lib), src)
)
stopifnot(status == 0L)

library(betanorm, lib.loc = lib)
stopifnot(
  normalizePath(dirname(getNamespaceInfo("betanorm", "path")), "/") ==
    normalizePath(lib, "/")
)

provenance <- list(
  head = system2("git", c("rev-parse", "HEAD"), stdout = TRUE),
  status = system2("git", c("status", "--porcelain"), stdout = TRUE),
  version = as.character(utils::packageVersion("betanorm", lib.loc = lib)),
  makevars = intersect(
    c("src/Makevars", "src/Makevars.win"),
    files
  ),
  r = R.version.string
)

# ---- helpers --------------------------------------------------------------

# Copies of tests/testthat/helper-synthetic.R, so the synthetic cases match
# the snapshot tests.
clip01 <- function(x, eps = 1e-4) {
  pmin(1 - eps, pmax(eps, x))
}

simulate_beta_mixture <- function(n, nL = 3L, seed = 1L) {
  set.seed(seed)
  if (nL == 2L) {
    n_u <- floor(n / 2)
    n_m <- n - n_u
    y <- c(rbeta(n_u, 2, 8), rbeta(n_m, 8, 2))
  } else {
    n_u <- floor(0.4 * n)
    n_m <- floor(0.4 * n)
    n_h <- n - n_u - n_m
    y <- c(rbeta(n_u, 2, 10), rbeta(n_h, 3, 3), rbeta(n_m, 10, 2))
  }
  clip01(y)
}

hard_responsibility <- function(y, thresholds, nL = 3L) {
  class <- betanorm:::class_by_thresh(y, thresholds)
  w <- matrix(0, nrow = length(y), ncol = nL)
  w[cbind(seq_along(y), class)] <- 1
  w
}

make_calibration_inputs <- function(n_probes = 1500L, nL = 3L, seed = 42L) {
  gold <- simulate_beta_mixture(n_probes, nL = nL, seed = seed)
  sample <- simulate_beta_mixture(n_probes, nL = nL, seed = seed + 7L)
  sample <- clip01(sample + 0.02)
  datM <- matrix(sample, nrow = 1L, dimnames = list("s1", NULL))
  list(datM = datM, gold = gold)
}

# Value, printed output, warnings and messages of f(); an error becomes its
# message and class.
run_case <- function(f) {
  warnings <- character()
  messages <- character()
  value <- NULL
  failed <- FALSE
  output <- withCallingHandlers(
    tryCatch(
      utils::capture.output(value <- f()),
      error = function(e) {
        failed <<- TRUE
        value <<- list(message = conditionMessage(e), class = class(e))
        character()
      }
    ),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    },
    message = function(m) {
      messages <<- c(messages, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )
  list(
    error = failed,
    value = value,
    output = output,
    warnings = warnings,
    messages = messages
  )
}

# Does f() leave a .Random.seed in a session that has none, and does it
# leave an existing one unchanged?
rng_footprint <- function(f) {
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    rm(".Random.seed", envir = globalenv())
  }
  suppressWarnings(suppressMessages(f()))
  created <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  set.seed(2026L)
  before <- get(".Random.seed", envir = globalenv())
  suppressWarnings(suppressMessages(f()))
  after <- get0(".Random.seed", envir = globalenv(), inherits = FALSE)
  c(created_when_absent = created, unchanged_when_set = identical(before, after))
}

# Wrap fit_mixture so sample 2's H map fails (its methylated mean is pushed
# above every probe), as in tests/testthat/test-ordering-h-policy.R.
real_fit_mixture <- get("fit_mixture", envir = asNamespace("betanorm"))
stopifnot(identical(
  names(formals(real_fit_mixture)),
  c(
    "beta", "thresholds", "nL", "nfit", "niter", "tol", "beta.maxit",
    "beta.score.tol", "context", "debug", "seed", "fit.idx"
  )
))
h_failure_fit_mixture <- function(
  beta,
  thresholds,
  nL,
  nfit,
  niter,
  tol,
  beta.maxit,
  beta.score.tol,
  context,
  debug = FALSE,
  seed = 1L,
  fit.idx = NULL
) {
  fit <- real_fit_mixture(
    beta = beta,
    thresholds = thresholds,
    nL = nL,
    nfit = nfit,
    niter = niter,
    tol = tol,
    beta.maxit = beta.maxit,
    beta.score.tol = beta.score.tol,
    context = context,
    debug = debug,
    seed = seed,
    fit.idx = fit.idx
  )
  if (identical(context, "Sample 2") && nL == 3L) {
    fit$component_means <- c(0.15, 0.5, 1.5)
  }
  fit
}
with_h_failure <- function(f) {
  function() {
    utils::assignInNamespace("fit_mixture", h_failure_fit_mixture, "betanorm")
    on.exit(
      utils::assignInNamespace("fit_mixture", real_fit_mixture, "betanorm")
    )
    f()
  }
}

# ---- inputs ---------------------------------------------------------------

gold.real <- betanorm::horvath_goldstandard$goldstandard2
datM.real <- betanorm::GPL21145_sample

syn3 <- make_calibration_inputs(n_probes = 1500L, nL = 3L, seed = 42L)
syn2 <- make_calibration_inputs(n_probes = 1500L, nL = 2L, seed = 43L)

block.gold <- simulate_beta_mixture(1200L, nL = 3L, seed = 500L)
block.datM <- t(vapply(
  seq_len(11L),
  function(i) clip01(simulate_beta_mixture(1200L, nL = 3L, seed = 500L + i)),
  numeric(1200L)
))
rownames(block.datM) <- paste0("s", seq_len(11L))

opt.gold <- simulate_beta_mixture(1500L, nL = 3L, seed = 11L)
opt.datM <- t(vapply(
  seq_len(3L),
  function(i) clip01(simulate_beta_mixture(1500L, nL = 3L, seed = 20L + i) + 0.02),
  numeric(1500L)
))
rownames(opt.datM) <- paste0("s", seq_len(3L))

em3.y <- simulate_beta_mixture(2000L, nL = 3L, seed = 1234L)
em3.w0 <- hard_responsibility(em3.y, c(0.25, 0.75), 3L)
em2.y <- simulate_beta_mixture(2000L, nL = 2L, seed = 4321L)
em2.w0 <- hard_responsibility(em2.y, 0.5, 2L)

set.seed(7L)
qn.obj <- matrix(rnorm(6L * 40L), nrow = 6L, ncol = 40L)
qn.target <- sort(rnorm(40L))
qn.target.short <- sort(rnorm(17L))
qn.tied <- round(datM.real, 2L)

# ---- cases ----------------------------------------------------------------

gold_fit <- function(gold, nL) {
  bmiq_gold_fit(goldstandard.beta = gold, nL = nL, verbose = FALSE)
}
calib_with_prefit <- function(prefit, nL) {
  bmiq_calibration(
    datM = datM.real,
    goldstandard.beta = prefit,
    nL = nL,
    debug = TRUE,
    verbose = FALSE
  )
}
calib_options <- function(h.policy, on.sample.error, failed.sample) {
  bmiq_calibration(
    datM = opt.datM,
    goldstandard.beta = opt.gold,
    doH = TRUE,
    nfit = 1500L,
    h.policy = h.policy,
    on.sample.error = on.sample.error,
    failed.sample = failed.sample,
    debug = TRUE,
    verbose = TRUE
  )
}

cases <- list(
  gold_fit_real_nL3 = function() gold_fit(gold.real, 3L),
  gold_fit_real_nL2 = function() gold_fit(gold.real, 2L),
  gold_fit_real_nL3_debug = function() {
    bmiq_gold_fit(gold.real, nL = 3L, debug = TRUE, verbose = TRUE)
  },
  calib_real_nL3_raw = function() {
    bmiq_calibration(datM.real, gold.real, nL = 3L, debug = TRUE, verbose = TRUE)
  },
  calib_real_nL2_raw = function() {
    bmiq_calibration(datM.real, gold.real, nL = 2L, debug = TRUE, verbose = TRUE)
  },
  calib_real_nL3_prefit = function() {
    calib_with_prefit(gold_fit(gold.real, 3L), 3L)
  },
  calib_real_nL2_prefit = function() {
    calib_with_prefit(gold_fit(gold.real, 2L), 2L)
  },
  calib_real_nL3_print = function() {
    print(bmiq_calibration(datM.real, gold.real, verbose = FALSE))
  },
  calib_syn_nL3 = function() {
    bmiq_calibration(
      datM = syn3$datM,
      goldstandard.beta = syn3$gold,
      nL = 3L,
      nfit = 1500L,
      niter = 5L,
      verbose = FALSE
    )
  },
  calib_syn_nL2 = function() {
    bmiq_calibration(
      datM = syn2$datM,
      goldstandard.beta = syn2$gold,
      nL = 2L,
      nfit = 1500L,
      niter = 10L,
      verbose = FALSE
    )
  },
  calib_syn_blocks = function() {
    bmiq_calibration(
      datM = block.datM,
      goldstandard.beta = block.gold,
      nL = 3L,
      nfit = 1200L,
      debug = TRUE,
      verbose = FALSE
    )
  },
  calib_syn_niter8 = function() {
    bmiq_calibration(
      datM = syn3$datM,
      goldstandard.beta = syn3$gold,
      nfit = 1500L,
      niter = 8L,
      verbose = FALSE
    )
  },
  h_optional = with_h_failure(function() {
    calib_options("optional", "stop", "NA")
  }),
  h_require_stop = with_h_failure(function() {
    calib_options("require", "stop", "NA")
  }),
  h_require_continue_na = with_h_failure(function() {
    calib_options("require", "continue", "NA")
  }),
  h_require_continue_original = with_h_failure(function() {
    calib_options("require", "continue", "original")
  }),
  h_require_continue_print = with_h_failure(function() {
    print(calib_options("require", "continue", "NA"))
  }),
  em_nL3 = function() {
    betanorm:::beta_mixture_em_cpp(
      y = em3.y,
      initial_responsibility = em3.w0,
      nL = 3L,
      maxiter = 25L,
      tol = 1e-6
    )
  },
  em_nL2 = function() {
    betanorm:::beta_mixture_em_cpp(
      y = em2.y,
      initial_responsibility = em2.w0,
      nL = 2L,
      maxiter = 40L,
      tol = 1e-8
    )
  },
  qn_equal = function() quantile_norm(qn.obj, qn.target),
  qn_unequal = function() quantile_norm(qn.obj, qn.target.short),
  qn_real = function() quantile_norm(datM.real, gold.real),
  qn_real_unequal = function() quantile_norm(datM.real, gold.real[1:5000]),
  qn_tied = function() quantile_norm(qn.tied, round(gold.real, 2L))
)

results <- lapply(cases, run_case)

rng <- rbind(
  calib_raw = rng_footprint(function() {
    bmiq_calibration(syn3$datM, syn3$gold, nfit = 1500L, verbose = FALSE)
  }),
  calib_prefit = rng_footprint(function() {
    bmiq_calibration(
      syn3$datM,
      bmiq_gold_fit(syn3$gold, nfit = 1500L, verbose = FALSE),
      nfit = 1500L,
      verbose = FALSE
    )
  }),
  gold_fit = rng_footprint(function() {
    bmiq_gold_fit(syn3$gold, nfit = 1500L, verbose = FALSE)
  }),
  quantile_norm = rng_footprint(function() quantile_norm(qn.obj, qn.target))
)

saveRDS(
  list(provenance = provenance, cases = results, rng = rng),
  file.path(outdir, "results.rds")
)

cat("\n== provenance ==\n")
str(provenance)
cat("\n== cases ==\n")
for (name in names(results)) {
  r <- results[[name]]
  cat(sprintf(
    "%-30s %s  warnings=%d messages=%d\n",
    name,
    if (r$error) "ERROR" else "ok",
    length(r$warnings),
    length(r$messages)
  ))
}
cat("\n== random-state footprint ==\n")
print(rng)

# ---- compare --------------------------------------------------------------

if (!is.null(baseline.dir)) {
  baseline <- readRDS(file.path(baseline.dir, "results.rds"))

  # Rerun on the gold fits the baseline saved: a prefit from the older
  # version must reproduce the older calibration.
  compat <- list(
    calib_real_nL3_prefit = run_case(function() {
      calib_with_prefit(baseline$cases$gold_fit_real_nL3$value, 3L)
    }),
    calib_real_nL2_prefit = run_case(function() {
      calib_with_prefit(baseline$cases$gold_fit_real_nL2$value, 2L)
    })
  )

  # The value (and whether the case errored) must be bit-identical. Printed
  # output, warnings and messages are reported apart, old then new, so a
  # deliberate wording change can be reviewed without masking a number.
  report <- function(label, new, old) {
    value.same <- identical(new$error, old$error) &&
      identical(new$value, old$value, num.eq = FALSE)
    text <- c("output", "warnings", "messages")
    text.diff <- text[
      !vapply(text, function(k) identical(new[[k]], old[[k]]), logical(1L))
    ]
    status <- if (!value.same) {
      "VALUE DIFFERS"
    } else if (length(text.diff)) {
      paste("value identical; text differs:", paste(text.diff, collapse = ", "))
    } else {
      "identical"
    }
    cat(sprintf("%-45s %s\n", label, status))
    if (!value.same) {
      cat(paste0("    ", utils::head(all.equal(old, new), 10L)), sep = "\n")
    }
    for (k in text.diff) {
      cat("    ", k, " was: ", sep = "")
      cat(setdiff(old[[k]], new[[k]]), sep = "\n             ")
      cat("\n    ", k, " now: ", sep = "")
      cat(setdiff(new[[k]], old[[k]]), sep = "\n             ")
      cat("\n")
    }
    c(value = value.same, text = !length(text.diff))
  }

  cat("\n== compare with", baseline.dir, "==\n")
  cat("baseline:", baseline$provenance$version, baseline$provenance$head, "\n")
  all.names <- union(names(baseline$cases), names(results))
  same <- vapply(
    all.names,
    function(name) report(name, results[[name]], baseline$cases[[name]]),
    logical(2L)
  )
  same.compat <- vapply(
    names(compat),
    function(name) {
      report(
        paste0(name, " (baseline gold fit)"),
        compat[[name]],
        baseline$cases[[name]]
      )
    },
    logical(2L)
  )
  cat("\nrandom-state footprint, baseline then now:\n")
  print(baseline$rng)
  print(rng)

  same <- cbind(same, same.compat)
  n.diff <- sum(!same["value", ])
  cat(sprintf(
    "\n%s: %d of %d case values differ; %d case(s) with text changes.\n",
    if (n.diff) "FAIL" else "PASS",
    n.diff,
    ncol(same),
    sum(!same["text", ])
  ))
  quit(status = if (n.diff) 1L else 0L)
}
