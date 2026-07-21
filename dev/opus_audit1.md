Overall assessment

The weighted Beta MLE, EM E-step, component canonicalization, and the nL = 2 truncated quantile map are mathematically sound.

The main qualifications are:

    The nL = 3 thresholds are adjacent pairwise density crossings, not necessarily global posterior/MAP class boundaries.
    The nL = 3 H transformation is a legacy empirical stitching heuristic, not a population-level three-component quantile map.
    A few numerical and EM-policy details can produce avoidable failures.
    For a 100 × 94,000 matrix, the biggest performance issue is not its total size but that samples are rows in R’s column-major layout. Blocked row gathering/scattering is the best memory/speed compromise.

1. Mathematical audit
Correct parts
Weighted Beta MLE

Your sufficient-statistic likelihood is correct:
[ \ell(a,b)

(a-1)\sum_iw_i\log y_i
+
(b-1)\sum_iw_i\log(1-y_i)
+
\left(\sum_iw_i\right)
\left[\log\Gamma(a+b)-\log\Gamma(a)-\log\Gamma(b)\right].
]

The score and expected-information matrix used in fit_beta_from_stats() are also correct. Solving

[
I(a,b)\Delta=s(a,b)
]

gives an ascent Newton direction because the Hessian is -weight * I.

The Armijo line search and positivity checks are appropriate.
EM E-step

The mixture log-density,

log(eta[k]) +
lgamma(a[k] + b[k]) - lgamma(a[k]) - lgamma(b[k]) +
(a[k] - 1) * log(y) +
(b[k] - 1) * log1p(-y)

and the log-sum-exp normalization are correct.
Component canonicalization

Sorting complete component tuples by increasing mean is correct. It is particularly important that you sort a, b, eta, mu, responsibilities, and statuses together rather than independently.
nL = 2 truncated map

The formulas are correct.

For the U side,

[
u(x)=\frac{F_{s,U}(x)}{F_{s,U}(t_s)},\qquad
y=F_{g,U}^{-1}\left[u(x)F_{g,U}(t_g)\right].
]

For the M side,

[
r(x)=\frac{\bar F_{s,M}(x)}{\bar F_{s,M}(t_s)},\qquad
y=\bar F_{g,M}^{-1}\left[r(x)\bar F_{g,M}(t_g)\right].
]

Both sides map the sample threshold to the gold threshold, so the map is continuous there. It remains possible to have a slope discontinuity, as your documentation states.
2. Important mathematical qualifications
2.1 The three-component thresholds are not necessarily MAP thresholds

For adjacent components (k) and (k+1), density_thresholds() solves

[
\eta_k f_k(x)=\eta_{k+1}f_{k+1}(x)
]

with the correct lower-to-higher orientation. That calculation is correct.

However, for nL = 3, a crossing between components 1 and 2 need not be a global posterior boundary if component 3 dominates both at that point. Likewise, component 2 can have pairwise crossings with components 1 and 3 while never being the global posterior winner.

Consequently:

full.class <- class_by_thresh(beta, posterior.thresholds)

is an ordinal classification induced by adjacent pairwise cuts, not necessarily

max.col(posterior_responsibilities)

for the fitted mixture.

This is especially relevant because you explicitly allow:

subset_map_counts[k] == 0

while still constructing a nonempty threshold-defined class for that component.
Recommended decision

Choose one of these semantics:
Option A: Keep the current behavior

Rename the concepts to make them precise:

    density_thresholds() → adjacent_density_thresholds()
    posterior.thresholds → adjacent.thresholds
    “posterior class” → “ordinal threshold class”

Document that these are adjacent weighted-density crossings used to impose an ordinal partition.
Option B: Require genuine posterior/MAP ordering

Compute or validate the upper envelope of all component weighted densities and require that the global winner sequence is exactly:

1 → 2 → 3

with two transitions. If component 2 never globally dominates, reject the three-component fit or fall back to a two-component model.

Merely checking the two adjacent equalities is insufficient for this.
2.2 Restricting crossings to the interval between component means is an assumption

You search only between:

means[k]
means[k + 1L]

A valid weighted-density crossing can occur outside the interval between the two means, especially with unequal mixture weights or strongly different concentration parameters.

Your function safely rejects such a fit rather than silently using the wrong crossing, which is good. But this is a model restriction, not a general property of Beta mixtures. It should be documented.
2.3 The three-component H step is empirical stitching

When H succeeds, you map:

    the minimum observed H/lower-M value to the maximum normalized U value;
    the maximum observed H/lower-M value to the minimum normalized upper-M value.

This makes the transformed observed order statistics join exactly. It does not define a population-level continuous map at the theoretical density thresholds.

Thus “continuous H stitching” is reasonable as a finite-sample description, but it is not the same mathematical construction as the continuous nL = 2 threshold map.

Also, when optional H normalization fails, both of these remain unchanged:

    the fitted H class;
    the lower half of the M class.

The documentation currently emphasizes lower-M observations but should also explicitly mention that H observations remain unchanged.

Suggested wording:

    When H normalization is skipped, U and upper-M observations are calibrated; the entire H class and lower-M observations are left unchanged.

2.4 The U tail split is mathematically redundant

For exact arithmetic,

qbeta(
  pbeta(x, a.sample, b.sample, lower.tail = FALSE),
  a.gold, b.gold,
  lower.tail = FALSE
)

is the same monotone quantile map as the corresponding lower-tail expression.

Therefore splitting U at its component mean does not change the mathematical map. It can still be useful numerically because it evaluates the smaller tail.

The same observation applies generally to component-wise full-distribution quantile mapping.
3. Numerical and fitting issues to fix
3.1 Endpoint clipping can still produce exact 0 or 1

This code can underflow or round back to an endpoint:

lower.clip <- min(positive) / 2
upper.clip <- 1 - (1 - max(below.one)) / 2

For example:

    half the smallest positive subnormal can become zero;
    1 - half_the_distance can round to exactly one.

Replace the clipping section in fit_mixture()

Replace:

lower.clip <- if (length(positive)) min(positive) / 2 else endpoint.eps
upper.clip <- if (length(below.one)) {
  1 - (1 - max(below.one)) / 2
} else {
  1 - endpoint.eps
}
y_fit <- pmin(upper.clip, pmax(lower.clip, y_fit))

with:

lower.clip <- if (length(positive)) {
  max(endpoint.eps, min(positive) / 2)
} else {
  endpoint.eps
}

upper.clip <- if (length(below.one)) {
  min(
    1 - endpoint.eps,
    1 - (1 - max(below.one)) / 2
  )
} else {
  1 - endpoint.eps
}

if (!(lower.clip < upper.clip)) {
  stop(
    context,
    " could not construct a valid open-interval clipping range.",
    call. = FALSE
  )
}

y_fit <- pmin(upper.clip, pmax(lower.clip, y_fit))

3.2 "usable" does not guarantee generalized-EM ascent

Each Beta M-step starts from a new method-of-moments estimate. If Newton stalls, you accept that estimate under fit.policy = "usable".

Although the accepted Newton steps ascend relative to the current method-of-moments initialization, they are not guaranteed to improve the EM auxiliary function relative to the previous outer iteration’s component parameters. Therefore the outer likelihood is not guaranteed to be monotone under "usable".
More rigorous alternatives

In increasing order of rigor:

    Record a log-likelihood trace and warn if it decreases materially.
    Warm-start each Beta fit from the previous outer iteration’s a and b.
    Compare the candidate component weighted log-likelihood against the previous component parameters and retain the previous parameters if the candidate is worse.
    Require fully converged component M-steps.

I recommend at least options 1 and 3.
3.3 Effective component size is not actually checked

positive_count counts observations whose responsibility is greater than zero. After a soft E-step, almost every observation can have a positive but negligible responsibility.

Thus:

stats.positive_count > 1

does not ensure a meaningful effective sample size.

Consider accumulating:

sum_w2 += effective_weight * effective_weight;

and defining:

[
n_{\mathrm{eff}} = \frac{(\sum_i w_i)^2}{\sum_iw_i^2}.
]

Reject or flag components with, for example:

effective sample size < 2

or a very small mixture weight.
3.4 Extreme-tail underflow in normalize_nl2()

pbeta() can return exactly zero for a mathematically nonzero extreme tail. Your current behavior is to reject the map.

For better numerical robustness, calculate the ratios in log-probability space:

log.FsU.ts <- pbeta(
  sample.threshold,
  sample.a[1L],
  sample.b[1L],
  lower.tail = TRUE,
  log.p = TRUE
)

log.FsM.ts <- pbeta(
  sample.threshold,
  sample.a[2L],
  sample.b[2L],
  lower.tail = FALSE,
  log.p = TRUE
)

Then, for example:

log.u <- pbeta(
  beta[u_idx],
  sample.a[1L],
  sample.b[1L],
  lower.tail = TRUE,
  log.p = TRUE
) - log.FsU.ts

log.u <- pmin(0, log.u)

out[u_idx] <- qbeta(
  log.u + log.FgU.tg,
  gold.a[1L],
  gold.b[1L],
  lower.tail = TRUE,
  log.p = TRUE
)

Use the analogous upper-tail expression for M. This also handles endpoint probabilities cleanly.
3.5 Kernel mode estimates can be outside [0,1]

Base density() is not boundary corrected. For beta values concentrated near zero or one, its estimated mode can fall outside the unit interval.

That mode is only used to shift initialization thresholds, so it does not directly invalidate the final fit, but it can create poor or impossible initial classes.

At minimum, constrain the result:

mode <- estimate$x[which.max(estimate$y)]
pmin(1, pmax(0, mode))

A boundary-corrected density estimate or component-based mode estimate would be more principled.
3.6 Random subsampling can fail despite adequate full-data classes

You randomly sample first and then demand two observations from every initial class. A rare class can therefore fail by chance even if the full dataset contains many observations from it.

A more robust approach is stratified sampling:

    classify the full vector using the initialization thresholds;
    sample at least two observations from every class;
    allocate the remaining sample budget proportionally or uniformly.

This changes the fit sample distribution slightly, so whether to do it depends on the desired legacy compatibility.

Also, repeated set.seed(seed) changes the caller’s global RNG state. Prefer withr::with_seed() or preserve and restore .Random.seed.
4. Easy memory improvement inside fit_mixture()

You currently allocate a full length(beta) × nL one-hot matrix and then retain only sampled rows.

For 94,000 probes and three classes that is about 2.3 MB per call, plus the subset copy. It is unnecessary.
Replace the start of fit_mixture()

Replace:

class0 <- class_by_thresh(beta, thresholds)
w0 <- matrix(0, nrow = length(beta), ncol = nL)
w0[cbind(seq_along(beta), class0)] <- 1

set.seed(seed)
rand.idx <- sample.int(
  length(beta),
  min(nfit, length(beta)),
  replace = FALSE
)

initial.counts <- require_all_classes(
  class = max.col(w0[rand.idx, , drop = FALSE], ties.method = "first"),
  nL = nL,
  context = paste0(context, " initial mixture"),
  min.count = 2L
)

with:

set.seed(seed)
rand.idx <- sample.int(
  length(beta),
  min(nfit, length(beta)),
  replace = FALSE
)

initial.class <- class_by_thresh(
  beta[rand.idx],
  thresholds
)

initial.counts <- require_all_classes(
  class = initial.class,
  nL = nL,
  context = paste0(context, " initial mixture"),
  min.count = 2L
)

w.init <- matrix(
  0,
  nrow = length(rand.idx),
  ncol = nL
)
w.init[cbind(seq_along(rand.idx), initial.class)] <- 1

Then replace:

initial_responsibility = w0[rand.idx, , drop = FALSE],

with:

initial_responsibility = w.init,

After the C++ call, these can be released:

rm(w.init, initial.class, y_fit)

5. Do not retain the EM responsibility matrix

After canonicalization, em$w is used only for:

subset.class <- max.col(em$w, ties.method = "first")
subset.counts <- tabulate(subset.class, nbins = nL)

Your diagnostics do not return the responsibility matrix.

Immediately after computing the counts, add:

em$w <- NULL

So this section becomes:

subset.class <- max.col(em$w, ties.method = "first")
subset.counts <- tabulate(subset.class, nbins = nL)

em$w <- NULL
rm(subset.class)

An even better version is to calculate MAP counts in C++ and not return w at all. That avoids copying the Armadillo responsibility matrix into an R matrix at function return.
6. Matrix-layout audit

Your input is:

100 samples × 94,000 probes

That contains 9.4 million doubles:

    about 75.2 MB in decimal units;
    about 71.7 MiB.

The size itself is manageable. The layout is the main issue.

R matrices are column-major. Therefore:

    all 100 samples for one probe are contiguous;
    one sample across 94,000 probes is strided by 100 doubles, or 800 bytes.

This operation:

datM[ii, ]

touches one value from many different cache lines. Repeating it for every sample causes substantially more memory traffic than the nominal 75 MB matrix size suggests.
Do not fully transpose if memory is the primary goal

A complete transpose makes each sample contiguous, but it introduces another approximately 75 MB matrix. Transposing back for the required result orientation introduces another full-size allocation.

The fastest layout would be probes in rows and samples in columns throughout the entire pipeline. If you can change the upstream and downstream API, that is ideal.

If you must retain samples × probes, use small blocked gathers and scatters.
7. Recommended blocked layout

Use blocks of eight samples:

    Read eight neighboring sample rows together.
    Store them internally as a probes × 8 matrix.
    Process each contiguous block column.
    Scatter the block back into the output matrix.

Eight samples correspond to eight doubles, matching a typical 64-byte cache line. The temporary block is only:

[
94{,}000 \times 8 \times 8 \approx 6\text{ MB}.
]

This avoids a full transpose while making input and output traffic cache-friendly.
8. C++ block helpers to copy and paste

Add these after scan_finite_unit_interval_cpp() in the C++ file:

// [[Rcpp::export]]
Rcpp::NumericMatrix gather_sample_block_cpp(
    const Rcpp::NumericMatrix &x,
    int first_sample,
    int sample_count)
{
    const int n_samples = x.nrow();
    const int n_probes = x.ncol();
    const int first0 = first_sample - 1;

    if (first0 < 0 ||
        sample_count < 1 ||
        first0 + sample_count > n_samples)
    {
        Rcpp::stop("Invalid sample block");
    }

    // Rows are probes; columns are samples. Each sample is therefore
    // contiguous in the returned column-major matrix.
    Rcpp::NumericMatrix block(n_probes, sample_count);

    const double *x_ptr = REAL(x);
    double *block_ptr = REAL(block);

    for (int probe = 0; probe < n_probes; ++probe)
    {
        const R_xlen_t x_offset =
            static_cast<R_xlen_t>(n_samples) * probe;

        for (int local_sample = 0;
             local_sample < sample_count;
             ++local_sample)
        {
            block_ptr[
                probe +
                static_cast<R_xlen_t>(n_probes) * local_sample
            ] =
                x_ptr[
                    x_offset +
                    first0 +
                    local_sample
                ];
        }
    }

    return block;
}


// [[Rcpp::export]]
void scatter_sample_block_cpp(
    Rcpp::NumericMatrix destination,
    const Rcpp::NumericMatrix &block,
    int first_sample)
{
    const int n_samples = destination.nrow();
    const int n_probes = destination.ncol();
    const int sample_count = block.ncol();
    const int first0 = first_sample - 1;

    if (block.nrow() != n_probes ||
        first0 < 0 ||
        first0 + sample_count > n_samples)
    {
        Rcpp::stop("Invalid destination or sample block");
    }

    double *destination_ptr = REAL(destination);
    const double *block_ptr = REAL(block);

    for (int probe = 0; probe < n_probes; ++probe)
    {
        const R_xlen_t destination_offset =
            static_cast<R_xlen_t>(n_samples) * probe;

        for (int local_sample = 0;
             local_sample < sample_count;
             ++local_sample)
        {
            destination_ptr[
                destination_offset +
                first0 +
                local_sample
            ] =
                block_ptr[
                    probe +
                    static_cast<R_xlen_t>(n_probes) * local_sample
                ];
        }
    }
}

scatter_sample_block_cpp() deliberately mutates its destination. This is safe here only if calibrated is a newly allocated private result matrix. Do not use it to mutate datM, because that could modify the caller’s matrix.

Regenerate the Rcpp exports afterward.
9. R-side locations to change for blocked processing
9.1 Replace the full-matrix aliases

Replace:

original.datM <- datM
calibrated <- datM

with:

calibrated <- matrix(
  NA_real_,
  nrow = number.of.samples,
  ncol = number.of.probes,
  dimnames = dimnames(datM)
)

The current assignments use copy-on-write, so they do not immediately create three matrices, but the first modification of calibrated creates a full output copy. Explicit allocation makes ownership clear and allows safe C++ scatter operations.
9.2 Pass beta values into process_sample()

Replace:

process_sample <- function(ii) {
  beta2.v <- as.numeric(original.datM[ii, ])

with:

process_sample <- function(ii, beta2.v) {
  beta2.v <- as.numeric(beta2.v)

9.3 Replace the single-sample loop with a block loop

Replace:

for (ii in seq_len(number.of.samples)) {
  ...
}

with this structure:

sample.block.size <- 8L

for (block.start in seq.int(
  1L,
  number.of.samples,
  by = sample.block.size
)) {
  block.count <- min(
    sample.block.size,
    number.of.samples - block.start + 1L
  )

  block <- gather_sample_block_cpp(
    datM,
    first_sample = block.start,
    sample_count = block.count
  )

  for (local.sample in seq_len(block.count)) {
    ii <- block.start + local.sample - 1L

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
      process_sample(ii, block[, local.sample]),
      bmiq_sample_error = function(error) error
    )

    if (inherits(attempt, "bmiq_sample_error")) {
      if (on.sample.error == "stop") {
        stop(attempt)
      }

      success[ii] <- FALSE

      if (failed.sample == "NA") {
        block[, local.sample] <- NA_real_
      }
      # For failed.sample == "original", leave the source block column
      # unchanged.

      failures[[length(failures) + 1L]] <- data.frame(
        sample_index = attempt$sample.index,
        sample_name = attempt$sample.name,
        stage = attempt$stage,
        message = attempt$original.message,
        stringsAsFactors = FALSE
      )

      if (debug) {
        sample.diagnostics[[ii]] <- attempt$diagnostics
      }

      next
    }

    block[, local.sample] <- attempt$beta
    success[ii] <- TRUE

    if (doH) {
      h.applied.vec[ii] <- attempt$h_applied
    }

    if (debug) {
      sample.diagnostics[[ii]] <- attempt$diagnostics
    }
  }

  scatter_sample_block_cpp(
    destination = calibrated,
    block = block,
    first_sample = block.start
  )

  rm(block)
}

This preserves the existing failure semantics while reducing cache-unfriendly row reads and writes.
10. Values to hoist out of the sample loop

Immediately after:

em1.o <- gold.fit$em
nth1.v <- gold.fit$thresholds

add:

gold.a <- as.numeric(em1.o$a[, 1L])
gold.b <- as.numeric(em1.o$b[, 1L])
gold.thresholds <- as.numeric(nth1.v)

Then use these inside process_sample() instead of repeatedly extracting:

as.numeric(em1.o$a[, 1L])
as.numeric(em1.o$b[, 1L])

For example:

gold.a = gold.a,
gold.b = gold.b,
gold.threshold = gold.thresholds[1L]

and in the three-component map:

gold.a[U]
gold.b[U]
gold.a[M]
gold.b[M]

The memory saving is small, but it makes the sample loop cleaner and avoids repeated allocations.
11. Small C++ hoisting improvements

These are currently allocated inside every EM iteration:

arma::vec log_prior(K);
arma::vec am1(K);
arma::vec bm1(K);

Move them above the outer loop, next to log_component and log_norm:

arma::vec log_component(K);
arma::vec log_norm(K);
arma::vec log_prior(K);
arma::vec am1(K);
arma::vec bm1(K);

Then remove their declarations inside the loop.

For K <= 3, the performance effect is minor. The matrix-layout changes and avoiding w0 are much more important.
12. Recommended testing before relying on the result

Add tests for:

    Likelihood monotonicity: record the outer log-likelihood trace, especially under "usable".
    Label invariance: permuting initial component columns should give the same canonicalized fit.
    Threshold roots: verify weighted log-density equality and negative crossing slope.
    Global dominance: for nL = 3, construct a case where H has pairwise crossings but is never the MAP component.
    nL = 2 map: check endpoints, threshold continuity, monotonicity, and identity when sample and gold parameters match.
    Extreme tails: use highly concentrated Beta components and thresholds near zero or one.
    Endpoint inputs: include exact 0, exact 1, and nearest representable values.
    Memory benchmark: compare the current row loop against block sizes 4, 8, and 16.

Bottom line

The statistical core is largely well implemented. I would make these changes before calling it fully mathematically robust:

    clarify or enforce the distinction between adjacent crossings and global posterior boundaries;
    harden endpoint and extreme-tail calculations;
    prevent accepted Beta M-steps from decreasing the EM objective;
    describe the H step as empirical finite-sample stitching;
    replace the full one-hot initialization matrix;
    process sample rows in small transposed blocks rather than fully transposing the matrix.
