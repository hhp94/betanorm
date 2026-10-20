# Design decisions

Settled decisions with the evidence that settled them. Don't re-litigate
without new evidence; append new decisions with a date.

## 2026-08-10 — Gold standard is a reusable prefit, not a per-call vector

The gold standard is consumed purely distributionally: the sample path uses
only the fitted component shapes, the density-crossing thresholds, and the
two class density modes. Nothing indexes the gold vector per probe.
Therefore `bmiq_gold_fit()` fits it once and `bmiq_calibration()` accepts
either a raw vector or the prefit object; the old `length(gold) == ncol(datM)`
assertion enforced an alignment that never existed and was removed.
Prefit-vs-vector paths verified identical (`test-gold-fit.R`).

## 2026-08-10 — Serial execution; no in-process threading

Profiled (450k probes x 4 samples, and 21k x 40):

- 450k scale: `qbeta` 76% self / 86% total, `pbeta` 11%, compiled EM ~4%.
- 21k (clock-style) scale: 0.09 s/sample warm -> 1.5 min per 1000 samples,
  3 min per 2000. That is the realistic cohort ceiling and it is fine serial.

An OpenMP port of the quantile map measured ~5-7x (bit-identical vs R's
vectorized `qbeta(pbeta(...))`, 1.75 s -> 0.24 s at 450k) but was rejected:

1. **Rmath thread hazard.** `Rf_qbeta`/`Rf_pbeta` are pure on the happy
   path, but their error paths route `ML_WARN -> Rf_warning` — R API from a
   worker thread is undefined behavior. Guaranteeing no edge input ever
   reaches a thread is a fragile invariant.
2. **macOS ships no OpenMP** (Apple clang); CRAN builds strip `-fopenmp`.
3. The thread-safe alternative (RcppThread + Boost `ibeta_inv`, see
   `dev/test-boost-math.R`) is **not bit-identical to Rmath**, so it would
   fork the numerics and break the snapshot contract; it could only ever be
   an opt-in mode.

Serial also keeps ordered progress messages, simple
`on.sample.error = "continue"` semantics, and deterministic warning order.

## 2026-08-10 — Parallelism, if ever needed, is process-level and user-side

The prefit gold object (tiny) makes chunked process parallelism trivial:
split `datM` into row chunks **before** shipping, `parLapply` over chunks so
each PSOCK worker serializes only its slice. Never export the full matrix.

Caveats recorded:

- Linux fork CoW is not free: the GC writes mark bits into SEXP headers in
  each child and refcount updates dirty pages on mere reads, so inherited
  heaps creep toward copies; child-side `datM[rows, ]` materializes fully
  regardless.
- ALTREP shared memory (Bioconductor `SharedObject`) is the genuine
  zero-copy option (works over PSOCK and on Windows; worker GC only touches
  a stub). Crossover is ~EPIC x 2000-sample scale (>10 GB). Footguns: any
  R-level slicing silently materializes; shm lifecycle management; heavy
  dependency. If adopted, the integration point is a worker entry that
  passes the shared matrix + row range straight into the existing
  `gather_sample_block_cpp` machinery (reads via `DATAPTR`, zero-copy).

## 2026-08-10 — No rcpptimer instrumentation (yet)

Profiling shows the compiled layer is ~4% of runtime at array scale;
instrumenting it would measure the wrong side. If EM internals ever need
tuning, wrap rcpptimer behind a `LOC_TIMER` compile flag (header design
already specified: `LOC_TIMER_OBJ`/`LOC_TIC`/`LOC_TOC` macros compiled out
by default) so release builds carry zero overhead.

## 2026-08-10 — Rejected cleanups (evidence-based)

- **`findInterval` for `class_by_thresh`**: tie-semantics risk vs snapshots;
  the loop runs at most twice (nL <= 3) and is already vectorized.
- **Single-`exp` E-step**: originally rejected here on an *assumed* ULP
  break. Superseded — measured and adopted 2026-08-10; see the single-exp
  entry below.
- **Gather/scatter block machinery**: deliberate cache-layout optimization
  (samples are strided rows in a column-major matrix), documented and
  contract-tested; a flat `datM[ii, ]` loop is simpler but was kept out.
- **Replacing `scan_finite_unit_interval_cpp` with checkmate**: the scanner
  is single-pass, copy-free, and contract-tested; not worth the churn.
- **Merging the nL = 2 and nL = 3 normalization paths**: genuinely different
  math (truncated conditional maps vs plain quantile maps + H stitching);
  merging would change legacy nL = 3 output.

## 2026-08-10 — Validation: one owner per invariant

Defensiveness audit outcome. Each invariant has exactly one owner; a second
check of the same invariant elsewhere is a bug (drift risk, misleading error
source), not extra safety:

- **Type/shape of user args** — checkmate at the two exported entry points.
  `datM` must already be storage-mode double (asserted, not coerced): an
  integer matrix reaching the C++ layer would silently coerce-copy per block
  call, so rejection is the honest contract.
- **Finite + [0, 1] range (incl. NA/NaN)** — `scan_finite_unit_interval_cpp`,
  single pass. checkmate's `any.missing` was dropped as a duplicate owner.
- **Sanity of the fitted mixture** — `canonicalize_em_components()` at the
  C++ boundary. `density_thresholds()` keeps only what the crossing search
  itself needs (dimensions + strictly increasing means).
- **EM's strict-interior (0, 1) precondition on y** — the C++ EM itself.
- Guards that look defensive but are reachable and stay: the clipping-range
  guard (all-0/1 fit subsets), `th2.initial` collapse on degenerate data,
  the +/-1e-12 output tolerance + clamp (H map ulp overshoot), and all
  statistical class-count / crossing checks.

## 2026-08-10 — Single-exp E-step adopted (deliberate numerics change)

Reverses the same-day rejection above, which assumed rather than measured
the ULP impact. The E-step now reuses the max-shifted exponentials for the
responsibilities — `w = exp(lc - max) * (1 / sum_exp)` — instead of
re-exponentiating `exp(lc - (max + log(sum_exp)))`. Measured on a verbatim
EM clone over an 80-config grid (nL 2/3 x n 5k/20k x maxiter 5/25 x
10 seeds; `scratchpad/estep_variants.cpp` bench):

- A single E-step pass differs by at most 1 ulp (2.2e-16 on
  responsibilities in [0, 1]).
- Through EM feedback: relative drift <= 8e-8 on shapes, <= 2e-9 on the
  log-likelihood; zero outer iteration-count or convergence flips across
  160 config-mode pairs.
- Speed: E-step kernel 4.6 -> 2.7 ms per pass at n = 20k (1.7x); a full
  EM call (niter = 5, nL = 3, n = 20k) 26 -> 16 ms. At clock scale the EM
  is roughly a third of the ~0.09 s/sample budget, so this saves
  ~10 ms/sample (~11% end-to-end); at 450k-probe scale the EM is only ~4%,
  so the gain there is ~2%.

Snapshots were regenerated in the same commit. Observed snapshot deltas:
calibrated output moved <= 7e-10; borderline diagnostics can flip
(inner-Newton `fit_status` converged <-> max_iter at the score_tol
boundary; the near-degenerate nL = 2 EM fixture now converges 2 iterations
earlier, `converged` FALSE -> TRUE). The 1e-12 exactness contract now pins
the single-exp numerics.

## 2026-08-10 — Profile with an -O2 DLL; load_all() builds at -O0 by default

pkgbuild (used by `devtools::load_all()`) injects `-g -O0` debug flags
unless `options(pkg.build_extra_flags = FALSE)` is set. R's own compiled
code (`qbeta`, `pbeta`) is always optimized, so profiles taken against a
load_all DLL overstate the package-C++ share. Corrected 21k stage
attribution with an -O2 DLL (40 samples, post single-exp):

- 72 ms/sample total (was 85.5 measured against the -O0 DLL);
- `map_beta_q` (qbeta+pbeta) **67.6%**, 3.65 us per probe eval,
  ~12.9k evals/sample;
- `fit_mixture` (EM + classing) **26.3%** (EM call itself 15.7 ms,
  matching the -O2 sourceCpp clone used in the single-exp bench);
- everything else 6.1%.

Earlier entries' *shares* for the compiled layer ("EM ~a third at 21k",
"~4% at 450k") were measured against -O0 DLLs and overstate the EM;
their conclusions stand (qbeta dominates even more, not less). Any future
profile or bench must set `options(pkg.build_extra_flags = FALSE)` and
rebuild (`pkgbuild::clean_dll()`) first, or compare against
`Rcpp::sourceCpp` output (which uses the standard -O2 Makeconf flags).

## 2026-10-05 — Kernel hardening taken back from methylCIPHERv2 (partly)

methylCIPHERv2 vendors these kernels and on 2026-08-10 ran a memory-safety
audit over them (its `dev/DECISIONS.md`, "The kernels bounds-check their own
arguments"). Taken back here, all bit-identical, none reachable from the
exported functions:

1. `gather_/scatter_sample_block_cpp`: block arithmetic in `R_xlen_t` with
   `NA_INTEGER` rejected first. In `int`, a large or `NA` `first_sample`
   wrapped negative, passed the existing bounds test, and segfaulted (the
   scatter's was a heap write).
2. `scatter_sample_block_cpp`: `destination` is `SEXP`, checked as a double
   matrix before wrapping. As `NumericMatrix`, a non-double matrix was
   coerced into a temporary, scattered into, and dropped silently.
3. `qnorm_target_rows_cpp`: comparator is a total order with NaN last. The
   bare `<` was UB in `std::sort` with NaN present.
4. `target_unequal_rank`: the two early returns clamp the index like the
   interpolating branch. Unreachable with ranks in [1, n]; kept for symmetry.

These fix guards that already existed (1, 3, 4) or a parameter type (2), so
they add no second validation owner.

**Not taken: C++-side validation of the EM scalars** (`maxiter`,
`beta_max_halving`, `min_shape`, `tol`, `armijo` NA/range checks). It adds
a second owner for arguments R already validates or defaults, against
"validation lives in R; C++ trusts its caller". Revisit only if those
arguments become reachable from user input.

The other methylCIPHERv2 divergences (no RNG/`nfit`, prefit-only gold, no
debug/verbose, per-sample value scan, no warnings) are downstream choices
for its own front door and stay there.

## 2026-10-20 — The vendored files follow the downstream code rules (0.1.2)

A downstream package vendors `R/bmiq_calibration.R`, `R/quantile_norm.R`,
`src/bmiq_norm.cpp`, `src/quantile_norm.cpp` and `inst/COPYRIGHTS`, and its
copy script strips roxygen lines and nothing else. Its code rules are
therefore kept here, so the files copy with no further edit (CLAUDE.md
invariants). One commit each:

1. **`rng = false` on every Rcpp export.** The default wrapper's
   `RNGScope` calls `GetRNGstate()`, which creates a time-seeded
   `.Random.seed` in a session that has none. So `bmiq_calibration()`,
   `bmiq_gold_fit()` and `quantile_norm()` each left one behind, although
   `draw_fit_indices()` restores the caller's state around its own draw.
   No kernel draws. Pinned by `test-random-state.R`.
2. **`[[`, never `$`, in R/.** `$` partial-matches on lists. Before the
   switch, the suite and the real-data calls ran under
   `warnPartialMatchDollar = TRUE` with no partial match, so `[[` reads
   the same fields. (`[[<-` on `NULL` builds an atomic vector where `$<-`
   builds a list; every `diagnostic[[...]] <-` runs on a list.)
3. **Exact matching, no `match.arg()`.** `h.policy`, `on.sample.error` and
   `failed.sample` default to one string, checked by
   `checkmate::assert_choice()`. The one user-visible change: an
   abbreviation such as `"opt"` is now an error
   (`test-option-matching.R`).
4. **Calls outside base are qualified**, so the vendored code does not lean
   on this package's `importFrom` lines. Those for stats and utils are gone
   (`dbeta` was imported with no caller).
5. **No Makevars.** They held the RcppArmadillo skeleton's flags
   (`-fopenmp`, LAPACK, BLAS, Fortran libraries), left over from before
   the port to plain Rcpp. Nothing in `src/` uses them.

Evidence (`dev/bitcheck-0.1.2.R`): it installs with R CMD INSTALL, runs 24
cases (real and synthetic calibration, raw gold and prefit, the H-failure
option paths, the EM kernel, `quantile_norm()`) plus two runs on 0.1.1's
saved gold fits, and compares with 0.1.1 by `identical(num.eq = FALSE)`.
After every step all 26 values were identical; the only text change is two
warnings, which now name `result[["failures"]]` and `result[["h.applied"]]`.
After the Makevars went, the same held with R's default -O2 against 0.1.1
built the same way, and 0.1.1 itself is identical at -O2 and -O3.
Snapshots unchanged throughout.

## 2026-10-20 — Code review, chunk 1: quantile_norm

The first `/code-review` (`high`) of `R/quantile_norm.R` and
`src/quantile_norm.cpp`; findings and dispositions in `dev/code-review.md`.
Changed, with every 0.1.1 number kept (bitcheck 26/26, snapshots
unchanged):

- **`obj` finiteness has an owner.** The docs required it but nothing
  checked it, so `Inf` was ranked and mapped to the ends of `target`. One
  `checkmate::assert_numeric()` call now owns missing and finite for `obj`,
  as one does for `target`; `assert_matrix()` owns the shape alone.
- **One column maps every value to `min(target)`: documented, not
  changed.** With `ncol(obj) == 1` a row has nothing to rank. Mapping to
  the median, or refusing the input, would move numbers computed with
  0.1.1.
- Ties, the one-column case, dimnames and the input checks are pinned in
  `test-quantile-norm.R`; the snapshots never reached them.

Kept as they are (a later review that raises them needs new evidence):

- **The two full transposes** (peak memory about 3x `obj`). A strided
  gather and scatter would avoid them, but `quantile_norm()` is a one-shot
  call on the user's matrix (nothing in the package calls it), and the gain
  is unmeasured. Revisit if memory is shown to bind.
- **The dead `fraction == 1.0` branch and the redundant `floor()`** in
  `target_unequal_rank()`. Kept so the port still reads line by line
  against preprocessCore.
- **The absolute `4 * DBL_EPSILON` snap.** It mirrors preprocessCore. At
  large `m` it stops snapping, and the interpolation it does instead is
  off by about 1e-11 times the gap between neighbouring target values. A
  relative tolerance would move numbers.

## 2026-10-20 — Code review, chunk 2: bmiq_norm.cpp

The `/code-review` (`high`) of `src/bmiq_norm.cpp`; findings and
dispositions in `dev/code-review.md`. None of the 2026-10-05 hardening is
reopened.

Changed, no number moved:

- **A retained component keeps its real fit status.** When the ascent
  guard keeps a component's previous shapes, `fit_status` was set to
  "converged" whatever that fit had ended as; it now stays the status of
  the fit the shapes came from. The nL = 2 EM snapshot takes this path at
  its last iteration ("max_iter" twice, was "converged"), so that one
  snapshot was redone. Against 0.1.1, every other field of that case and
  every other bitcheck case are bitwise identical.

Adjudicated against "Validation: one owner per invariant" (2026-08-10):

- **The kernel's `nL` in {2, 3} check stays.** It duplicates R's
  `checkmate::assert_int(nL, lower = 2L, upper = 3L)`, and K = 4 would be
  memory-safe (the responsibility shape check is the bounds guard). But no
  R code downstream of the kernel handles more than three components, and
  `test-cpp-contract.R` pins the check as the kernel's own contract. It is
  an exception to the rule, now listed in CLAUDE.md beside the other
  deliberate C++ checks.

Known and accepted:

- **`log_beta_norm()` uses `std::lgamma`, from the platform libm**, while
  digamma and trigamma come from Rmath. On this machine (mingw, gcc 14.2)
  `std::lgamma` differs from R's `lgamma` at 58% of 200k points in
  (0.001, 100), by up to 2.8e-14 (2e-11 relative near its zeros at 1 and
  2), so another libm gives slightly different EM numbers. Whether they stay
  inside the 1e-12 snapshots is unmeasured (no CI on another platform).
  Moving to `R::lgammafn` or `R::lbeta` would make the EM portable but is a
  deliberate numerics change, like the single-exp E-step, and would break
  results computed with 0.1.1. The cancellation in the three-term form does
  not matter at real shapes (1.5 to 38 in the bitcheck cases: about 6e-14
  from `lbeta`). If EM threading is ever tried: glibc's `lgamma` writes the
  global `signgam`, so it is not thread-safe.

Kept as they are:

- **The non-finite-concentration fallback to 1.0** in
  `fit_beta_from_stats()`. Reachable only with a subnormal weighted
  variance; a zero variance is refused one check earlier.
- **The Welford update's division per observation.** Plain weighted sums
  would change the arithmetic, and the EM is about 26% of runtime (the -O2
  profile); the quantile map is where the time goes.
- **The redundant `!std::isfinite(p)`.** Bit-identical either way; not
  worth churn in a vendored file.

Observed, not acted on: the Newton fit often stops at `beta_maxit = 50`
("max_iter") rather than at the 1e-10 scaled-score tolerance (14 of 38 EM
iterations in the nL = 2 snapshot case, for one component or both), with
every step still accepted, i.e. at the noise floor. The cost is per
component, not per observation.

## 2026-10-20 — Code review, chunk 3: bmiq_calibration.R

The `/code-review` (`high`) of `R/bmiq_calibration.R`; findings and
dispositions in `dev/code-review.md`. Every fix went through the gate:
suite passes, snapshots as after chunk 2, bitcheck values as after chunk 2
(only `em_nL2`'s `fit_status` differs from 0.1.1, from that fix).

Changed, no number moved in a default session:

- **The fit subset is drawn under R's default generator.**
  `draw_fit_indices()` seeded the caller's `RNGkind()`, so the gold and
  sample subsets, and every calibrated value, followed it: up to 0.021
  under "L'Ecuyer-CMRG" and 0.0056 under `sample.kind = "Rounding"` on a
  synthetic sample. Legacy code had the same dependence. `set.seed()` now
  names Mersenne-Twister / Inversion / Rejection, and the restore already
  in place puts the caller's kind back (the kind lives in `.Random.seed`).
  Pinned in `test-random-state.R`. A session that removed its
  `.Random.seed` by hand after choosing another kind is left on the
  default kind; nothing else can reach that state.
- **A prefitted gold has an owner.** `check_gold_fit()` checks a
  user-supplied `bmiq_gold_fit` (`nL`, `a`, `b`, `thresholds`, the two
  modes) at the entry point. A fresh `bmiq_gold_fit()` result satisfies it
  by construction, so the raw-vector path does not run it and there is
  still one owner per invariant.
- **The empty-H check is gone** (with the no-op `unique()`). It was a
  second owner of "class 2 is not empty", which `require_all_classes()`
  already owns. "Validation: one owner per invariant" keeps the
  *reachable* class-count checks; this one was not reachable.
- **The gold is fitted before the output is allocated**, so a bad gold
  fails before an n x p matrix exists.
- Docs: the defaults are a legacy configuration, not legacy numbers (see
  "Legacy compatibility" below); a skipped H step also leaves the H class
  unchanged; a prefitted gold ignores the call's EM settings.

Kept as they are:

- **`result[["call"]]` is `match.call()`.** Under `do.call()` it holds the
  argument values (1.79 MB for a 24 KB `datM`), but a plain call stores
  symbols, nothing reads the field, and removing it is a visible change.
- **Each sample's debug diagnostics hold the shared fit indices.** One
  vector in memory; only `saveRDS()` repeats it per sample. Moving it
  would change the diagnostics shape for a debug-only cost.

## Legacy compatibility (predates this log; do not drift)

- `nL = 3`, `niter = 5` is the intentional legacy BMIQ configuration
  (Horvath-clock-style pipelines). The cap is a compatibility choice, not a
  convergence claim; `LEGACY_BMIQ_NITER` is the single source for both the
  default and the drift warning.
- C++ returns keep the old RcppArmadillo shapes (K x 1 matrices, bare `eta`
  vector) — pinned by `test-cpp-contract.R`.
- **The configuration is legacy, the numbers are not.** The rewrite
  (`f536c1d`, 2026-07-11) departs from legacy BMIQ code
  (`dev/bmiq-calibration-legacy.R`, MEAT) on purpose, in three steps.
  Recorded 2026-10-20 (code review, chunk 3); the roxygen Details say the
  same.
  1. Class thresholds are weighted-density crossings
     (`density_thresholds()`). Legacy code takes midpoints between
     neighbouring classes' extreme values in the fit subset.
  2. Components are fitted by Newton's method under the ascent guard.
     Legacy code uses BFGS (`blc`) or Nelder-Mead (Horvath's `blc2`).
  3. The H map is glued straight to the calibrated U maximum and upper-M
     minimum. Legacy code offsets both anchors by the raw gaps (`deltaUH`,
     `deltaHM`), so its output jumps at each join. Here the map has no
     jump; the cost is that the outermost H probes tie with the probes
     across each join (3000 distinct inputs, 2998 distinct outputs, in
     the review check).
  Changing any of them moves numbers: a feature change, not a cleanup.
