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
  default kind; nothing else can reach that state. (Closed the same day:
  see "The fit-subset draw also restores the RNG kind" below.)
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

## 2026-10-20 — Code review, chunk 4: package glue, data, README

The last chunk of the review (`high`): `R/betanorm-package.R`,
`R/data.R`, `DESCRIPTION`, and the README. Findings and dispositions are
in `dev/code-review.md`. All nine fixed; no package number moved, and no
data value changed.

- **The datasets credit their sources.** `horvath_goldstandard` is
  Horvath (2013)'s `probeAnnotation21kdatMethUsed.csv` (CC BY 2.0), and
  `GPL21145_sample` comes from NCBI GEO, platform GPL21145 (the user's
  account; a link is enough). Both are in the help pages' `@source`. The
  COPYRIGHTS catch-all excepts "the example datasets" without naming
  `data/` paths, because the downstream copies COPYRIGHTS.
- **`GPL21145_sample` keeps its 1062 zero-filled probes.** They are the
  gold-standard probes the EPIC array does not measure: none is in the
  EPIC annotation, and every other probe is. Removing them would move
  the bitcheck's real-data baselines, so the data stays as it is. The
  help page and the README say so and drop them before calibrating.
- **`horvath_goldstandard` is a plain data frame.** It was a readr
  `spec_tbl_df` with a `problems` pointer that is dead once loaded. Every
  column is identical, and so are the bitcheck's real-data cases.
- **`README.md` was first synced by hand**, as there is no pandoc here,
  with outputs from knitting `README.Rmd`. A pandoc render afterwards
  changed only one line wrap and two apostrophes. It no longer passes the
  removed `fit.policy`, and it no longer claims legacy numbers or
  presents `nL = 2` as the way to run EM longer.

## 2026-10-20 — dev/ is opt-in

`.gitignore` ignores everything under `dev/` and lists each tracked file.
A file is opted in only when CLAUDE.md, this log, `dev/to-do.md` or
`inst/COPYRIGHTS` cites it, or the gate runs it. Scratch work in `dev/`
stays out of git by default. Tracked:

- records: `DECISIONS.md`, `to-do.md`;
- the gate: `bitcheck-0.1.2.R`;
- reference copies (COPYRIGHTS, "Legacy compatibility"):
  `bmiq-calibration-legacy.R`; MEAT 1.24.0 cut to its calibration
  (`MEAT_1.24.0/MEAT/R/data_calibration.R`); and `preprocessCore/` cut
  to what the quantile normalizer derives from (`R/quantile_extensions.R`,
  `src/qnorm.c`, `src/qnorm.h`). Each copy keeps its `DESCRIPTION`, and
  MEAT its `LICENSE`, so the upstream licence notice stays with the code;
- prior art for the Boost bench in `to-do.md`: `boost_math_exports.cpp`,
  `inst/include/bmiqpp/boost_math.hpp`, `test-boost-math.R`.

Deleted: the finished plans (`0.1.2.md`, whose outcome is the 0.1.2
entry above, and `relicense-plan.md`, whose reasons are the next entry);
the July audit `opus_audit1.md`, superseded by the code review; the rest
of the MEAT 1.24.0 and preprocessCore copies, and a loose duplicate of
MEAT's calibration (`meat_data_calibration.R`); the pre-rewrite
`with_nl2/`; the July benchmarks and their output (`benchmark.R`,
`bench-qnorm.R`, `tmp/`); `Exported Items.bib`, whose three papers
DESCRIPTION and the help pages already cite; and
`slideimp_description.txt`. Older entries and `code-review.md` still cite
some of them. Each is in git history, for example
`git show 7110321:dev/with_nl2/bmiq-calibration.R`.

`code-review.md` is the exception to the rule: the chunk entries above
cite it, but it is kept locally, untracked (user, 2026-10-20). Those
entries hold the outcomes; the full log is in git history, as
`git show ae2f1ba:dev/code-review.md`.

## 2026-10-20 — `bmiq_calibration_fixed()`: the legacy configuration as one file (0.1.3)

Amended the same day by the gold quantile table entry below: the gold
side of `map_beta_q` now reads a table, and the output pin is 1e-14, not
`identical()`. Amended again by the "internal, quiet" entry below: the
function is not exported, has no `verbose`, and its warnings no longer
name samples or errors. The rest stands.

`R/bmiq_calibration_fixed.R` is `bmiq_calibration()` cut to the Horvath
legacy configuration, for the downstream package to vendor in place of
`R/bmiq_calibration.R` (user, 2026-10-20). It stands alone: it shares no
R code with the general file, and its helpers carry the `fixed_` prefix
so both files can load into one package (a second `fit_mixture` would
silently replace the first). The vendored-file code rules apply to it.

Fixed, not arguments: `nL = 3` with H, `niter = 5`, `nfit = 20000`,
`tol = 0.001`, `beta.maxit = 50`, `beta.score.tol = 1e-10`, the optional
H policy, a failed sample becomes an `NA` row (samples are independent,
so one cannot sink the run) and its index and error go in a warning. The
gold is chosen by name (`gold = "horvath"`) from `FIXED_GOLD`, the fitted
shapes, density crossings and class modes of `bmiq_gold_fit()` at its
defaults, written
by `data-raw/bmiq_gold_fixed.R` as hexadecimal floating-point literals.
Hex, because R's decimal parser is not correctly rounded on every
platform, so a 17-digit decimal could land one ulp off; a hex literal
parses exactly. The decimal stands beside each as a comment.

Cut: the nL = 2 path (`normalize_nl2`, the log-space maps), the raw gold
vector and the prefit object with their checks, `th1.v`, `doH`,
`h.policy`, `on.sample.error`, `failed.sample`, `debug` and every
diagnostic, the result object (a plain matrix comes back; the samples
that ran without H are named in the warning), `validate_em_settings`,
`check_gold_fit`, `check_thresholds`, the dimension checks in the
crossing search, and the responsibility-count diagnostic (so the
responsibility matrix is never reordered). Kept call for call: the fit
subset draw, the endpoint clip, the EM call, canonicalization by mean,
the density crossings, the tail split, `map_beta_q`, the H glue, the
±1e-12 clamp, and the gather/scatter blocks.

Pinned by `tests/testthat/test-bmiq-fixed.R` with `identical()`: the
stored Horvath numbers against the live `bmiq_gold_fit()`; the output
against `bmiq_calibration(datM, horvath_goldstandard$goldstandard2)` on
the six EPIC samples (20306 measured probes), on 11 synthetic samples
(two gather/scatter blocks), and on the skipped-H path (both fitters
mocked). The general function is untouched; the suite and the snapshots
pass unchanged. The function is not in `dev/bitcheck-0.1.2.R`, which
compares with 0.1.1, where it did not exist.

Known: the EM's `std::lgamma` is the platform libm (chunk 2 above), so
on another platform the live gold fit may differ from the stored numbers
at the last ulp, and the identity test would say so. The stored numbers
(fitted here, mingw gcc 14.2) are then the contract: with them, the gold
side of `bmiq_calibration_fixed()` is the same on every platform, which
the general function cannot promise.

Measured on the six EPIC samples (load_all at -O2, median of 3): general
with the raw gold 0.39 s, general with a prefit 0.32 s, fixed 0.35 s. The
gold fit is all that is saved; the quantile map still dominates
(`dev/to-do.md`, gold quantile table).

## 2026-10-20 — Gold quantile table in `bmiq_calibration_fixed()` (0.1.3, deliberate numerics change)

Amended the same day by the stored table entry below: the table is
built ahead of time into `R/sysdata.rda`, and the **API** paragraph
(the `qtable` argument, the exported builder, the passed-table checks)
is superseded. The table itself, the contract and the measurements
stand.

The gold side of the quantile map in `bmiq_calibration_fixed()` reads a
precomputed table in place of `stats::qbeta()`; the sample side stays
`stats::pbeta()` in probability space, and `bmiq_calibration()` is
untouched. The gold shapes are constant for a run (the gold is consumed
distributionally), so `qbeta(p, a.gold, b.gold)` is one fixed function
per component and tail, and Rmath's `qbeta` is a crude initial guess
plus Newton iterations on `pbeta`: a table is a better initial guess,
and at a fine enough grid it is the answer.

**The table** (`bmiq_gold_qtable(gold)`): per gold component, two
halves, Beta(a, b) and Beta(b, a), each the lower-tail quantile for p in
[1e-30, 1/2] on 8192 nodes uniform in u = log p, stored as per-interval
quintic Hermite coefficients of v = log x in t = (u - u_i) / h, from v,
dv/du = p / (x f(x)) and the analytic second derivative, evaluated by
Horner. The reflection x -> 1 - x swaps the shapes and the tails, so the
upper tail reads the swapped half and p above 1/2 reads 1 - p (as
`log(-expm1(u))`) from the swapped half. p = 0 and 1 are exact. p at or
below 1e-30, and any element the table does not place in [0, 1], fall
back to `stats::qbeta(p, ...)` in probability space, bit-identical to the
general function (no real input hits it). Nodes are anchored at
u1 = log(1/2), node i = u1 - (n - i) h, and the evaluator rebuilds a
node with that same expression: nodes from `seq()` reconstructed as
u0 + i h disagreed by 2.5e-11 at n = 4096, which put a 1e-12 floor under
every Hermite variant until found.

**Why quintic alone at n = 8192** (Phase 1, `dev/qtable-exp.R`,
2026-10-20). Reference: the production path `qbeta(pbeta(x))`. Real
inputs: 88,292 gold-side evaluations of the six EPIC samples at the
three call sites, reproduced bit for bit from `bmiq_calibration_fixed()`
first. Sweep: 81,015 p per (component, tail), all four combinations:
nodes, 20k interval interiors, 20k random, their mirrors above 1/2, the
ends, the neighbours of 1/2 and of p.min, p next to 1. Timings: sample 1
(14,778 evaluations), `load_all` at -O2, `bench::mark` median; "map" is
the full production map including the sample-side `pbeta`. mono =
monotonicity violations beyond 1 ulp over the four sweeps (Rmath's own
`qbeta` has 1).

| variant | n | real max | bit-id % | <=1 ulp % | sweep max | mono | us/elem | map ms | build ms | size KB |
|---|---|---|---|---|---|---|---|---|---|---|
| rmath |  | 0 | 100 | 100 | 0 | 1 | 2.07 | 34.7 | 0 | 0 |
| cubic | 1024 | 3.1e-08 | 0.0 | 0.0 | 3.1e-08 | 2 | 0.16 | 6.9 | 10 | 128 |
| cubic + Newton | 1024 | 6.3e-15 | 39.4 | 48.4 | 6.2e-15 | 2 | 0.58 | 13.1 | 10 | 128 |
| quintic | 1024 | 3.9e-11 | 2.8 | 3.4 | 3.8e-11 | 2 | 0.17 | 7.0 | 10 | 192 |
| quintic + Newton | 1024 | 3.3e-16 | 44.3 | 55.5 | 2.2e-16 | 2 | 0.58 | 13.1 | 10 | 192 |
| cubic | 2048 | 2.1e-09 | 0.0 | 0.1 | 2.1e-09 | 2 | 0.16 | 6.9 | 16 | 256 |
| cubic + Newton | 2048 | 2.2e-16 | 43.8 | 55.1 | 2.2e-16 | 1 | 0.57 | 13.0 | 16 | 256 |
| quintic | 2048 | 6.9e-13 | 11.6 | 14.4 | 6.9e-13 | 2 | 0.17 | 7.1 | 16 | 384 |
| quintic + Newton | 2048 | 3.3e-16 | 45.9 | 57.4 | 3.3e-16 | 1 | 0.58 | 13.1 | 16 | 384 |
| cubic | 4096 | 1.4e-10 | 0.2 | 0.2 | 1.4e-10 | 2 | 0.16 | 6.9 | 34 | 512 |
| cubic + Newton | 4096 | 3.3e-16 | 43.8 | 54.9 | 2.2e-16 | 2 | 0.58 | 13.0 | 34 | 512 |
| quintic | 4096 | 1.2e-14 | 27.9 | 35.5 | 1.2e-14 | 2 | 0.17 | 7.0 | 34 | 768 |
| quintic + Newton | 4096 | 3.3e-16 | 48.5 | 61.5 | 3.3e-16 | 2 | 0.58 | 13.1 | 34 | 768 |
| cubic | 8192 | 8.9e-12 | 0.8 | 1.0 | 8.9e-12 | 2 | 0.17 | 7.0 | 66 | 1024 |
| cubic + Newton | 8192 | 2.2e-16 | 44.0 | 55.2 | 3.3e-16 | 2 | 0.58 | 13.1 | 66 | 1024 |
| **quintic** | **8192** | **3.3e-16** | 39.8 | 51.4 | **3.3e-16** | 2 | **0.17** | **7.1** | 66 | 1536 |
| quintic + Newton | 8192 | 3.3e-16 | 51.1 | 65.1 | 2.2e-16 | 2 | 0.59 | 13.3 | 66 | 1536 |
| cubic | 16384 | 5.6e-13 | 3.5 | 4.3 | 5.6e-13 | 2 | 0.17 | 7.2 | 132 | 2048 |
| cubic + Newton | 16384 | 2.2e-16 | 44.5 | 55.8 | 2.2e-16 | 2 | 0.59 | 13.2 | 132 | 2048 |
| quintic | 16384 | 3.3e-16 | 40.4 | 52.0 | 2.2e-16 | 2 | 0.17 | 7.2 | 132 | 3072 |
| quintic + Newton | 16384 | 3.3e-16 | 51.2 | 65.5 | 2.2e-16 | 2 | 0.60 | 13.3 | 132 | 3072 |

Quintic alone converges as h^6 (3.9e-11, 6.9e-13, 1.2e-14, then the
3.3e-16 floor from n = 8192), where it is as accurate as any Newton
variant on the real inputs and the whole sweep; cubic alone is h^4 and
never reaches the bar by 16384. Any Newton step lands at the floor from
n = 2048 regardless of the interpolant, and costs 0.4 us per element
(map 13.0 ms against 7.1 ms per sample). Polishing the nodes with Newton
at build time changed nothing: the `qbeta(log.p = TRUE)` nodes already
solve log F = u to 1.4e-14. The two interpolants cost the same per
element; the lookup, not the arithmetic, is the cost. Newton is dropped,
not kept as an option: one path to test and vendor. Cost of the choice:
66 ms to build and 1.5 MB per table.

**API** (user, 2026-10-20): an explicit `qtable = NULL` argument after
`verbose` and the exported builder, not a package-level cache. Reuse is
visible in caller code, nothing lives in a hidden environment, and the
vendored file stays self-contained. `NULL` builds the table in the call;
a passed table is validated: layout, schema version, gold name, every
shape against `FIXED_GOLD`, the grid constants, finite coefficient
vectors of the right length, and a spot check at sixteen interval
midpoints per half against `stats::qbeta()` at 1e-12 (built tables sit at
1e-16; shifted, foreign or differently built coefficients miss by far
more; review, 2026-10-20). The version is bumped whenever what a table
holds changes, the coefficient construction included. A table saved with
`saveRDS()` and read back is accepted.

**Contract.** `bmiq_calibration_fixed()` agrees with the `calibrated`
element of `bmiq_calibration(datM, gold, on.sample.error = "continue")`
to within 1e-14 max-abs, with identical dimnames and NA pattern
(measured 3.3e-16 on the six EPIC samples, 41.8% of values
bit-identical). Not bit for bit, so a deliberate numerics change under
CLAUDE.md, as the single-exp E-step was: NEWS and this entry. The
snapshots at 1e-12 pass unchanged, since the general function is
untouched; the gold-fit pin stays `identical()`.
`tests/testthat/test-bmiq-fixed.R` sweeps the table against
`stats::qbeta()` for every gold set in `names(FIXED_GOLD)`, which is what
carries the guarantee to a gold set added later, and on the real
gold-side inputs harvested from the general function's diagnostics;
checks the nodes, the fallback, nineteen malformed tables, the saved
table round trip and the random state.

Measured on the six EPIC samples (20306 probes; installed at -O2, median
of 5, 10 ms timer): general 0.30 s, fixed with `qtable = NULL` 0.25 s,
fixed with a passed table 0.16 s, build alone 0.07 s. `devtools::check()`
clean.

Closed in `dev/to-do.md`: "Gold quantile table + Newton polish" (landed,
Newton dropped as measured above) and "Boost + RcppThread quantile-map
bench" (superseded: the table reaches most of that goal serially, with no
dependency and without forking the numerics beyond the last ulp; its
prior art, `boost_math_exports.cpp`, `inst/include/bmiqpp/boost_math.hpp`
and `test-boost-math.R`, stays tracked for a later threading thread).
`dev/qtable-exp.R` is opted in: the code comments, the tests and this
entry cite it. The hand-off plan `dev/bmiq_calibration_fixed.md` stays
untracked; this entry is its outcome.

## 2026-10-20 — The gold quantile table ships in `R/sysdata.rda` (0.1.3)

The table is built once, ahead of time, by `data-raw/bmiq_gold_qtable.R`
and stored as `FIXED_QTABLE` in `R/sysdata.rda`: one entry per gold set
in `FIXED_GOLD`, components `U` and `M`, each two half-tables with their
grid and coefficients. `bmiq_calibration_fixed()` reads it; nothing
builds a table at call time or at install (user, 2026-10-20). At
vendoring time the downstream pins a commit, copies `FIXED_QTABLE` from
this `R/sysdata.rda` into its own under the same name, and copies
`R/bmiq_calibration_fixed.R` stripping only roxygen lines: no edit after
the copy, so the two R files never differ.

Cut: the `qtable` argument, the exported `bmiq_gold_qtable()`,
`fixed_validate_qtable()` with its spot check, `FIXED_QTABLE_VERSION`
and the grid constants. The builder moved to the data-raw script
unchanged (same arithmetic, grid anchored at log(1/2)), so the vendored
file carries only the evaluator. Added: `fixed_gold_qtable()`, run on
every call, refuses a stored table whose component and half shapes are
not `identical()` to `FIXED_GOLD`. The table travels apart from the R
file, and one from another commit would otherwise map through the wrong
quantiles silently. It costs microseconds; the coefficients themselves
are checked by the tests, not per call.

Why ahead of time (measured 2026-10-20 on the six EPIC samples,
installed, `bench::mark` median of 15):

- **Per call** (the old `qtable = NULL` default): 65 ms per call to
  build, 183 ms against 117 ms for the six samples. It is per call, not
  per sample, so it dominates a caller that calibrates a few samples at
  a time.
- **At install** (one top-level line building the constant into the
  namespace): the same call time and nothing in the repo, but every
  installing machine rebuilds the table with its own Rmath, the
  downstream included. Rejected: the downstream takes the numbers, it
  does not rebuild them.
- **Hex literals in the R file**, as for `FIXED_GOLD`: 196,632 lines,
  4.6 MB of source in the vendored file. Rejected.

As with `FIXED_GOLD`, the stored table (built here, R 4.6.0, mingw gcc
14.2) is the contract on every platform, instead of each platform's
`qbeta`.

Costs: 1.08 MB in the repo and the tarball (xz), and another blob in git
history at each rebuild; `sysdata.rdb` is 1.09 MB installed. First use
in a session costs 23 ms, because R CMD INSTALL compresses a
`sysdata.rda` over 1e6 bytes with xz; `SysDataCompression: gzip` in
DESCRIPTION would make it about 3 ms at 1.28 MB, but it is not set,
since the downstream does not copy DESCRIPTION. `R/sysdata.rda` holds
`FIXED_QTABLE` and nothing else; the script overwrites it whole.

Output: `identical()` to the previous `qtable = NULL` path here on the
six EPIC samples and 11 synthetic ones, so the 1e-14 pin to
`bmiq_calibration()` and the snapshots are unchanged. Fixed 0.12 s,
general 0.29 s.

`tests/testthat/test-bmiq-fixed.R` checks the stored table against
`FIXED_GOLD` for every gold set and against the documented grid (8192
nodes from log(1/2) down to 1e-30, finite coefficient vectors of length
n - 1). The qbeta sweeps, the node test, the fallback and the real-input
test read the stored table, and a stored table with its components or
halves swapped, a shape one ulp off, or no entry for the gold stops the
call. The malformed passed-table cases and the saved-table round trip
went with the argument. To regenerate: `data-raw/bmiq_gold_fixed.R`
(paste its output), then `data-raw/bmiq_gold_qtable.R`.

## 2026-10-20 — `bmiq_calibration_fixed()` is internal, quiet, and plainly commented (0.1.3)

The downstream copies the file word for word, so any change it would
want is a commit here. The file therefore carries as little as it can
(user, 2026-10-20):

- **Internal.** Not exported: `@keywords internal` and `@noRd`, so no
  help page. The downstream calls its own copy, the tests reach it in the
  namespace, and a user can call `betanorm:::bmiq_calibration_fixed()`.
  The README no longer lists it.
- **No `verbose`.** The signature is `(datM, gold = "horvath")`. The
  per-sample messages (progress, H skipped, sample failed) are gone: for a
  few dozen samples the run takes seconds and they were noise. Above 200
  samples the run reports progress on one line, "Calibrating n samples:
  20%...40%...60%...80%...100%" (`message(appendLF = FALSE)` per step; a
  failed sample counts). The function is only used at the Horvath gold's ~21k
  probes (user), where a sample takes 20 ms here and about 50 ms on the
  user's 2021 machine, so 200 samples are 4 to 10 s of silence. A line per
  100 samples, and a fixed "may take a couple of minutes" notice, were
  considered: the notice was wrong at that scale, and one line is all the
  progress a run of minutes needs. Steps by count, not by time, so the
  output is the same on every machine.
- **Plain conditions.** Each `stop()` and `warning()` is one fixed
  message, with no sample index, no values and no `call. = FALSE`, and
  the `context` argument the helpers threaded through for their messages
  is gone. A sample that cannot be calibrated is still an `NA` row, which
  is the signal. One warning says that some samples are `NA`, another that
  some were calibrated without H (which the output does not show). The
  sample lists in those warnings went, and with them their cap at ten
  (`fixed_list_samples()`, from the review of the same day).
  `fixed_gold_qtable()` stays; its message is now "FIXED_QTABLE does not
  match FIXED_GOLD."
- **Comments cut to what the code needs** (the licence header, the
  component order, the table layout, the node placement, the H map, the
  block loop). The reasons they used to give are recorded here or in the
  entries named:
  - the endpoint clip is half the distance to the nearest interior value,
    as in legacy BMIQ; the length guards cover an all-0 or all-1 fit
    subset ("Validation: one owner per invariant"), and the floor at
    `sqrt(.Machine$double.eps)` keeps a subnormal value from rounding
    back to 0 or 1 and feeding `log(0)` to the EM;
  - the sample modes are clamped to [0, 1] because `density()` is not
    boundary-corrected, so its Gaussian tails can put the mode just
    outside; the clamp only moves the initial thresholds toward a valid
    value;
  - the mode-shifted initial thresholds are EM starting cuts, not fitted
    roots, so sorting a crossed pair is safe;
  - two Beta densities can cross twice, so the crossing taken is the one
    where component k stops dominating and k + 1 takes over; the
    log-density ratio has at most one interior turning point, so cutting
    there brackets every sign change;
  - one fit subset is drawn for all samples (same probe count, `nfit` and
    seed);
  - every U and M value goes to exactly one tail, so a value at the
    component mean is mapped too;
  - the H set is never empty (code review, chunk 3), and the H map is
    glued to the calibrated U maximum and upper-M minimum without the
    legacy gaps ("Legacy compatibility", step 3);
  - `datM` is asserted double, not coerced ("Validation: one owner per
    invariant"); the output is a fresh matrix, so the scatter writes into
    it without touching `datM`; a block of eight samples is one 64-byte
    cache line per probe ("Rejected cleanups");
  - the table's layout, fallback and node placement, and the 1e-14
    measurement: the gold quantile table entry above;
  - the hex literals: the one-file entry above.

Output `identical()` to before on the six EPIC samples and 11 synthetic
ones. `test-bmiq-fixed.R` drops the message and sample-list checks, pins
the two warnings' text and the progress line, which only a run of more
than 200 samples gets (201 samples that fail at once), and mocks the
fitter with its new signature.

## 2026-10-20 — The fit-subset draw also restores the RNG kind with no seed (0.1.3)

In both `draw_fit_indices()` and `fixed_draw_fit_indices()` (user,
2026-10-20). The draw keeps its save-and-restore rather than the plain
`set.seed(1)` of legacy BMIQ (`dev/bmiq-calibration-legacy.R`): a fixed
seed left in place replaces the caller's stream, so every draw after
the call, in every iteration of a loop, would follow seed 1.

The gap closed: with no `.Random.seed`, the caller's kind lives only in
R's internal state, and the draw's `set.seed(kind = "Mersenne-Twister",
...)` changes it. Removing the seed afterwards left the session on
Mersenne-Twister / Inversion / Rejection. The no-seed branch now saves
`RNGkind()` (which, called with no arguments, writes no seed) and sets
it back on exit, then removes the seed that call writes. Reachable only
by choosing a kind and then removing `.Random.seed` by hand, since
`RNGkind()` and `set.seed()` both write one. Pinned in
`test-random-state.R` for `bmiq_calibration()`, `bmiq_gold_fit()` and
`bmiq_calibration_fixed()` under two non-default kinds (6 of 6 checks
fail without the fix). No number moves: the draw itself is unchanged.

Checked by a probe against a run without the call: a fresh session, a
`set.seed(42)` session with or without draws since, "L'Ecuyer-CMRG" and
`sample.kind = "Rounding"` all keep the same seed, kind and next draws,
and the call leaves no other object in the global environment.

Known, and not fixable from R: under `normal.kind = "Box-Muller"`, R
keeps the second normal of each pair in an internal variable that
`.Random.seed` does not hold, and every `set.seed()` clears it. After an
odd number of `rnorm()` draws, the caller's next normals differ from a
session without the call. Only a draw that bypasses R's generator (the
subset reproduced in C++) would avoid it, which is not worth the code
for a normal kind few sessions use.

## 2026-10-20 — A kernel-only review of `bmiq_norm.cpp` (0.1.3)

An agent shown `src/bmiq_norm.cpp` alone, without the R callers or this
log, reviewed it for CRAN and memory safety and for speed. What it
raised that bears on this package, and the outcome of each, is recorded
here; the report itself is not tracked (user, 2026-10-20). It found no
memory-management fault: the kernels write only through live Rcpp
objects, nothing needs `PROTECT` by hand, and the copy of
`initial_responsibility` keeps the E-step off the caller's matrix. One
test was added; no code changed and no number moved.

Taken:

- **A test that `datM` comes back as found.** The report's main point:
  `scatter_sample_block_cpp()` writes into its destination in place,
  past copy-on-modify, so it is safe only while every caller passes a
  freshly allocated private matrix, and a test should hold the public
  functions to leaving their input unchanged, compared with an
  independent copy rather than another binding. Both functions allocate
  a fresh output, but nothing pinned it: an edit such as
  `calibrated <- datM`, to keep the dimnames, would write the calibrated
  values into the caller's matrix. `test-input-preserved.R` runs
  `bmiq_calibration_fixed()` and `bmiq_calibration()` (`failed.sample`
  "NA" and "original") over 11 samples, two blocks with a failed sample
  in the second, and compares `datM` with a copy made by `serialize()`
  (`saved <- datM` would be a second binding of the same vector). With
  `calibrated <- datM` in either file, it fails. The rest of the point
  (the wrapper stays internal; a `MAYBE_SHARED()` check would refuse the
  legitimate caller) is already so, and in the kernel's comment.
- **When a CRAN submission is planned:** one run under valgrind and the
  sanitizers, and a small case under `gctorture(TRUE)`. Only
  `devtools::check()` has been run so far.

Declined:

- **Interrupt checks inside the passes.** The report: one
  `checkUserInterrupt()` per EM iteration leaves the scan, the gather
  and scatter and the E-step without one; check every 65536 elements.
  No call runs long enough: an EM call takes 10 ms with a check at each
  iteration, about 2 ms apart, and gathering six samples takes 39 µs
  (measured below).
- **Validating the EM controls and initial responsibilities in C++.**
  The report: `maxiter` NA or below 1 returns placeholder shapes,
  `BetaStats::add()` skips negative and NaN responsibilities silently,
  and the rows need not sum to 1. "Not taken: C++-side validation of the
  EM scalars" (2026-10-05) stands. The fixed path passes constants and the kernel defaults; the
  general path asserts `niter >= 1` and passes the defaults for the
  halving limit, `min_shape` and `armijo`. Both build `w.init` one-hot
  with at least two members per class, so each component starts with an
  effective size of at least 2, and the E-step's own responsibilities
  are finite and non-negative by construction. Even `maxiter = 0`, which
  returns placeholder shapes with `eta = 0`, would stop at the fixed
  path's "Invalid mixture parameters." check.
- **Subtraction-based block bounds.** The report: `first0 + count >
  n_samples` can overflow; test `count > n_samples - first0` instead.
  It can overflow only where `R_xlen_t` is 32 bits; R dropped 32-bit
  Windows in 4.2.0 and CRAN no longer checks i386. Callers pass a start in 1..nrow
  and a count of at most 8, and chunk 2 does not reopen the 2026-10-05
  hardening. The comments' "wraps negative" is loose (signed overflow is
  undefined; it wrapped in practice); the next edit of the file can say
  so.
- **The halving loop's counter.** The report:
  `for (int h = 0; h <= max_halving; ++h)` overflows `h` at
  `max_halving = INT_MAX` when no step is accepted. The loop runs to its
  limit only when a shape starts exactly at
  `min_shape` and Newton points below it, and neither caller passes
  `beta_max_halving`, so the limit is the default 30. The 2026-10-05
  rule applies: revisit if the argument becomes reachable from user
  input.
- **Tolerance-based instead of exact regression tests.** The report:
  `std::lgamma`, `log` and `exp` differ across platforms. Against
  "Numeric exactness is the contract". The snapshots already hold to
  1e-12, not bit for bit, and the portability point is chunk 2's
  `std::lgamma` entry.
- Agreed with what is here: the extra memory is ordinary (about 40
  bytes per observation at K = 3, the cached logs and the copied
  responsibilities), and the numerical safeguards stay (Welford
  variance, log-sum-exp, the cached logs, no `-ffast-math`).

**Speed.** The report proposed fewer full passes over the data and warm
starts for Newton, without measuring. Measured with R CMD build and INSTALL at R's
default -O2 (R 4.6.0, `bench::mark` medians) on the six EPIC samples
(20306 probes), the EM calls on the inputs captured from the fixed path
(n = 20000): `bmiq_calibration_fixed()` 113.7 ms, 19 ms per sample; one
EM call 9.99 ms at `maxiter = 5` (every call runs all five), 2.84 ms at
`maxiter = 1`, and 0.55 ms on a 600-point subset of the same inputs. The
EM is about 53% of the fixed path: since the gold quantile table, chunk
2's "the quantile map is where the time goes" holds for
`bmiq_calibration()` but not for the fixed path. The report's proposals
still miss:

- Warm-starting each M-step's Newton fit from the previous shapes
  instead of the method of moments: the per-component work is at most
  0.55 ms per call, about 3% of the run, and the change moves numbers.
- Accumulating the next M-step's statistics inside the E-step loop, so
  the responsibilities are not read again: bit-identical, since each
  component gets the same values in the same order, but it saves only
  rereading a working set of about 1 MB that sits in cache; the Welford
  divisions and the E-step's `exp` and `log` stay. A bench candidate if
  the EM ever needs to be faster.
- Making `w` optional changes the return list pinned by
  `test-cpp-contract.R`, and the general path reads `w` for its
  responsibility-count diagnostic. Fixed-size arrays and templates on
  K, and building the status strings once at return, save microseconds
  per call. Tiling the gather and scatter: the gather takes 39 µs and
  is already one cache line per probe ("Rejected cleanups").

Noted, not measured: per observation the E-step calls `exp` K times and
`log` once, likely most of the EM's time. For the largest component
`log_component[k] - maximum` is exactly 0 and `exp(0)` is exactly 1, so
skipping that call would be bit-identical.

Known, for a CRAN submission: testthat skips snapshot expectations on
CRAN (`cran = FALSE`), so `test-numeric-regression.R` does not run
there, but the `identical()` pin of the stored gold fit against a live
`bmiq_gold_fit()` in `test-bmiq-fixed.R` does, and on another libm it is
expected to fail (the one-file entry above). It would need
`skip_on_cran()`.

## 2026-09-25 — GPL (>= 3), with the upstream authors as `cph` (0.1.0)

Recorded 2026-10-20, from the deleted `dev/relicense-plan.md`.

- betanorm was written with preprocessCore (LGPL >= 2) and RPMM (GPL >= 2)
  open as the reference, and the BMIQ calibration with Teschendorff's
  BMIQ and Horvath's BMIQcalibration. Code translated from those sources
  is derived from them. LGPL allows a move to GPL and GPL >= 2 does not
  allow MIT, so the MIT licence was wrong as declared. The user accepted
  GPL; a clean-room rewrite to keep MIT was not wanted.
- `GPL (>= 3)`, not GPL-2 alone: it covers LGPL >= 2, GPL >= 2, and the
  GPL-3 copies of BMIQ in ChAMP and wateRmelon.
- The upstream authors are `cph` only, not `aut` or `ctb`: they did not
  write this package. A paper's authors are cited, not credited as `cph`,
  unless their code was the reference.
- No copy of the GPL text: R ships it. Each derived file names its origin
  in a header comment, and `inst/COPYRIGHTS` keeps the per-file record.

## Legacy compatibility (predates this log; do not drift)

- `nL = 3`, `niter = 5` is the intentional legacy BMIQ configuration
  (Horvath-clock-style pipelines). The cap is a compatibility choice, not a
  convergence claim; `LEGACY_BMIQ_NITER` is the single source for both the
  default and the drift warning.
- C++ returns keep the old RcppArmadillo shapes (K x 1 matrices, bare `eta`
  vector) — pinned by `test-cpp-contract.R`.
- **The configuration is legacy, the numbers are not.** The rewrite
  (`f536c1d`, 2026-07-11) departs from legacy BMIQ code
  (`dev/bmiq-calibration-legacy.R`, and MEAT 1.24.0's
  `dev/MEAT_1.24.0/MEAT/R/data_calibration.R`) on purpose, in three steps.
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
