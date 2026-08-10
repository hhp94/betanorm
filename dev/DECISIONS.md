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
- **Single-`exp` E-step** (`exp(a)/s` vs `exp(a - log s)`): last-ULP change,
  breaks bit-exact snapshots.
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

## Legacy compatibility (predates this log; do not drift)

- `nL = 3`, `niter = 5` is the intentional legacy BMIQ configuration
  (Horvath-clock-style pipelines). The cap is a compatibility choice, not a
  convergence claim; `LEGACY_BMIQ_NITER` is the single source for both the
  default and the drift warning.
- C++ returns keep the old RcppArmadillo shapes (K x 1 matrices, bare `eta`
  vector) — pinned by `test-cpp-contract.R`.
