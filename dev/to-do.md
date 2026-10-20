# To-do

Work in flight. One section per thread; delete a section when it lands
(its outcome goes to DECISIONS.md).

## Code review in chunks (branch: boost-investigate)

Started 2026-10-20 at `dce688a` (0.1.2). First `/code-review` of the
package, one chunk at a time: quantile_norm, `src/bmiq_norm.cpp`,
`R/bmiq_calibration.R`, then optional package glue. Plan, per-chunk
process, status and findings: `dev/code-review.md`. Fixes go through
the 0.1.2 bit-identity gate; a finding that contradicts DECISIONS.md is
adjudicated and the log updated.

## One up-front scan of the input matrix (not started)

Raised 2026-10-20 during chunk 2. **Goal:** every value check and
per-sample fact the code derives from the input matrix comes from one pass
over it, done before any sample is processed, with one owner per invariant
(DECISIONS.md, 2026-08-10). Bit-identical, like every other refactor.

What touches the input values today:

- `bmiq_calibration(datM)`: `checkmate::assert_matrix(mode = "double")`
  (type and shape only, no value pass), then
  `scan_finite_unit_interval_cpp(datM)` (one pass: finite and [0, 1]).
- Per sample, in `fit_mixture()`: the endpoint clip finds the smallest
  value above 0 and the largest below 1 in the fit subset, which builds
  two temporary vectors (`positive`, `below.one`). The fit subset
  (`sample.fit.idx`) is the same for every sample, so these two numbers
  per sample could come out of the up-front scan. Min and max are exact,
  so this would be bit-identical.
- `debug = TRUE`: `range(beta2.v)` per sample, another pass over the whole
  sample (`input_range` in the diagnostics).
- The gold vector, when raw: `as.numeric()` copy, `assert_numeric()`, then
  its own `scan_finite_unit_interval_cpp()`. A prefit skips all of this.
- `quantile_norm(obj, target)`: `assert_matrix()` for the shape, one
  `assert_numeric(any.missing = FALSE, finite = TRUE)` pass over `obj`
  (since chunk 1), and one over `target`.
- The downstream does its value scan per sample instead (DECISIONS.md,
  2026-10-05, "divergences"), so a change here also bears on what it
  copies.

Open questions: return the per-sample clip bounds (and the debug ranges)
from an extended scan kernel, or keep the scan as a pure check and add a
second kernel; whether the gold and `quantile_norm()` checks belong in the
same design; and the speed gain, which is unmeasured and probably small
next to the quantile map. Measure before building.

## Boost + RcppThread quantile-map bench (branch: boost-investigate)

**Goal:** decide whether an opt-in `engine = "boost"` mode (threaded
Boost `ibeta`/`ibeta_inv` quantile maps) is worth forking the numerics.

**Why:** corrected -O2 profile at 21k (see DECISIONS.md, -O2 entry):
`map_beta_q` is **67.6%** of runtime (3.65 us per qbeta+pbeta eval,
~12.9k evals/sample), EM 26.3%, other 6.1%. Amdahl ceiling for threading
the map alone: ~3.1x (~2.5x at 8 threads; 2000 samples ~2.4 min -> ~1 min).

**Plan — three gates, each can kill it cheaply:**

1. **Fidelity on the real workload.** Harvest fitted
   `(a.sample, b.sample, a.gold, b.gold, tail)` triples from ~50 real 21k
   fits; evaluate the *composed* map on those samples' actual probe values
   plus adversarial points (near thresholds, near 0/1). Max + quantile
   |diff| vs Rmath per call site. Run Boost under BOTH policies: default
   (promotes to long double internally — 80-bit x87 on Windows gcc, plain
   double on macOS arm64, i.e. NOT cross-platform stable) and
   `promote_double<false>` (same arithmetic everywhere; the likely
   contract if adopted, since macOS is the motivation).
2. **Speed + thread scaling at the real grain.** Per-eval cost:
   R-vectorized Rmath / C++-loop Rmath / Boost serial (both policies) /
   Boost + `RcppThread::parallelFor` at 2/4/8/14 threads. Sizes: n = 13k
   (one sample's map — thread-pool overhead vs ~47 ms serial work) and
   n = 450k (array scale; earlier OpenMP experiment saw 5-7x there).
3. **End-to-end prototype.** `assignInNamespace` swap of `map_beta_q` to
   the Boost/threaded path (no package surgery): wall s/sample on 40x21k,
   projected 2000-sample time, max |delta calibrated| vs the Rmath run.

**Pre-registered adoption bar:** opt-in engine only if workload max
|diff| <= 1e-8 AND end-to-end >= 2x at 8 threads; otherwise record in
DECISIONS.md as measured-and-rejected.

**Mechanics:** bench lives in `dev/boost_bench/` (tracked, not shipped);
compile via `Rcpp::sourceCpp` (standard -O2 Makeconf flags — NEVER -O0
for header-only Boost; have the harness print its effective flags).
Threading stays inside pure C++ math (Boost `ignore_error` policies, no
R API from workers — the Rmath thread hazard in DECISIONS.md does not
apply). Prior art: `dev/inst/include/bmiqpp/boost_math.hpp` (tail-variant
mapping done), `dev/boost_math_exports.cpp`, `dev/test-boost-math.R`
(1e-8/1e-10 agreement on small grids). Installed: BH 1.90.0.1,
RcppThread 2.4.0, rcpptimer 1.2.1.

**Deliberately out of scope (noted for later):**

- Threading the EM E-step (pure std::exp/log, no Rmath — but the loglik
  reduction order would change bits; EM is now only ~26%).
- The nL = 2 log-space path: Boost has no `log.p` incomplete beta, so an
  engine swap initially covers the legacy nL = 3 path only (which is what
  the clock pipeline runs).
