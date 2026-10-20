# To-do

Work in flight. One section per thread; delete a section when it lands
(its outcome goes to DECISIONS.md).

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

## Gold quantile table + Newton polish (branch: bmiq-fixed; next)

**Goal:** replace the gold-side `stats::qbeta()` of the quantile map with
a precomputed table plus one Newton step, first in
`bmiq_calibration_fixed()` (whose gold is fixed, so its tables can be
built once at load or even stored), then decide about the general
function.

**Why:** the gold shapes are constant for a whole run (the gold is
consumed distributionally), so `qbeta(p, a.gold, b.gold)` is one fixed
function per component and tail; Rmath's `qbeta` is a crude initial
guess plus Newton iterations on `pbeta`, so a table is only a better
initial guess. The sample-side `pbeta` varies per sample and is cheap.

**Prototype (2026-10-20, `dev/qtable-proto.R`, `dev/qtable-extra.R`,
untracked):** per gold component, two halves (Beta(a, b) and Beta(b, a)),
lower-tail quantile for p in [1e-30, 1/2] on a grid uniform in
u = log p (n = 4096), storing v = log x and the slope
dv/du = p / (x f(x)) = exp(u - v - dbeta(x, a, b, log = TRUE)); cubic
Hermite in u, then one Newton step on log F(x) = u with
`pbeta(log.p = TRUE)` and the closed-form log density
(a-1) log x + (b-1) log1p(-x) - lbeta(a, b). The mirror
Q_upper(p; a, b) = 1 - Q_lower(p; b, a), and log(1 - p) = log(-expm1(u))
for p > 1/2 with the swapped half, cover both tails and the nL = 2
log-p path (which Boost could not). Exact ends at p = 0 and 1; fall back
to `qbeta` below p.min (never hit).

Measured against Rmath on every gold-side `qbeta` input of the six EPIC
samples (210k evaluations, 5 call sites, nL = 2 and 3) and an
adversarial shape sweep (0.5..40 x 0.5..40, both tails, p from 1e-28 to
1 - 1e-12):

- max |table + 1 Newton - qbeta| = 3.3e-16 (6.7e-16 in the sweep);
  Hermite alone ~1e-10 at n = 4096, 3e-12 at 8192, so the Newton step is
  required and one suffices (grid >= 2048).
- 43.6% of outputs bit-identical, 99.5% within 1 ulp: not bit-identical,
  so a deliberate numerics change under CLAUDE.md (DECISIONS.md entry,
  NEWS, snapshots at 1e-12 pass unchanged, the bitcheck will differ at
  the last ulp), as the single-exp E-step was.
- per eval: Rmath `qbeta` 1.97 us, `pbeta` 0.29, `dbeta` 0.20, table
  Hermite 0.11, table + 1 Newton 0.67 (0.55 with the closed-form
  density). One sample's map (14778 evals): 34 ms -> 12.5 ms (2.7x).
  Table build 36 ms per run (two gold components). Projected ~1.8x end
  to end at 21k probes, ~2.5-3x at 450k.

**Adoption bar:** max |diff| vs Rmath <= 1e-14 on the real workload and
the sweep; the suite and snapshots pass unchanged; -O2 end-to-end bench
before and after. Serial, no dependency, no threads.

## Boost + RcppThread quantile-map bench (branch: boost-investigate)

The table thread above reaches most of this goal serially and without
forking the numerics beyond the last ulp; keep this section until the
table lands, then close it as superseded or run it for the remainder.

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

**Mechanics:** bench lives in `dev/boost_bench/` (not shipped; opt each
file worth keeping into `.gitignore`);
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
