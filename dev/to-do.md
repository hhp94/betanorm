# To-do

Work in flight. One section per thread; delete a section when it lands
(its outcome goes to DECISIONS.md).

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
