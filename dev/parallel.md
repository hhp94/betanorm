# Per-sample parallelism — status (not implemented)

This is a **future design note only**. There is no sample-level parallel
BMIQ path in production.

## Current state

| Piece | Status |
|-------|--------|
| Per-sample `RcppThread` workers | **Not implemented** |
| Gold fit + sample loop | Still serial R (`BMIQcalibration`) |
| EM / Beta Newton in C++ | Serial; uses R math APIs where needed from the R side for quantile maps / modes |
| Boost.Math wrappers | Ready for workers: `inst/include/bmiqpp/boost_math.hpp` + `boost_*_cpp` equality tests |
| `LinkingTo: RcppThread` | Present for a future worker build; not used by the algorithm today |

Production calibration still calls R `density()` / `pbeta()` / `qbeta()` for
modes and quantile maps. Those are not thread-safe for worker threads.

## Intended architecture (when implemented)

1. R validates inputs, fits gold once, builds subsample indices serially.
2. R calls one batched C++ entry (probes × samples layout).
3. Each `RcppThread` worker owns one sample column; pure C++ only (no
   `Rcpp::stop`, no R RNG, no R special functions).
4. Workers use Boost.Math digamma/trigamma/pbeta/qbeta and return status
   structs; main thread assembles the R result.

No nested parallelism: sample-parallel first, EM remains serial per sample.

Do not treat this file as a checklist of completed work.
