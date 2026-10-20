# betanorm 0.1.2

* `bmiq_calibration()`, `bmiq_gold_fit()` and `quantile_norm()` leave the
  random state as they found it. Before, a session with no `.Random.seed`
  had one afterwards, created by Rcpp's default wrapper around the C++
  kernels; every kernel is now exported with `rng = false`.
* `h.policy`, `on.sample.error` and `failed.sample` need an exact value:
  an abbreviation such as `h.policy = "opt"` is now an error instead of
  being completed by `match.arg()`. The defaults are unchanged.
* `src/Makevars` and `src/Makevars.win` are gone. They only added OpenMP,
  LAPACK and BLAS flags that no code used, so the package now builds with
  R's default flags.
* Results are bit-identical to 0.1.1 and the numeric snapshots are
  unchanged. Two warnings are reworded: they now point to
  `result[["failures"]]` and `result[["h.applied"]]`.

# betanorm 0.1.1

* Hardened the internal C++ kernels against malformed arguments, none of
  which the exported functions can produce: the BMIQ block gather/scatter no
  longer overflow `int` on an `NA` or very large block start, the scatter
  refuses a destination it would otherwise silently coerce, and the quantile
  kernel sorts `NaN` with a defined order and clamps every target index.
  Results are bit-identical; the numeric snapshots are unchanged.

# betanorm 0.1.0

* Relicensed from MIT to GPL (>= 3). The quantile normalizer is derived from
  preprocessCore and the beta-mixture EM from RPMM, the latter via Steve
  Horvath's BMIQcalibration, which also underlies the BMIQ calibration
  (itself adapted from Andrew Teschendorff's BMIQ). Their authors are now
  credited as copyright holders; see `inst/COPYRIGHTS`. No code change.
