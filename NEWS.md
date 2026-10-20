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
* `quantile_norm()` refuses an `obj` holding `Inf` or `-Inf`, as its
  documentation already required. Before, they were ranked like any other
  value and mapped to the ends of `target`.
* When the mixture EM keeps a component's previous shapes because the new
  fit did not improve the likelihood, the component status in the
  `debug = TRUE` diagnostics is now that of the fit the shapes came from
  (for example `"max_iter"`). Before, it always said `"converged"`.
* `bmiq_calibration()` and `bmiq_gold_fit()` draw their fit subsets under
  R's default random-number generator, so the results no longer depend on
  the session's `RNGkind()`. Before, a session using another kind (for
  example `"L'Ecuyer-CMRG"`, which parallel workers often set) got
  different calibrated values. The caller's generator is restored as
  before.
* `bmiq_calibration()` checks a prefitted `bmiq_gold_fit` (its `nL`,
  shapes, thresholds and modes) before calibrating. Before, a malformed
  one made every sample fail with a misleading error, or gave an all-`NA`
  result under `on.sample.error = "continue"`.
* With a raw gold vector, `bmiq_calibration()` fits the gold standard
  before allocating the output matrix, so a bad gold fails sooner.
* The documentation now says that the defaults follow the legacy BMIQ
  configuration but do not reproduce its numbers, and why; that a skipped
  H step leaves the H class unchanged as well as lower-M (the per-sample
  message says so too); and which settings a prefitted gold ignores.
* Results are bit-identical to 0.1.1 under R's default random-number
  generator. One numeric snapshot is redone, for the component status
  above only; its numbers are unchanged. Two warnings are reworded: they
  now point to `result[["failures"]]` and `result[["h.applied"]]`.

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
