# betanorm 0.1.4

* The mixture EM is about twice as fast: per E-step it takes one `log`
  per 256 values instead of one per value, and it skips the `exp` of each
  value's largest term, which is exactly 1. On the six EPIC samples,
  `bmiq_calibration()` takes 0.26 s instead of 0.29 s, the internal
  `bmiq_calibration_fixed()` 0.09 s instead of 0.12 s, and
  `bmiq_gold_fit()` of the Horvath gold 8 ms instead of 14. The calibrated
  values and the fitted shapes are unchanged on those samples. Only the
  log-likelihood moves, in its last digits: its convergence criterion by
  at most 4e-11 relative. So a fit could stop one iteration earlier or
  later only if that criterion sits that close to `tol`.

# betanorm 0.1.3

* New internal `bmiq_calibration_fixed()`, not exported:
  `bmiq_calibration()` at its defaults (the legacy three-state
  configuration) against a built-in gold-standard fit chosen by name, for
  pipelines that copy the code. It fits every value instead of a random
  subset of 20000, so it draws no random numbers: `gold = "horvath"` is
  the fit of all 21,368 values of `horvath_goldstandard$goldstandard2`,
  stored so nothing is refitted for the gold, and each sample is fitted on
  all of its probes. It takes `datM` and `gold` and returns the calibrated
  matrix. A
  sample that cannot be calibrated becomes a row of `NA`; one warning says
  that some samples did, another that some ran without H. Above 200
  samples it reports progress on one line, in steps of 20%.
* The gold side of its quantile map reads a table of the gold quantiles
  (quintic Hermite on 8192 nodes in log p, within 1e-14 of
  `stats::qbeta()`), built ahead of time and stored with the package, in
  place of `stats::qbeta()`, which adds at most 1e-14. Because it fits
  every value, its output differs from that of
  `bmiq_calibration(datM, horvath_goldstandard$goldstandard2)`, which fits
  a random subset, by up to 9e-4 on the six EPIC samples. On those samples
  it takes 0.12 s against 0.29 s for `bmiq_calibration()`, which still
  calls `stats::qbeta()` and is unchanged.
* Its source file and the stored table (`FIXED_QTABLE` in
  `R/sysdata.rda`) stand alone.
* `bmiq_calibration()` and `bmiq_gold_fit()` keep the session's
  random-number generator kind when the session has no `.Random.seed`.
  Before, a session that chose another kind (for example
  `"L'Ecuyer-CMRG"`) and then removed `.Random.seed` was left on R's
  default kind. Results are unchanged.

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
* `horvath_goldstandard` is a plain data frame. It was a readr tibble
  carrying a column spec and a `problems` pointer that no longer worked
  once loaded. The columns and values are unchanged.
* The datasets' help pages give their sources: `horvath_goldstandard`
  comes from the additional files of Horvath (2013), under CC BY 2.0, and
  `GPL21145_sample` from NCBI GEO. `?GPL21145_sample` says that the 1062
  probes the EPIC array does not measure are filled with 0, and shows how
  to drop them.
* The README no longer passes the removed `fit.policy` argument. It now
  shows `bmiq_gold_fit()` and drops the zero-filled probes from its
  example.
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
