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
  so this would be bit-identical. In `bmiq_calibration_fixed()` the fit
  is every value (2026-10-21), so the two are each row's own smallest
  value above 0 and largest below 1.
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
