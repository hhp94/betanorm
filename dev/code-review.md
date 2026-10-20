# Code review in chunks

Started 2026-10-20 on branch boost-investigate, at `dce688a` (0.1.2).
Plain ASCII throughout.

## Why

No `/code-review` has ever run on this package (searched every local
session transcript, 2026-08-07 on). Review it once, in chunks, so each
chunk's findings get checked and decided before the next one starts.

## Chunks

A review target is a path, so the smallest chunk is one file.

| # | Target | Lines | Focus |
|---|---|---|---|
| 1 | `R/quantile_norm.R`, `src/quantile_norm.cpp` | 280 | Self-contained; the trial run that calibrates noise and effort |
| 2 | `src/bmiq_norm.cpp` | 826 | EM kernel, block gather/scatter, input scan: memory safety, numerics. Check against the 2026-10-05 hardening entry |
| 3 | `R/bmiq_calibration.R` | 1564 | Thresholds, nL = 2 map, H stitching, failure options, gold fit |
| 4 | `R/betanorm-package.R`, `R/data.R`, `DESCRIPTION` | small | Optional, last |

Out of scope: `RcppExports.*` (generated) and the tests as a target (they
are the evidence; a coverage gap found in a chunk gets a test there).

## Per chunk

1. Run `/code-review high <paths>`, no `--fix`.
2. Verify every finding against the code, the CLAUDE.md invariants and
   DECISIONS.md. Classify: bug / moves a number / settled / not a bug /
   cleanup.
3. A finding that contradicts DECISIONS.md is adjudicated on its merits:
   if it stands, DECISIONS.md is updated (dated entry); if not, the
   reason is logged here.
4. Record each finding and its disposition below.
5. The user picks the fixes. Each fix is its own commit through the 0.1.2
   gate: `devtools::document()`, `devtools::test()` with snapshots
   unchanged, `dev/bitcheck-0.1.2.R` against 0.1.1. A fix that moves a
   number stops at a report: it breaks the downstream prefit and is a
   separate decision.
6. Pause before the next chunk.

## Status

- [x] 1. quantile_norm: reviewed 2026-10-20 (`high`), 6 findings. Fixed 1
  (`4c73b99`), 2 (`24aab96`) and 3 (`5fe3b4e`); kept 4, 5 and 6 as they
  are. Outcome in DECISIONS.md (2026-10-20, chunk 1).
- [x] 2. bmiq_norm.cpp: reviewed 2026-10-20 (`high`), 7 findings. The
  user took every proposal: fixed 2 and 4 (`1186b03`, one snapshot redone
  for a status string, numbers identical), logged 1 and 6, kept 3, 5 and
  7. Outcome in DECISIONS.md (2026-10-20, chunk 2).
- [x] 3. bmiq_calibration.R: reviewed 2026-10-20 (`high`), 9 findings.
  The user took every proposal: fixed 1, 3, 8 and 9, documented 2, 4 and
  5, kept 6 and 7. No number moves in a default session. Outcome in
  DECISIONS.md (2026-10-20, chunk 3, and "Legacy compatibility").
- [ ] 4. package glue (optional). Add `README.Rmd`/`README.md`: see the
  chunk 3 side observation.

## Findings

### Chunk 1: `R/quantile_norm.R`, `src/quantile_norm.cpp`

No finding contradicts DECISIONS.md. Class in brackets, then the
disposition. Every fix went through the gate: suite passes, snapshots
unchanged, bitcheck 26/26 values identical against 0.1.1.

1. **`obj` finiteness not validated** (`R/quantile_norm.R:30`) [bug].
   Docs say "Must be finite"; only `any.missing = FALSE` is checked, so
   `Inf`/`-Inf` are ranked and mapped silently (confirmed by running it).
   **Fixed** (`4c73b99`): one `checkmate::assert_numeric(obj,
   any.missing = FALSE, finite = TRUE)` owns the values, as for `target`;
   `assert_matrix()` keeps the shape. The new test failed on the old code
   (2 of 3 expectations). NEWS.md 0.1.2 bullet.
2. **`ncol(obj) == 1` maps to `min(target)`** (`src/quantile_norm.cpp:188`)
   [moves a number if changed]. Undocumented; confirmed (`0.5, 0.7` vs
   target `0.1, 0.5, 0.9` -> `0.1, 0.1`). **Fixed by documenting it**
   (`24aab96`): Rd details, a comment at the C++ branch (the downstream
   copy keeps comments, not roxygen), and a test with an unsorted target.
   No number moves.
3. **Tests don't pin ties, `ncol == 1`, validation or dimnames**
   (`test-numeric-regression.R:207`) [test gap]. Both snapshots are
   tie-free rnorm. **Fixed** (`5fe3b4e`, plus the tests in the two commits
   above): `test-quantile-norm.R`, hand-computed expectations with
   `expect_identical()`. Tests only, so no bitcheck for that commit.
4. **Two full transposes, peak ~3x `obj`** (`src/quantile_norm.cpp:179`)
   [cleanup, efficiency]. **Kept**: one-shot call, no internal caller, gain
   unmeasured.
5. **Dead `fraction == 1.0` branch, redundant `floor(target_floor + 0.5)`**
   (`src/quantile_norm.cpp:101`) [cleanup]. **Kept**: the port stays
   line-comparable with preprocessCore.
6. **Absolute `4 * DBL_EPSILON` snap ineffective at large `m`**
   (`src/quantile_norm.cpp:85`) [moves a number if changed]. **Not a
   bug**: mirrors preprocessCore; the error is about 1e-11 times the gap
   between neighbouring target values.

### Chunk 2: `src/bmiq_norm.cpp`

Checked against the 2026-10-05 hardening entry: none of its four fixes is
reopened, and nothing asks for the C++ scalar validation it declined.
Finding 6 bears on "Validation: one owner per invariant" (2026-08-10); it
agrees with that entry, so there is nothing to overturn, only a call to
make. Checks were run on the -O3 build the bitcheck installed. Class in
brackets, then the proposal.

1. **`std::lgamma` in `log_beta_norm()`** (`src/bmiq_norm.cpp:80`) [moves
   a number if changed]. The EM's -log B(a, b) comes from the platform
   libm, while digamma and trigamma come from Rmath. On this machine
   mingw's `std::lgamma` differs from R's `lgamma` at 58% of 200k points
   in (0.001, 100), by up to 2.8e-14 (2e-11 relative near the zeros at 1
   and 2). So another libm (macOS, glibc) will give slightly different EM
   numbers. Whether those numbers stay inside the 1e-12 snapshots is
   unknown: there is no CI on another platform. The cancellation part of
   the finding does not matter in practice: real fits have shapes from 1.5
   to 38 (the bitcheck cases), where the three-term form is off from
   `lbeta` by about 6e-14. Moving to `R::lgammafn` or `R::lbeta` would
   make the numbers portable but moves them, and would break the 0.1.1
   downstream prefit. Proposed: leave it, and record the platform
   dependence in DECISIONS.md. It matters again if EM threading is ever
   tried, since glibc's `lgamma` writes the global `signgam`.
2. **Retained shapes reported as "converged"** (`src/bmiq_norm.cpp:648`)
   [bug, output only, no number moves]. Confirmed on the nL = 2 snapshot
   case (2000 draws, seed 4321, `maxiter = 40`, `tol = 1e-8`). At iteration
   37 both components end at "max_iter" (the Newton fit hit `beta_maxit`).
   At iteration 38 both keep those shapes and are labelled "converged".
   The roxygen says the per-component convergence status is reported in
   `diagnostics`, so users are shown a wrong status. Fix: in the retained
   branch, keep the status of the iteration whose shapes are kept (it is
   still in `fit_status[k]`) and say in the reason that the shapes were
   kept. No number moves, but `fit_status` in the nL = 2 EM snapshot goes
   from "converged" to "max_iter", so that snapshot has to be redone,
   which is against the 0.1.2 rule. That makes it the user's call.
   `test-em-convergence.R:50` asserts `all(fit_status == "converged")`
   for its case and has to be re-checked under the fix.
3. **Non-finite concentration falls back to 1.0** (`src/bmiq_norm.cpp:151`)
   [not a bug in practice]. It is only reachable when the weighted
   variance is subnormal, and a variance of exactly 0 is already refused
   one check earlier. Proposed: leave it.
4. **Ascent guard's retained path has no test** (`src/bmiq_norm.cpp:591`)
   [mostly wrong]. The nL = 2 EM snapshot takes that path at its last
   iteration, so the path is pinned, but only through the serialized
   snapshot. Proposed: add an explicit `fit_reason` assertion only if 2 is
   fixed, as that fix's test.
5. **A division per observation in the Welford update**
   (`src/bmiq_norm.cpp:40`) [efficiency, moves a number if changed]. The
   update costs n x K divisions per EM iteration. Plain weighted sums would
   change the arithmetic. EM is about 26% of runtime (DECISIONS.md, the
   -O2 profile), and the quantile map is where the time goes. Proposed:
   leave it.
6. **The C++ `nL` in {2, 3} check duplicates R's `assert_int`**
   (`src/bmiq_norm.cpp:464`) [settled by contract]. By the one-owner rule
   it is a second owner: K is sized from `nL` and the shape check is the
   bounds guard, so `nL = 4` would be safe. But `test-cpp-contract.R:129`
   pins the message, and CLAUDE.md says to change the contract
   "deliberately or not at all". Removing the check moves no number.
   Proposed: keep it, and add it to the CLAUDE.md list of deliberate C++
   checks so a later review does not raise it again.
7. **Redundant `!std::isfinite(p)`** (`src/bmiq_norm.cpp:137`) [cleanup].
   NaN and Inf already fail `p > 0 && p < 1`. Removing it is bit-identical.
   Proposed: leave it; the churn is not worth it in a vendored file.

Side observation, not a finding: the Newton fit often stops at
`beta_maxit = 50` ("max_iter") instead of at the 1e-10 scaled-score
tolerance. In the snapshot case this happens at 14 of the 38 iterations
for one component or both. Each step is still accepted, so the fit sits
at the noise floor. The cost is per component, not per observation, so
it is small.

Dispositions (the user took every proposal above):

- 2 and 4 **fixed** (`1186b03`): the retained branch no longer writes
  `fit_status`, so it keeps the status of the fit the shapes came from. A
  new test in `test-em-convergence.R` checks that the kept shapes and the
  status match a run one iteration short. The nL = 2 EM snapshot is redone
  (`fit_status` "converged" x2 -> "max_iter" x2). Bitcheck: 25/26 cases
  identical; in `em_nL2` every field except `fit_status` is bitwise
  identical to 0.1.1. `test-em-convergence.R:50` still passes. NEWS.md
  0.1.2 updated; its "snapshots unchanged" line now names the one redone
  snapshot.
- 1 **logged**: the platform dependence is in DECISIONS.md, with the
  measurements and the thread-safety note.
- 6 **adjudicated**: the check stays, recorded in DECISIONS.md as an
  exception to the one-owner rule and added to the CLAUDE.md list of
  deliberate C++ checks.
- 3, 5 and 7 **kept as they are**, with reasons in DECISIONS.md.
- At the user's request, a new to-do: one up-front scan of the input
  matrix (`dev/to-do.md`).

### Chunk 3: `R/bmiq_calibration.R`

Two findings touch DECISIONS.md, and neither overturns it. Finding 2 sits
next to "Legacy compatibility (do not drift)". That entry pins the
`nL = 3`, `niter = 5` configuration. The H anchors differ from legacy
code, but they already did in the rewrite (`f536c1d`), before the log
began, so nothing has drifted since; the divergence was simply never
written down. Finding 8 sits next to "Validation: one owner per
invariant", which keeps the *reachable* class-count checks. The empty-H
check is unreachable, and it is a second owner of "class 2 is not
empty", so the rule asks for its removal. Checks ran on the -O3 build the
chunk 2 bitcheck installed. Class in brackets, then the proposal.

1. **The fit subset depends on the caller's RNG kind**
   (`R/bmiq_calibration.R:412`) [bug; moves numbers only in a session
   with a non-default RNG kind]. `draw_fit_indices()` calls
   `set.seed(seed)` with the caller's kind, so both the gold subset and
   the sample subset follow `RNGkind()`. Measured on one synthetic
   3000-probe sample: under `RNGkind("L'Ecuyer-CMRG")` (which parallel
   and future workers set) the calibrated values move by up to 0.021,
   and under `sample.kind = "Rounding"` by up to 0.0056. A chunked
   parallel run then disagrees with the serial run. Legacy code has the
   same dependence (`set.seed(1)`). Fix: `set.seed(seed, kind =
   "Mersenne-Twister", normal.kind = "Inversion", sample.kind =
   "Rejection")`. That is R's default, so default sessions, the bitcheck
   and existing prefits stay bit-identical. The restore that is already
   there puts the caller's kind back, because the kind is stored in
   `.Random.seed`. Proposed: fix, with a test that a non-default kind
   gives the default result and that `RNGkind()` is unchanged afterwards.
2. **The H anchors drop the legacy gap patching**
   (`R/bmiq_calibration.R:1256`) [settled in the rewrite, not recorded;
   moves a number if changed]. Both legacy references
   (`dev/bmiq-calibration-legacy.R:309-314`, MEAT) offset the anchors by
   `deltaUH`/`deltaHM`. Here H is glued straight to the calibrated U
   maximum and the upper-M minimum. The rewrite did this on purpose: the
   comment in `dev/with_nl2/bmiq-calibration.R:1494` reads "no raw-gap
   carry into output". The comment was later dropped, and DECISIONS.md
   never recorded the change. One side effect: the outer H probes tie
   with the probes on the far side of each join (3000 distinct inputs
   gave 2998 distinct outputs). This is one of several deliberate
   departures from legacy numbers: density-crossing thresholds instead of
   legacy's subset midpoints, a Newton EM with an ascent guard instead of
   BFGS / Nelder-Mead, and the continuous H join. Meanwhile the roxygen
   says "defaults target legacy BMIQ compatibility", and the README says
   the defaults "reproduce" legacy implementations. Proposed: keep the
   map. Record the divergences in DECISIONS.md, restore a one-line
   comment at the anchors (the downstream copy keeps comments), and
   reword that roxygen sentence to say the *configuration* is legacy, not
   the numbers. The README goes to chunk 4.
3. **A prefit gold is checked only for its class and `nL`**
   (`R/bmiq_calibration.R:998`) [validation gap; no number moves]. `a`,
   `b`, `thresholds` and the two modes are used per sample unchecked. An
   object saved by another version, or edited by hand, makes every
   sample fail with a misleading threshold error. Under `"continue"` the
   result is an all-NA matrix with one warning. A `NULL` `nL` gives
   "argument is of length zero". By the one-owner rule, the user-facing
   entry point owns the shape of user arguments, and nothing owns this
   one today. Proposed: fix with one check in the prefit branch. It would
   require: `nL` in {2, 3}; `a` and `b` of length `nL`, finite and
   positive; `thresholds` through `check_thresholds(require.unit.interval
   = TRUE)`; modes finite in [0, 1]. Add a test. The downstream is
   prefit-only, so it benefits too.
4. **With a prefit gold the call's EM settings don't reach the gold**
   (`R/bmiq_calibration.R:1513`) [documentation; no number moves]. The
   `niter` doc says it governs "the gold-standard fit and each sample
   fit", and `result[["settings"]]` records the call's values as though
   they applied to the gold. The gold's own settings are in the prefit's
   `settings`. Only `th1.v` is documented as ignored. Proposed: a doc
   fix. Say once, in the `goldstandard.beta` doc, that with a prefit only
   the sample fits use `nfit`, `niter`, `tol` and `beta.*`, and that the
   gold's own values are in its `settings`. Leave the result as it is.
5. **An optional H skip also leaves the H class raw**
   (`R/bmiq_calibration.R:1297`) [documentation; no number moves]. The
   `h.policy` doc and the verbose message mention only lower-M. The code,
   and legacy `doH = FALSE`, also leave every class-2 probe unchanged.
   `dev/opus_audit1.md` said the same. Proposed: fix the `h.policy` doc,
   the details paragraph and the message ("H class and lower-M left
   unchanged"). The warning and `print()` say "(U plus upper-M)", which
   is accurate, so they stay.
6. **`result[["call"]]` carries the data under `do.call()`**
   (`R/bmiq_calibration.R:911`) [efficiency; no number moves].
   `match.call()` inlines argument values when the function is called
   through `do.call()`. Measured: a 1.79 MB call object for a 24 KB
   `datM`. The field is undocumented and nothing reads it. Proposed:
   leave it. A plain call stores symbols, and removing a returned field
   is a visible change.
7. **`debug = TRUE` repeats the fit indices for every sample**
   (`R/bmiq_calibration.R:555`) [efficiency; debug only]. Every sample's
   diagnostics hold the same `sample.fit.idx`. In memory that is one
   shared vector, but `saveRDS()` writes it once per sample (about 160 MB
   at 2000 x 20000). Proposed: leave it; moving it would change the
   diagnostics shape for a debug-only cost.
8. **`unique()` on disjoint sets, and an unreachable empty-H check**
   (`R/bmiq_calibration.R:1245`) [cleanup; bit-identical]. Class 2 and
   lower-M (class 3) never overlap, so `unique(c(...))` is `c(...)` in
   the same order. `require_all_classes()` already guarantees class 2 is
   not empty. Proposed: fix by dropping both, per the one-owner rule
   (see above).
9. **The output matrix is allocated before the gold fit**
   (`R/bmiq_calibration.R:978`) [efficiency; bit-identical]. A bad raw
   gold vector, or a failed gold EM, stops the call only after the full
   n x p NA matrix is allocated and filled. At EPIC x 2000 samples that is
   about 14 GB. Proposed: fix by moving the gold block above the
   allocation. The datM scan and the probe-count check stay first. The
   gold draw restores the RNG state, so the sample subset is unchanged.

Side observation, for chunk 4: `README.Rmd`/`README.md` use
`fit.policy`, an argument that no longer exists. Line 131 passes it to a
call, and with no `...` that call errors with "unused argument". The
README also advises switching to `nL = 2` to run EM longer, which the
`nL` doc argues against, and it says the defaults "reproduce" legacy
implementations (see 2).

Dispositions (the user took every proposal above). Each fix is its own
commit through the gate: suite passes (264), snapshots as after chunk 2,
bitcheck values as after chunk 2 (only `em_nL2`'s `fit_status` differs
from 0.1.1, from `1186b03`).

- 1 **fixed** (`21a0be1`): `set.seed()` names R's default kinds. The new
  test in `test-random-state.R` failed on the old code (both value
  checks) and checks that the caller's kind is back afterwards. The one
  behaviour change outside a default session: a session that removed its
  `.Random.seed` by hand after picking another kind is left on the
  default kind (confirmed; logged in DECISIONS.md).
- 3 **fixed** (`86aa144`): `check_gold_fit()` in the prefit branch. The
  nL-mismatch message is unchanged. The new test in `test-gold-fit.R`
  failed on the old code. 0.1.1's saved gold fits pass the check and give
  the 0.1.1 numbers.
- 8 **fixed** (`bf86de6`): `c()` in place of `unique(c())`, and the empty-H
  stop removed, with a comment saying why the set cannot be empty.
- 9 **fixed** (`414d711`): the gold block now runs right after the `datM`
  checks; the messages come in the same order.
- 5 **documented** (`c23dd19`): `h.policy` doc, Details and the per-sample
  message. The bitcheck shows that message as the only new text change.
- 4 **documented** (`dc41ac6`): in the `goldstandard.beta` doc of
  `bmiq_calibration()` only. `bmiq_gold_fit()` inherits `niter` and the
  rest, so their docs stay as they are.
- 2 **documented** (`126df97`): the title paragraph and a new Details
  paragraph say the defaults are a legacy configuration, not legacy
  numbers, and name the three steps that differ. A comment is back at the
  H anchors. DECISIONS.md "Legacy compatibility" records the three steps.
- 6 and 7 **kept as they are**, with reasons in DECISIONS.md.
- NEWS.md 0.1.2: a bullet each for 1, 3 and 9 and one for the docs. The
  bit-identity line now says "under R's default random-number generator".
