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
- [ ] 3. bmiq_calibration.R
- [ ] 4. package glue (optional)

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
