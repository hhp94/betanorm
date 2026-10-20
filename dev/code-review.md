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
- [ ] 2. bmiq_norm.cpp
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
