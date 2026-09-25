# betanorm: relicense to GPL and credit the reference code

Written 2026-09-25 for an agent working in the hhp94/betanorm repository. Self-contained.
Plain ASCII throughout. Single hyphens only in prose a person reads.

## Why

betanorm (MIT, version 0.0.1, copyright "betanorm authors") was pair-coded with two other
packages' source open as the reference:

- `preprocessCore` (Ben Bolstad, LGPL >= 2) for `quantile_norm.R` / `src/quantile_norm.cpp`.
  The roxygen and README already say "port of normalize.quantiles.use.target()", and the
  C++ cites `get_ranks()` for tie handling.
- `RPMM` (E. Andres Houseman and Devin C. Koestler, GPL >= 2) for the beta-mixture EM in
  `src/bmiq_norm.cpp` and the fitting half of `R/bmiq_calibration.R` (the `blc()` /
  `betaEst()` logic).

Code written by reading and translating those sources is a derivative work. LGPL allows
conversion to GPL; GPL >= 2 does not allow MIT. So betanorm's license is wrong as declared
and the fix is to declare what it actually is: GPL. CRAN policy also requires that where
code is derived from another project, that project's copyright holders are listed.

The user has stated the derivation and has no objection to GPL. Do not spend effort on a
clean-room rewrite to keep MIT; that is not the ask.

## Target

- `License: GPL (>= 3)`. GPL-3 satisfies LGPL >= 2 (section 3 of LGPL 2.1 permits
  relicensing under GPL 2 or later), GPL >= 2, and GPL-3 (ChAMP, wateRmelon) in case any
  BMIQ.R-derived line is found (see Open question). Do not pick GPL-2 alone.
- Copyright holders credited in `Authors@R`, with the file each one is credited for.
- A per-file provenance record the reader can verify.

## Changes

1. `DESCRIPTION`
   - `License: GPL (>= 3)`
   - `Authors@R`: keep Hung Pham as aut, cre, cph. Add:
     ```
     person("Ben", "Bolstad", role = "cph",
            comment = "quantile normalization, derived from preprocessCore (LGPL >= 2)"),
     person("E. Andres", "Houseman", role = "cph",
            comment = "beta-mixture EM, derived from RPMM (GPL >= 2)"),
     person("Devin C.", "Koestler", role = "cph",
            comment = "beta-mixture EM, derived from RPMM (GPL >= 2)")
     ```
     `cph` only. They did not author this package, so not `aut` and not `ctb`.
   - Remove the `LICENSE` file reference: `GPL (>= 3)` needs no `+ file LICENSE`.
2. Delete `LICENSE` and `LICENSE.md` (the MIT texts). Do not add a copy of the GPL text;
   R ships it and CRAN asks packages not to.
3. Add `inst/COPYRIGHTS`, plain text, one block per derived file:
   ```
   R/quantile_norm.R, src/quantile_norm.cpp
     Derived from preprocessCore (normalize.quantiles.use.target, get_ranks),
     copyright Ben Bolstad, LGPL (>= 2). Ported to a no-missing-data, row-wise
     form in C++.

   src/bmiq_norm.cpp, R/bmiq_calibration.R (mixture fitting)
     Derived from RPMM (blc, betaEst), copyright E. Andres Houseman and
     Devin C. Koestler, GPL (>= 2). Rewritten in C++ with Armijo line search
     and per-component convergence reporting.

   All other files: copyright Hung Pham, GPL (>= 3).
   ```
   Adjust the function names to whatever the code actually followed; do not invent a
   derivation that is not there. If a file was written from the paper alone, say so.
4. File headers. At the top of each derived file, one comment naming the origin, its
   author and its license, matching `inst/COPYRIGHTS`. Keep it to three lines.
5. `README.md`: the "focused port of preprocessCore" wording is accurate and stays.
   Add a Licensing section of three sentences: the package is GPL-3, why (derived from
   two GPL-family packages), and a pointer to `inst/COPYRIGHTS`. Add the RPMM and
   preprocessCore citations beside the Teschendorff one.
6. `NEWS.md`: one entry. "Relicensed from MIT to GPL (>= 3). The quantile normalizer is
   derived from preprocessCore and the beta-mixture EM from RPMM; both are now credited
   as copyright holders. No code change."
7. Bump `Version` to 0.1.0. A license change is not a patch.

## Do not

- Do not touch any function body. This is a metadata change; a behaviour diff would
  confuse the relicense with a fix.
- Do not add the BMIQ paper's authors as `cph`. A paper is not code. Teschendorff stays
  a citation unless the Open question below says otherwise.
- Do not run R CMD check unprompted. `devtools::test()` and `devtools::document()`.

## Open question, ask the user before step 3

Was any part of `R/bmiq_calibration.R` (the threshold, class-mode and transformation
steps, as opposed to the mixture fitting) written with BMIQ.R (Teschendorff) or its copy
in ChAMP / wateRmelon open as the reference? If yes, add Andrew Teschendorff as `cph`
with a comment naming the file, and note that ChAMP / wateRmelon are GPL-3, which the
chosen `GPL (>= 3)` already covers. If no, say so in `inst/COPYRIGHTS`.

## Consequence elsewhere, recorded here so nobody is surprised

- methylCIPHERv2 vendors these files (`R/normalize_bmiq.R`, `R/normalize_quantile.R`,
  `src/bmiq_norm.cpp`) under BSD-3. That is now a known license defect in that package.
  The user is not contributing to it; the fact is recorded in dev/handoff-meeting.md.
- The successor package (dev/successor-plan.md) either goes GPL (>= 3) itself, which is
  the simple answer and the one the plan now takes, or depends on betanorm from CRAN
  rather than vendoring it, which needs betanorm on CRAN first.
