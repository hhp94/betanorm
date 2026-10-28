# Blast radius: `skip` and in-place for the fixed functions

Scoped 2026-10-28 on branch `mod-inplace`; nothing is built yet. Two
requests for `bmiq_calibration_fixed()` and `quantile_norm_fixed()`:

1. a `skip` vector, one logical per row of the matrix (per sample);
2. writing the result into the input matrix in place.

**Decided (user, 2026-10-28): `quantile_norm()` gets both too.** The
`identical()` clone pin (`test-quantile-norm.R`, DECISIONS.md 2026-10-21)
stays, so the exported function takes the same arguments as its copy. The
rest is open; see "Open questions" at the end.

## Summary

No number should move. Rows are calibrated independently, and both kernels
finish reading a row before they write it, so a processed row is bit for
bit what it is today. The snapshots, and the 1e-14 pin of the fixed path to
`bmiq_calibration()`, stand. What moves is contracts: a public function
that writes into its argument, the "input comes back as found" rule, and
how skipped rows are validated.

**Untouched:** `bmiq_calibration()`, `bmiq_gold_fit()`, the EM kernel,
`FIXED_GOLD`, `FIXED_QTABLE` and its `.rds`, `data-raw/`, the numeric
snapshots, and the gather/scatter kernel code.

## 1. `skip`

Checked with `checkmate::assert_logical(skip, any.missing = FALSE, len =
nrow(...))` after the matrix assertion. Not a choice, so the
`assert_choice()` rule does not apply.

### `bmiq_calibration_fixed()`: R only

- In the block loop a skipped sample does not go through
  `process_sample()`. Its block column keeps the gathered values and the
  scatter writes them back, so in copy mode a skipped row passes through
  unchanged at no cost. Returning NA instead needs an explicit fill.
- **Bug to avoid:** the "calibrated without H" warning tests
  `any(!h.applied & !failed)`. A skipped sample has both `FALSE`, so it
  would fire for every skip. It needs `& !skip`.
- The progress line counts `nrow(datM)` (text pinned in
  `test-bmiq-fixed.R`: "Calibrating 201 samples"). Decide whether it counts
  all rows or only the samples processed, and whether the 200-sample
  threshold does.
- Pin: `fixed(datM, skip)[!skip, ]` is `identical()` to
  `fixed(datM[!skip, ])`. Exact, because samples do not interact (the gold
  is fixed, nothing is drawn, and a block is only a copy).

### `quantile_norm_fixed()` and `quantile_norm()`: R and C++

- `qnorm_target_rows_cpp()` loops over every row, so it takes a mask.
  `Rcpp::Nullable<Rcpp::LogicalVector> skip = R_NilValue` keeps the
  two-argument calls in `test-cpp-contract.R` working.
- The result is allocated `no_init` ("every cell is written below"). With
  a skip that no longer holds: skipped rows must be filled explicitly, or
  they come back holding whatever was in memory. The kernel comment
  changes with it.
- Both R files get the same edit, and the clone pin enforces it.

### Validation of skipped rows (both functions)

Today's checks cover the whole matrix: `scan_finite_unit_interval_cpp()`
for the calibration, `checkmate::assert_numeric(finite = TRUE)` for
`obj`. If rows are skipped *because* they hold NA or out-of-range values,
the call is refused before anything runs.

Exempting skipped rows needs a row-masked scan: a change to
`scan_finite_unit_interval_cpp()` in `bmiq_norm.cpp` (the general function
and the gold fit call it too, and `test-cpp-contract.R` pins it) or a new
kernel, and a C++ finite scan for `obj`, since checkmate cannot mask rows
without copying them. That overlaps the "One up-front scan" thread in
`dev/to-do.md`. A per-sample check as each row comes up is ruled out once
writes are in place (see "Not atomic" below).

## 2. In place

This reverses a pinned rule. All of these say the input comes back as
found:

- `test-input-preserved.R`, its header and both tests (the calibration
  functions; `quantile_norm()`, `quantile_norm_fixed()` and the kernel);
- the scatter kernel's comment in `bmiq_norm.cpp`: "safe only when it is a
  freshly allocated private matrix";
- DECISIONS.md, "`bmiq_calibration_fixed()` is internal, quiet...": "the
  output is a fresh matrix, so the scatter writes into it without touching
  `datM`".

So in place is opt-in (a flag, `checkmate::assert_flag()`), the default
path keeps every test it has now, and the change gets a DECISIONS.md entry.

### `bmiq_calibration_fixed()`: no kernel change

Scatter into `datM` instead of `calibrated`, and do not allocate
`calibrated`. Each block is copied out (gathered) before anything is
written, and blocks do not overlap, so the arithmetic is the same. Saves
one `datM`-sized allocation. `datM` is already asserted double, so the
scatter's type check passes. The scatter comment (a copied file) is
rewritten to name the in-place caller.

With `skip`, the scatter writes a skipped row's own gathered bits back to
it: no change, but a write. A masked scatter would avoid it; not needed.

### `quantile_norm_fixed()` and `quantile_norm()`: kernel change

In place is safe for the arithmetic: `normalize_row()` copies the whole row
into `items` before it writes anything. The kernel can take the
destination as the input itself, either by a flag on
`qnorm_target_rows_cpp()` or by a second export. Either way the
exported signatures pinned by `test-cpp-contract.R` change, and the new
export is `rng = false`.

**The coercion trap** (the scatter kernel's, DECISIONS.md 2026-10-05, item
2): with `const Rcpp::NumericMatrix &obj`, an integer matrix is converted
into a temporary, the writes go there, and they are lost with no error.
`quantile_norm()` accepts integer `obj` today (`mode = "numeric"`), so in
place must require double: in R (the assertion), and by `TYPEOF()` on a
`SEXP` parameter in C++, as the scatter does. Copy mode keeps accepting
integer.

Saves one `obj`-sized allocation: 343 MB at 100 x 450k (NEWS 0.1.5).

### Because `quantile_norm()` is exported

- Its signature gains `skip` and the in-place flag: roxygen `@param` for
  each, `@return` for the in-place case, `man/quantile_norm.Rd`
  regenerated with `roxygen2::roxygenise()`. `NAMESPACE` is unchanged.
- An exported R function that writes into its argument is unusual
  (data.table's `set*()` functions are the known precedent). The help page
  has to say what it does to other references (below). Default off.
- README.Rmd / README.md list `quantile_norm()` and show a call; a line on
  the new arguments may belong there. NEWS gets a user-visible entry.

### Risks, both functions

- **Other references change too.** `y <- x; f(x, in.place = TRUE)` changes
  `y`. This cannot be detected: `MAYBE_SHARED()` is always true by the
  time the kernel sees the matrix, since the R wrapper binds it as a
  formal (the scatter comment says so already).
- **Temporaries are lost.** `f(as.matrix(df), ...)` or `f(x[keep, ], ...)`
  writes into a copy nobody holds. Returning the matrix (invisibly) makes
  `x <- f(x, in.place = TRUE)` right in every case.
- **Failed samples (calibration).** Copy mode makes them NA rows. In
  place, NA overwrites the caller's original row, and the warnings do not
  say which rows failed. The alternative is to leave the row as it was
  (the general function's `failed.sample = "original"`). Either way the
  caller may need the failed rows back, which a returned matrix does not
  carry without an attribute or a different return type.
- **Not atomic.** Writes go out per block of 8 samples (calibration) or per
  row (quantile). A Ctrl-C mid-run leaves part of the matrix done, with no
  record of which part. So every check must run before the first write,
  which is true today and must stay true.
- **ALTREP matrices** (memory-mapped, or wrappers): writing through
  `REAL()` can land in read-only or shared memory. An `ALTREP()` refusal
  in C++ would be a memory-safety guard, like the scatter's type check,
  not a second validation owner.

## Files touched

- Copied files: `R/bmiq_calibration_fixed.R`, `R/quantile_norm_fixed.R`,
  `src/quantile_norm.cpp`, `src/bmiq_norm.cpp` (the scatter comment only,
  unless the masked scan goes there). Every comment in them stays true in
  the downstream; the downstream code rules apply (`rng = false`, `[[`,
  qualified calls, no `match.arg()`).
- Not copied: `R/quantile_norm.R`, `man/quantile_norm.Rd`, `R/RcppExports.R`
  and `src/RcppExports.cpp` (regenerated), README.Rmd / README.md.
- Tests: `test-bmiq-fixed.R` (skip pins, the H warning, progress text),
  `test-quantile-norm.R` (skip, in place, the pin unchanged),
  `test-cpp-contract.R` (kernel signature, mask, in-place, the integer
  refusal), `test-input-preserved.R` (default paths kept; in-place paths
  added: the mutated matrix is `identical()` to the copy-mode result,
  skipped rows untouched).
- Records: `dev/DECISIONS.md` (new entry; amend the "internal, quiet"
  entry's fresh-output line), `NEWS.md`, `DESCRIPTION` to 0.1.6,
  probably a line in CLAUDE.md's invariants on the opt-in mutation.
- Downstream: the new arguments have defaults, so its calls keep working;
  it moves its pinned commit and regenerates its `RcppExports`.

## Open questions

1. In copy mode, does a skipped row come back as the original values or
   as NA?
2. Must skipped rows still be finite (and in [0, 1] for the calibration)?
   If not, a masked scan is needed (section 1, validation).
3. In place as a flag on the same functions, or separate functions? With
   the clone pin, a separate function would also exist in both quantile
   files.
4. In place, does a failed calibration sample become NA or stay as it was?
   Should the failed and skipped rows be returned, and how (attribute,
   list)?
5. Argument names (`skip`, `in.place`) and the `skip` default (`NULL`, or
   `rep.int(FALSE, nrow(...))`).
6. Refuse ALTREP matrices in place?
