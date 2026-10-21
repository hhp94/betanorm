# Derived from preprocessCore (normalize.quantiles.use.target, get_ranks),
# copyright Ben Bolstad, LGPL (>= 2). Distributed here under GPL (>= 3).
# See inst/COPYRIGHTS.

# betanorm's quantile_norm() under another name, which betanorm does not
# export; betanorm's tests pin the two identical(). Each row of obj (samples
# in rows) is ranked and mapped onto target. obj is a finite numeric matrix,
# target a finite numeric vector of any length, in any order.
quantile_norm_fixed <- function(obj, target) {
  checkmate::assert_matrix(
    obj,
    mode = "numeric",
    min.rows = 1L,
    min.cols = 1L
  )
  checkmate::assert_numeric(
    obj,
    any.missing = FALSE,
    finite = TRUE,
    .var.name = "obj"
  )
  checkmate::assert_numeric(
    target,
    any.missing = FALSE,
    min.len = 1L,
    finite = TRUE,
    .var.name = "target"
  )
  qnorm_target_rows_cpp(obj, target)
}
