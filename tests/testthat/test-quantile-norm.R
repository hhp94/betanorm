# quantile_norm(): what it refuses, and the cases the snapshots in
# test-numeric-regression.R do not reach.

test_that("obj must be finite", {
  target <- c(0.2, 0.4, 0.6, 0.8)
  obj <- matrix(c(0.5, 0.1, 0.9, 0.3), nrow = 1L)

  # The finite matrix runs, so each error below comes from the one value
  # changed.
  expect_no_error(quantile_norm(obj, target))

  obj.inf <- obj
  obj.inf[1L, 1L] <- Inf
  expect_error(quantile_norm(obj.inf, target))

  obj.inf[1L, 1L] <- -Inf
  expect_error(quantile_norm(obj.inf, target))
})

test_that("one column maps every value to min(target)", {
  # A row with one value has nothing to rank against. target is unsorted, so
  # this pins the minimum, not the first element.
  obj <- matrix(c(0.5, 0.7, 0.1), ncol = 1L)
  expect_identical(
    quantile_norm(obj, c(0.9, 0.1, 0.5)),
    matrix(0.1, nrow = 3L, ncol = 1L)
  )

  # With a length-1 target the lengths are equal; the answer is the same.
  expect_identical(
    quantile_norm(obj, 0.4),
    matrix(0.4, nrow = 3L, ncol = 1L)
  )
})
