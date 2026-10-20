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

test_that("missing values, empty input and wrong types are refused", {
  target <- c(0.2, 0.4, 0.6, 0.8)
  obj <- matrix(c(0.5, 0.1, 0.9, 0.3), nrow = 1L)
  expect_no_error(quantile_norm(obj, target))

  obj.na <- obj
  obj.na[1L, 2L] <- NA
  expect_error(quantile_norm(obj.na, target))
  obj.na[1L, 2L] <- NaN
  expect_error(quantile_norm(obj.na, target))

  expect_error(quantile_norm(obj[0L, , drop = FALSE], target))
  expect_error(quantile_norm(obj[, 0L, drop = FALSE], target))
  expect_error(quantile_norm(as.vector(obj), target))
  expect_error(quantile_norm(obj > 0.4, target))

  expect_error(quantile_norm(obj, c(target, NA)))
  expect_error(quantile_norm(obj, c(target, Inf)))
  expect_error(quantile_norm(obj, numeric(0L)))
  expect_error(quantile_norm(obj, as.character(target)))
})

test_that("tied values share their mean rank", {
  target <- c(0.1, 0.2, 0.3, 0.4)

  # Equal lengths. A pair tied at ranks 1 and 2 takes the mean of those two
  # target values; three tied at ranks 2 to 4 take the middle one. Each row
  # is ranked on its own.
  obj <- rbind(c(1, 1, 2, 3), c(5, 5, 5, 1))
  expect_identical(
    quantile_norm(obj, target),
    rbind(
      c(0.5 * (0.1 + 0.2), 0.5 * (0.1 + 0.2), 0.3, 0.4),
      c(0.3, 0.3, 0.3, 0.1)
    )
  )

  # Unequal lengths. Mean rank 1.5 of 3 is percentile 0.25, which lands
  # exactly on the second of five target values.
  expect_identical(
    quantile_norm(matrix(c(1, 1, 2), nrow = 1L), c(0, 0.25, 0.5, 0.75, 1)),
    matrix(c(0.25, 0.25, 1), nrow = 1L)
  )
})

test_that("dimnames are kept and do not change the values", {
  target <- c(0.2, 0.4, 0.6)
  obj <- matrix(
    c(3, 1, 2, 1, 2, 3),
    nrow = 2L,
    dimnames = list(c("s1", "s2"), c("cg1", "cg2", "cg3"))
  )
  out <- quantile_norm(obj, target)
  expect_identical(dimnames(out), dimnames(obj))
  expect_identical(unname(out), quantile_norm(unname(obj), target))

  colnames(obj) <- NULL
  expect_identical(
    dimnames(quantile_norm(obj, target)),
    list(c("s1", "s2"), NULL)
  )
})
