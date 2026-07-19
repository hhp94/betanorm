test_that("Beta-matrix has been input in the correct format in clean_beta", {
  # Check that SE is a SummarizedExperiment object
  data("GSE121961", envir = environment())
  expect_error(
    clean_beta(SE = GSE121961),
    "Please make sure SE is a SummarizedExperiment object."
  )
})
