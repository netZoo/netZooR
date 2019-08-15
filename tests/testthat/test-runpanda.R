context("test-runpanda")

test_that("runPanda function works", {
  
  expect_error(runPanda())
  treated_expression_file_path <- system.file("extdata", "expr4_matched.txt", package = "netZooR", mustWork = TRUE)
  motif_file_path <- system.file("extdata", "chip_matched.txt", package = "netZooR", mustWork = TRUE)
  ppi_file_path <- system.file("extdata", "ppi_matched.txt", package = "netZooR", mustWork = TRUE)
  expect_message(testPanda1 <- runPanda(treated_expression_file_path))
  expect_message(testPanda2 <- runPanda(treated_expression_file_path,motif_file_path,rm_missing = TRUE),"")
})
