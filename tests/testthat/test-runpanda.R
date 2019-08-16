context("test-runpanda")

test_that("runPanda function works", {
  
  expect_error(runPanda())
  T4_expression_file_path <- system.file("extdata", "expr4_matched.txt", package = "netZooR", mustWork = TRUE)
  T10_expression_file_path <- system.file("extdata", "expr10_matched.txt", package = "netZooR", mustWork = TRUE)
  motif_file_path <- system.file("extdata", "chip_matched.txt", package = "netZooR", mustWork = TRUE)
  ppi_file_path <- system.file("extdata", "ppi_matched.txt", package = "netZooR", mustWork = TRUE)
  expect_message(testPanda1 <- runPanda(T4_expression_file_path))
  expect_message(testPanda2 <- runPanda(T4_expression_file_path,motif_file_path,rm_missing = TRUE),"")
  actual_T4pandaNet <- runPanda(T4_expression_file_path, motif_file_path,ppi_file_path)$panda
  expected_T4pandaNet <- read.csv("./testdata/expr4_PandaNet.txt",sep = " ",header = F)
  expect_equal(sum(abs(actual_T4pandaNet[,4]-expected_T4pandaNet[,4])),0)
})



