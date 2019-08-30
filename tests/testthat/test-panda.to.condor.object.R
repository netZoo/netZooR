context("test CONDOR result")

test_that("panda.to.condor.object function works", {
  load("./testDataset.RData")
  
  T4_expression_file_path <- system.file("extdata", "expr4_matched.txt", package = "netZooR", mustWork = TRUE)
  motif_file_path <- system.file("extdata", "chip_matched.txt", package = "netZooR", mustWork = TRUE)
  ppi_file_path <- system.file("extdata", "ppi_matched.txt", package = "netZooR", mustWork = TRUE)
  actual_T4pandaNet <- runPanda(T4_expression_file_path, motif_file_path,ppi_file_path)$panda
  
  # test panda.to.condor.object
  actual_T4condor.object <- panda.to.condor.object(actual_T4pandaNet)
  expect_equal(actual_T4condor.object$modularity, expected_T4condor.object$modularity)
  
  # test when threshold is invalid
  expect_error(test1_condor.object <- panda.to.condor.object(actual_T4pandaNet,threshold=100))
  expect_error(test2_condor.object <- panda.to.condor.object(actual_T4pandaNet,threshold=-100))
})

