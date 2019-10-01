context("test PANDA result")

test_that("runPanda function works", {
  
  load("./testDataset.RData")
  
  # test error message when empty inputs
  expect_error(runPanda())
  
  # file path
  T4_expression_file_path <- system.file("extdata", "expr4_matched.txt", package = "netZooR", mustWork = TRUE)
  motif_file_path <- system.file("extdata", "chip_matched.txt", package = "netZooR", mustWork = TRUE)
  ppi_file_path <- system.file("extdata", "ppi_matched.txt", package = "netZooR", mustWork = TRUE)
  
  # test message when only expression data input
  expect_message(testPanda1 <- runPanda(T4_expression_file_path))
  
  # test message when PPI is not provided
  expect_message(testPanda2 <- runPanda(T4_expression_file_path,motif_file_path,rm_missing = TRUE),"")
  actual_T4pandaList <- runPanda(T4_expression_file_path, motif_file_path,ppi_file_path)
  
  # panda network compare
  actual_T4pandaNet <- actual_T4pandaList$panda
  expect_equal(sum(abs(actual_T4pandaNet[,4]-expected_T4pandaNet[,4])),0)
  
  # panda indegree
  actual_T4panda_indegree <- actual_T4pandaList$indegree
  expect_equal(sum(abs(actual_T4panda_indegree[,1] - expected_T4panda_indegree[,1])),0)
  
  # panda outdegree
  actual_T4panda_outegree <- actual_T4pandaList$outdegree
  expect_equal(sum(abs(actual_T4panda_outegree[,1] - expected_T4panda_outdegree[,1])),0)

})



