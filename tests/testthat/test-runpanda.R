context("test-runpanda")

test_that("runPanda function works", {
  
  expect_error(runPanda())
  T4_expression_file_path <- system.file("extdata", "expr4_matched.txt", package = "netZooR", mustWork = TRUE)
  T10_expression_file_path <- system.file("extdata", "expr10_matched.txt", package = "netZooR", mustWork = TRUE)
  motif_file_path <- system.file("extdata", "chip_matched.txt", package = "netZooR", mustWork = TRUE)
  ppi_file_path <- system.file("extdata", "ppi_matched.txt", package = "netZooR", mustWork = TRUE)
  expect_message(testPanda1 <- runPanda(T4_expression_file_path))
  expect_message(testPanda2 <- runPanda(T4_expression_file_path,motif_file_path,rm_missing = TRUE),"")
  actual_T4pandaList <- runPanda(T4_expression_file_path, motif_file_path,ppi_file_path)
  # panda net
  actual_T4pandaNet <- actual_T4pandaList$panda
  expected_T4pandaNet <- read.csv("./testdata/T4_PandaNet.txt",sep = " ",header = F)
  expect_equal(sum(abs(actual_T4pandaNet[,4]-expected_T4pandaNet[,4])),0)
  # panda indegree
  actual_T4panda_indegree <- actual_T4pandaList$indegree
  expected_T4panda_indegree <- read.csv("./testdata/T4_PandaIndegree.txt",sep = "\t",header = T,row.names = 1)
  expect_equal(sum(abs(actual_T4panda_indegree[,1] - expected_T4panda_indegree[,1])),0)
  

})



