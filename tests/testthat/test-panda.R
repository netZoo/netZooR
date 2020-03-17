context("test PANDA result")

test_that("panda function works", {
  
  load("./testDataset.RData")
  
  # test error message when empty inputs
  expect_error(panda.py())
  
  # file path
  T4_expression_file_path <- system.file("extdata", "expr4_matched.txt", package = "netZooR", mustWork = TRUE)
  motif_file_path <- system.file("extdata", "chip_matched.txt", package = "netZooR", mustWork = TRUE)
  ppi_file_path <- system.file("extdata", "ppi_matched.txt", package = "netZooR", mustWork = TRUE)
  
  # test message when only expression data input
  #expect_message(testPanda1 <- panda.py(T4_expression_file_path))
  
  # test message when PPI is not provided
  # expect_message(testPanda2 <- panda.py(T4_expression_file_path,motif_file_path,rm_missing = TRUE),"")
  actual_T4pandaList <- panda.py(T4_expression_file_path, motif_file_path,ppi_file_path,mode_process = "legacy")
  
  # access PANDA network
  actual_T4pandaNet <- actual_T4pandaList$panda
  # check data type in PANDA network
  expect_equal(class(actual_T4pandaNet$TF), "character")
  expect_equal(class(actual_T4pandaNet$Gene), "character")
  expect_equal(class(actual_T4pandaNet$Motif), "numeric")
  expect_equal(class(actual_T4pandaNet$Score), "numeric")
  
  # panda network compare
  expect_equal(sum(abs(actual_T4pandaNet[,4]-expected_T4pandaNet[,4])),0)
  
  # panda indegree & data type check
  actual_T4panda_indegree <- actual_T4pandaList$indegree
  expect_equal(class(actual_T4panda_indegree$`Target`), "character")
  expect_equal(class(actual_T4panda_indegree$`Target_Score`), "numeric")
  
  expect_equal(sum(abs(actual_T4panda_indegree[,2] - expected_T4panda_indegree[,1])),0)
  
  # panda outdegree & data type check
  actual_T4panda_outdegree <- actual_T4pandaList$outdegree
  expect_equal(class(actual_T4panda_outdegree$`Regulator`), "character")
  expect_equal(class(actual_T4panda_outdegree$`Regulator_Score`), "numeric")
  
  expect_equal(sum(abs(actual_T4panda_outdegree[,2] - expected_T4panda_outdegree[,1])),0)

})



