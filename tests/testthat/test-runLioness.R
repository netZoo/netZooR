context("test LIONESS result")

test_that("runLioness function works", {
  
  # load expected output dataset
  load("./testDataset.RData")
  
  # test error message when empty inputs
  expect_error(runLioness())
  
  # file path
  T4_expression_file_path <- system.file("extdata", "expr4_matched.txt", package = "netZooR", mustWork = TRUE)
  motif_file_path <- system.file("extdata", "chip_matched.txt", package = "netZooR", mustWork = TRUE)
  ppi_file_path <- system.file("extdata", "ppi_matched.txt", package = "netZooR", mustWork = TRUE)
  
  # test message when PPI is not provided
  expect_message(testPanda2 <- runLioness(T4_expression_file_path,motif_file_path,rm_missing = TRUE),"")
  
  # run LIONESS with T4 dataset
  actual_T4lioness<- runLioness(T4_expression_file_path, motif_file_path,ppi_file_path)
  
  # Lioness sample 1 network compare
  expect_equal(sum(abs(actual_T4lioness[,3]-expected_T4lioness_samp1_list$T4_Samp1)),0)
  
})
