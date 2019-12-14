context("test ALPACA result")


test_that("panda.to.alpaca works", {
  load("./testDataset.RData")
  
  T4_expression_file_path <- system.file("extdata", "expr4_matched.txt", package = "netZooR", mustWork = TRUE)
  T10_expression_file_path <- system.file("extdata", "expr10_matched.txt", package = "netZooR", mustWork = TRUE)
  motif_file_path <- system.file("extdata", "chip_matched.txt", package = "netZooR", mustWork = TRUE)
  ppi_file_path <- system.file("extdata", "ppi_matched.txt", package = "netZooR", mustWork = TRUE)
  actual_T4pandaNet <- panda.py(T4_expression_file_path, motif_file_path,ppi_file_path,mode_process = "legacy")$panda
  actual_T10pandaNet <- panda.py(T10_expression_file_path, motif_file_path,ppi_file_path,mode_process = "legacy")$panda
  
  # test panda.to.alpaca()
  actual_alpaca <- panda.to.alpaca(actual_T4pandaNet, actual_T10pandaNet)
  expect_equal(actual_alpaca, expected_alpaca)

})
