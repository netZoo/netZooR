context("test LIONESS result")

test_that("lionessPy() function works", {
  
  
  # test 1: check test error message when empty inputs
  expect_error(lionessPy())
  
  # input file path
  system("curl -O  https://netzoo.s3.us-east-2.amazonaws.com/netZooR/example_datasets/expr4.txt")
  T4_expression_file_path <- "./expr4.txt"
  motif_file_path <- system.file("extdata", "chip.txt", package = "netZooR", mustWork = TRUE)
  ppi_file_path <- "./ppi.txt"
  
  # test 2: check message when only expression data input
  # To do 1: error occurred when only expression as input dataset
  # expect_message(lionessPy(T4_expression_file_path))
  
  # test 3: check message when PPI is not provided, to do 2.
  # expect_message(lionessPy(T4_expression_file_path,motif_file_path),"")
  
  # test 4: when all arguments are default, except end_sample = 1 to expedite computing.
  # computing="cpu", precision="double", save_tmp=TRUE, modeProcess="union", remove_missing=FALSE, start_sample=1, end_sample=1, save_single_network=FALSE
  # check lioness result
  test1Lioness <- lionessPy(T4_expression_file_path, motif_file_path,ppi_file_path, end_sample=1)
  expect_equal(test1Lioness[[1,3]],-0.01674112, tolerance=1e-5)
  
  # test 5: check if LIONESS result is correct when arguments set as following:
  # i.e computing = "cpu", save_memory =T , precision="single", save_tmp=F, keep_expression_matrix = T, modeProcess = 'intersection',remove_missing=FALSE, start_sample=1, end_sample=1, save_single_network=FALSE
  test2Lioness <- lionessPy(T4_expression_file_path, motif_file_path,ppi_file_path,precision = "single", save_tmp = F, modeProcess = "intersection", 
                             remove_missing=FALSE, start_sample=1, end_sample=1, save_single_network=FALSE)
  expect_equal(test2Lioness[[1,3]],-0.2456164, tolerance=1e-5)
  
  # test 6: when processMode = legacy, remove_missing=FALSE
  test3Lioness <- lionessPy(T4_expression_file_path, motif_file_path,ppi_file_path, modeProcess = "legacy", remove_missing = FALSE,start_sample=1, end_sample=1, save_single_network=FALSE)
  expect_equal(test3Lioness[[1,3]],-0.165805,tolerance=1e-5)
  
  # test 18: when processMode = legacy, remove_missing=TRUE
  test4Lioness <- lionessPy(T4_expression_file_path, motif_file_path,ppi_file_path, modeProcess = "legacy", remove_missing = TRUE,start_sample=1, end_sample=1, save_single_network=FALSE)
  expect_equal(test4Lioness[[1,3]],-0.1407832,tolerance=1e-5)
  

})



