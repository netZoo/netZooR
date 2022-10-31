context("test PANDA result")

test_that("panda function works", {

   
  # test 1: check test error message when empty inputs
   expect_error(pandaPy())
   
   # input file path
   system("curl -O  https://netzoo.s3.us-east-2.amazonaws.com/netZooR/example_datasets/expr4_200_L.txt")
   T4_expression_file_path <- "./expr4_200_L.txt"
   motif_file_path <- system.file("extdata", "chip_medium.txt", package = "netZooR", mustWork = TRUE)
   ppi_file_path <- "./ppi_medium.txt"
  
  # test 2: check message when only expression data input
  # To do 1: error occurred when only expression as input dataset
  expect_message(pandaPy(T4_expression_file_path))

  # test 3: check message when PPI is not provided, to do 2.
  expect_message( pandaPy(T4_expression_file_path,motif_file_path))

  # test 4: when all arguments are default
  # computing="cpu", precision="double",save_memory=FALSE, save_tmp=TRUE, keep_expression_matrix=FALSE, modeProcess="union", remove_missing=FALSE
  test1Panda<- pandaPy(T4_expression_file_path, motif_file_path, ppi_file_path)$panda
  expect_equal(test1Panda[1,4],-0.08132568,tolerance=1e-7)
   
   # test 5: check if PANDA result is correct when arguments settiing like below:
   # i.e computing = "cpu", save_memory =T , precision="single", save_memory = T, save_tmp=F, keep_expression_matrix = T, modeProcess = 'intersection'
   test2Panda <- pandaPy(T4_expression_file_path, motif_file_path,ppi_file_path,precision = "single", save_memory = T, save_tmp = F,keep_expression_matrix = TRUE, modeProcess = "intersection" )$WAMpanda
   expect_equal(test2Panda[1,1],2.229422, tolerance=1e-5)
  
   # test 6: when processMode = legacy, remove_missing=FALSE
   test3Panda <- pandaPy(T4_expression_file_path, motif_file_path,ppi_file_path, modeProcess = "legacy", remove_missing = FALSE)$panda
   expect_equal(test3Panda[1,4],5.86307,tolerance=1e-7)
  
   # test 7: when processMode = legacy, remove_missing=TRUE
   test4List  <- pandaPy(T4_expression_file_path, motif_file_path,ppi_file_path, modeProcess = "legacy", remove_missing = TRUE)
   test4Panda <-  test4List$panda
   
   # test 8-11: check data type in PANDA network
   expect_equal(class(test4Panda$TF), "character")
   expect_equal(class(test4Panda$Gene), "character")
   expect_equal(class(test4Panda$Motif), "numeric")
   expect_equal(class(test4Panda$Score), "numeric")
   
   # test 12: check PANDA network result
   expect_equal(test4Panda[1,4],-0.4123061,tolerance=1e-7)
   
   
   # test 13-15: check if PANDA indegree network is correct
   test4Indegree <- test4List$indegree
   
   expect_equal(class(test4Indegree$`Target`), "character")
   expect_equal(class(test4Indegree$`Target_Score`), "numeric")
   expect_equal(test4Indegree[1,2], 0.1914127,tolerance=1e-7)
   
   # test 15-18: check if PANDA outdegree network is correct
   test4Outdegree <- test4List$outdegree
   
   expect_equal(class(test4Outdegree$`Regulator`), "character")
   expect_equal(class(test4Outdegree$`Regulator_Score`), "numeric")
   expect_equal(test4Outdegree[1,2], -12.23053,tolerance=1e-6)
  

})
