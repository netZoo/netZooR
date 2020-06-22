context("test PANDA result")

test_that("panda function works", {

  
  # test 1: check test error message when empty inputs
  expect_error(panda.py())
  
  # input file path
  T4_expression_file_path <- system.file("extdata", "expr4.txt", package = "netZooR", mustWork = TRUE)
  motif_file_path <- system.file("extdata", "chip.txt", package = "netZooR", mustWork = TRUE)
  ppi_file_path <- system.file("extdata", "ppi.txt", package = "netZooR", mustWork = TRUE)
  
  # test 2: check message when only expression data input
  # To do 1: error occurred when only expression as input dataset
  # expect_message(panda.py(T4_expression_file_path))
  
  # test 3: check message when PPI is not provided, to do 2.
  # expect_message( panda.py(T4_expression_file_path,motif_file_path),"")
  
  # test 4: when all arguments are default
  # computing="cpu", precision="double",save_memory=FALSE, save_tmp=TRUE, keep_expression_matrix=FALSE, modeProcess="union", remove_missing=FALSE
  test1 <- panda.py(T4_expression_file_path, motif_file_path,ppi_file_path)
  test1Panda <- test1$panda
  # test 5-8: check data type in PANDA network
  expect_equal(class(test1Panda$TF), "character")
  expect_equal(class(test1Panda$Gene), "character")
  expect_equal(class(test1Panda$Motif), "numeric")
  expect_equal(class(test1Panda$Score), "numeric")
  
  # test 9: check if PANDA result is correct when processMode = union, precision = double
  # 
  expect_equal(test1Panda[1,4],-0.23212458160041557,tolerance=1e-7)
  
  # # test 10-12: check if PANDA indegree network is correct
  # test1Indegree <- test1$indegree
  # 
  # expect_equal(class(test1Indegree$`Target`), "character")
  # expect_equal(class(test1Indegree$`Target_Score`), "numeric")
  # expect_equal(test1Indegree[1,2], -272.9487716751912,tolerance=1e-7)
  # 
  # # test 13-15: check if PANDA outdegree network is correct
  # test1Outdegree <- test1$outdegree
  # 
  # expect_equal(class(test1Outdegree$`Regulator`), "character")
  # expect_equal(class(test1Outdegree$`Regulator_Score`), "numeric")
  # expect_equal(test1Outdegree[1,2], -567.4816444521832,tolerance=1e-7)
  # 
  # 
  
  # # test 16: check if PANDA result is correct when processMode = union and save_memory =F ,and the rest arguments are opposite to the default values: 
  # # precision="single",save_memory = F, save_tmp=F, keep_expression_matrix = T, modeProcess = 'union'
  # 
  # test2Panda <- panda.py(T4_expression_file_path, motif_file_path,ppi_file_path,precision = "single", save_memory = F, save_tmp = F,keep_expression_matrix = TRUE,modeProcess = "union" )$panda
  # expect_equal(test2Panda[1,4],-0.2321256, tolerance=1e-5)
  # 
  # # test 17: check if PANDA result is correct when save_memory = TRUE:
  # 
  # test3Panda <- panda.py(T4_expression_file_path, motif_file_path,ppi_file_path,save_memory = T)$WAMpanda
  # expect_equal(test3Panda[1,1],-0.23212458160041557,tolerance=1e-7)
  # 
  # # test 18：when processMode = interaction
  # test4Panda <- panda.py(T4_expression_file_path, motif_file_path,ppi_file_path, modeProcess = "intersection")$panda
  # expect_equal(test4Panda[1,4],-0.26648627625773935,tolerance=1e-7)
  # 
  # # test 19: when processMode = legacy, remove_missing=FALSE
  # test5Panda <- panda.py(T4_expression_file_path, motif_file_path,ppi_file_path, modeProcess = "legacy", remove_missing = FALSE)$panda
  # expect_equal(test5Panda[1,4],-0.06200149888611282,tolerance=1e-7)
  # 
  # # test 20: when processMode = legacy, remove_missing=TRUE
  # test6Panda <- panda.py(T4_expression_file_path, motif_file_path,ppi_file_path, modeProcess = "legacy", remove_missing = TRUE)$panda
  # expect_equal(test6Panda[1,4],-0.11907165898368743,tolerance=1e-7)
  # 

})



