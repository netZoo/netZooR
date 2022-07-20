context("test LIONESS result")

test_that("lionessPy() function works", {
  
  
  # test 1: check test error message when empty inputs
  expect_error(lionessPy())
  
  # input file path
  system("curl -O  https://netzoo.s3.us-east-2.amazonaws.com/netZooR/example_datasets/expr4_800.txt")
  T4_expression_file_path <- "./expr4_800.txt"
  motif_file_path <- system.file("extdata", "chip.txt", package = "netZooR", mustWork = TRUE)
  ppi_file_path <- "./ppi.txt"
  
  # test 2: check message when only expression data input
  #expect_message(lionessPy(T4_expression_file_path, end_sample=1, save_fmt='no', save_single_network=TRUE), regexp="motif network", fixed=TRUE)
  
  # test 3: check message when PPI is not provided
  #expect_message(lionessPy(T4_expression_file_path,motif_file_path, end_sample=1, save_fmt='no', save_single_network=TRUE), regexp="No PPI", fixed=TRUE)
  
  # test 4: when all arguments are default, except end_sample = 1 to expedite computing.
  # computing="cpu", precision="double", save_tmp=TRUE, modeProcess="union", remove_missing=FALSE, start_sample=1, end_sample=1, save_single_network=FALSE
  test1Lioness <- lionessPy(T4_expression_file_path, motif_file_path,ppi_file_path, end_sample=1, save_fmt='no', save_single_network=TRUE)
  expect_equal(test1Lioness[[1,3]],-0.02283242, tolerance=1e-5)
  
  # test 5: check if LIONESS result is correct when arguments set as following:
  # i.e computing = "cpu", save_memory =T , precision="single", save_tmp=F, keep_expression_matrix = T, modeProcess = 'intersection',remove_missing=FALSE, start_sample=1, end_sample=1, save_single_network=FALSE
  test2Lioness <- lionessPy(T4_expression_file_path, motif_file_path,ppi_file_path,precision = "single", save_tmp = F, modeProcess = "intersection", 
                            remove_missing=FALSE, start_sample=1, end_sample=1, save_single_network=TRUE, save_fmt='no')
  expect_equal(test2Lioness[[1,3]],13.00665, tolerance=1e-5)
  
  # test 6: when processMode = legacy, remove_missing=FALSE
  test3Lioness <- lionessPy(T4_expression_file_path, motif_file_path,ppi_file_path, 
                            modeProcess = "legacy", remove_missing = FALSE,start_sample=1, end_sample=1, save_single_network=TRUE, save_fmt='no')
  expect_equal(test3Lioness[[1,3]],11.4918,tolerance=1e-5)
  
  # test 7: when processMode = legacy, remove_missing=TRUE
  test4Lioness <- lionessPy(T4_expression_file_path, motif_file_path,ppi_file_path, 
                            modeProcess = "legacy", remove_missing = TRUE,start_sample=1, end_sample=1, save_single_network=TRUE, save_fmt='no')
  expect_equal(test4Lioness[[1,3]],-0.5589099,tolerance=1e-5)
  
})

test_that("lioness() function works for network.inference.method = 'panda'", {
  data(pandaToyData)
  test5Lioness <- lioness(expr = pandaToyData$expression[,1:4], 
                        motif = pandaToyData$motif, ppi = pandaToyData$ppi, network.inference.method = 'panda')
  expect_equal(test5Lioness[[1]][1],-0.6704147,tolerance=1e-5)
})

test_that("lioness() function works for network.inference.method = 'pearson'", {

  # create toy data
  x1 = c(1,1,2,3)
  x2 = c(1,2,3,4)
  X = matrix(c(x1,x2),ncol=2,byrow=F)
  
  # extract correlation networks needed for LIONESS
  fullNet = cor(X)
  subNet1 = cor(X[-1,])
  subNet2 = cor(X[-2,])
  subNet3 = cor(X[-3,])
  subNet4 = cor(X[-4,])
  
  # construct LIONESS sample-specific correlation networks manually
  N = 4
  ssNet1 = N*(fullNet-subNet1) + subNet1
  ssNet2 = N*(fullNet-subNet2) + subNet2
  ssNet3 = N*(fullNet-subNet3) + subNet3
  ssNet4 = N*(fullNet-subNet4) + subNet4
  
  # run LIONESS on the toy data
  lionessNets = lioness(expr=t(X),network.inference.method="pearson")
  
  # assert equal
  # there is no file I/O or rounding here
  # so differences should be machine zero in R
  expect_equal(lionessNets[[1]],ssNet1,tolerance=1e-15) 
  expect_equal(lionessNets[[2]],ssNet2,tolerance=1e-15) 
  expect_equal(lionessNets[[3]],ssNet3,tolerance=1e-15) 
  expect_equal(lionessNets[[4]],ssNet4,tolerance=1e-15) 
  
})

test_that("lioness() function throws appropriate error for undefined network.inference.method", {

  # create toy data
  x1 = c(1,1,2,3)
  x2 = c(1,2,3,4)
  X = matrix(c(x1,x2),ncol=2,byrow=F)
  
  # run LIONESS on the toy data w/undefined method
  expect_error(lioness(expr=t(X),network.inference.method="divination"))
  
})

