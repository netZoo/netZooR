context("test LIONESS result")

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

