context("test ALPACA result")

test_that("ALPACA works", {
  
  load("./testDataset.RData")
  
  simp.alp <- alpaca(simp.mat,NULL,verbose=F)
  
  expect_equal(as.vector(simp.alp[[1]]),simp.memb)
  
  
})


