context("test panda.to.alpaca function")

test_that(" panda.to.alpaca works or not", {
  load("./testDataset.RData")
  alpaca <- panda.to.alpaca(T4Panda_subset,T10Panda_subset,NULL,verbose = F)
  expect_equal(alpaca, T4T10Alpaca)
}) 
