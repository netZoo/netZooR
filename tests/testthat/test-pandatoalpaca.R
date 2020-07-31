context("test panda.to.alpaca function")

test_that(" panda.to.alpaca works or not", {
  load("./testDataset.RData")
  alpaca <- panda.to.alpaca(T4Panda,T10Panda,NULL,verbose = F)
  expect_equal(alpaca, T4T10Alpaca)
}) 
