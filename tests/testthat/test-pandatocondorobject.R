context("test panda.to.condor.object function")

test_that(" panda.to.condor.object works or not", {
  load("./testDataset.RData")
  condor.obj <- panda.to.condor.object(T4Panda)
  expect_equal(condor.obj$modularity, T4Condor$modularity)
}) 
