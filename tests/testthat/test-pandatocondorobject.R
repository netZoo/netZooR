context("test panda.to.condor.object function")

test_that(" panda.to.condor.object works or not", {
  load("./testDataset.RData")
  condor.obj <- pandaToCondorObject(T4Panda_subset)
  expect_equal(condor.obj$edges, T4Condor$edges)
}) 
