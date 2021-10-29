context("test pandaToCondorObject function")

test_that(" pandaToCondorObject works or not", {
  load("./testDataset.RData")
  condor.obj <- pandaToCondorObject(T4Panda_subset)
  expect_equal(condor.obj$edges, T4Condor$edges)
}) 
