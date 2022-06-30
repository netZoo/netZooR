context("test pandaDiffEdges() result")

test_that("pandaDiffEdges() works or not", {
  
  load("./testDataset.RData")
  
  diffPanda <- pandaDiffEdges(T4Panda_subset, T10Panda_subset)
  expect_equal(diffPanda,T4T10diffPanda)
  
})