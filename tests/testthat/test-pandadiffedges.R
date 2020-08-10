context("test panda.diff.edges() result")

test_that("panda.diff.edges() works or not", {
  
  load("./testDataset.RData")
  
  diffPanda <- panda.diff.edges(T4Panda_subset, T10Panda_subset)
  expect_equal(diffPanda,T4T10diffPanda)
  
})