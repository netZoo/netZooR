context("test source PPI")

test_that("sourcePPI works", {
  load("./testDataset.RData")
  actual_PPI <- sourcePPI(tf,83332)
  expect_equal(actual_PPI,PPI.form.StringDB)
})
