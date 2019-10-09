context("test source PPI")

test_that("source.PPI works", {
  load("./testDataset.RData")
  actual_PPI <- source.PPI(tf,83332)
  expect_equal(actual_PPI,PPI.form.StringDB)
})
