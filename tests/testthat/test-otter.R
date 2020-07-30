context("test OTTER result")

test_that("otter function works", {
   load("./testDataset.RData")
   
   # Run OTTER algorithm
   otterW <- otter(W, P, C, Iter = 1)
   expect_equal(as.integer(otterW*10**10)/10**10, as.integer(gt*10**10)/10**10)  
})

