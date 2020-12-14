context("test OTTER result")

test_that("otter function works", {
   load("./testDataset.RData")
   
   # Run OTTER algorithm
   otterW <- otter(W, P, C, Iter = 1, lambda = 0.0035, gamma = 0.335)
   expect_equal(as.integer(otterW*10**10)/10**10, as.integer(gt*10**10)/10**10)  
})

