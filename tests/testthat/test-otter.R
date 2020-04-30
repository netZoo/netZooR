context("test OTTER result")

test_that("otter function works", {
   W=as.matrix(read.csv('../../data/w.csv', header = FALSE))
   C=as.matrix(read.csv('../../data/c.csv', header = FALSE))
   P=as.matrix(read.csv('../../data/p.csv', header = FALSE))
   gt=as.matrix(read.csv('../../data/test_otter.csv', header = FALSE))

   # Run OTTER algorithm
   W <- otter(W, P, C)
   assert_that((all((as.integer(W*10**2)/10**2) == (as.integer(gt*10**2)/10**2))))  
})
