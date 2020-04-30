context("test OTTER result")

test_that("otter function works", {
   W=as.matrix(read.csv('./w.csv', header = FALSE))
   C=as.matrix(read.csv('./c.csv', header = FALSE))
   P=as.matrix(read.csv('./p.csv', header = FALSE))
   gt=as.matrix(read.csv('./test_otter.csv', header = FALSE))

   # Run OTTER algorithm
   W <- otter(W, P, C, Iter = 1)
   assert_that((all((as.integer(W*10**2)/10**2) == (as.integer(gt*10**2)/10**2))))  
})
