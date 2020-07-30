context("test SAMBAR result")

test_that("SAMBAR function works", {
  
  load("./testDataset.RData")
  data("exon.size")
  data("mut.ucec")
  data("genes")
  expect_equal(sambar(mutdata=mut.ucec, esize=exon.size, signatureset=system.file("extdata", "h.all.v6.1.symbols.gmt", package="netZooR", mustWork=TRUE), 
                     cangenes=genes, kmin=2, kmax=4),subtypes)
  
  
})

