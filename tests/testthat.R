library(testthat)
library(netZooR)

#download test data
system('curl -o inst/extdata/ppi.txt https://netzoo.s3.us-east-2.amazonaws.com/netZooR/unittest_datasets/ppi.txt')

test_check("netZooR")
