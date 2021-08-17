library(testthat)
library(netZooR)

#download test data
system('curl -o ../inst/extdata/ppi.txt https://netzoo.s3.us-east-2.amazonaws.com/netZooR/unittest_datasets/ppi.txt')
system('curl -o testthat/testDataset.RData https://netzoo.s3.us-east-2.amazonaws.com/netZooR/unittest_datasets/testDataset.RData')

test_check("netZooR")
