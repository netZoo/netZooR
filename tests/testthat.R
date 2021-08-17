library(testthat)
library(netZooR)

#download test data
system('cd ../inst/extdata/ && { curl -O https://netzoo.s3.us-east-2.amazonaws.com/netZooR/unittest_datasets/ppi.txt ; cd -; }')
system('curl -O testDataset.RData https://netzoo.s3.us-east-2.amazonaws.com/netZooR/unittest_datasets/testDataset.RData')

test_check("netZooR")
