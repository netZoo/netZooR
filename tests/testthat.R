library(testthat)
library(netZooR)

#download test data
getwd()
system('cd ../inst/extdata/ && { curl -O https://netzoo.s3.us-east-2.amazonaws.com/netZooR/unittest_datasets/ppi.txt ; cd -; }')
system('cd testthat/ && { curl -O https://netzoo.s3.us-east-2.amazonaws.com/netZooR/unittest_datasets/testDataset.RData ; cd -; }')

test_check("netZooR")
