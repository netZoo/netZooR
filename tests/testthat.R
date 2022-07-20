library(testthat)
library(netZooR)

#download test data
download.file('https://netzoo.s3.us-east-2.amazonaws.com/netZooR/unittest_datasets/ppi_medium.txt','testthat/ppi_medium.txt')
download.file('https://netzoo.s3.us-east-2.amazonaws.com/netZooR/unittest_datasets/testDataset.RData','testthat/testDataset.RData')

test_check("netZooR")
