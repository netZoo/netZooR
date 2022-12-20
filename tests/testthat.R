library(testthat)
library(netZooR)

#download test data
download.file('https://netzoo.s3.us-east-2.amazonaws.com/netZooR/example_datasets/ppi_medium.txt','testthat/ppi_medium.txt')
download.file('https://netzoo.s3.us-east-2.amazonaws.com/netZooR/unittest_datasets/testDataset.RData','testthat/testDataset.RData')
download.file('https://netzoo.s3.us-east-2.amazonaws.com/netZooR/example_datasets/dragon_test_get_shrunken_covariance.csv','testthat/dragon_test_get_shrunken_covariance.csv')

test_check("netZooR")
