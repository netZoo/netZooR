library(testthat)
library(netZooR)

#download test data
download.file('https://netzoo.s3.us-east-2.amazonaws.com/netZooR/example_datasets/ppi_medium.txt','testthat/ppi_medium.txt')
download.file('https://netzoo.s3.us-east-2.amazonaws.com/netZooR/unittest_datasets/testDataset.RData','testthat/testDataset.RData')
download.file('https://netzoo.s3.us-east-2.amazonaws.com/netZooR/example_datasets/dragon/dragon_test_get_shrunken_covariance.csv','testthat/dragon_test_get_shrunken_covariance.csv')
download.file('https://netzoo.s3.us-east-2.amazonaws.com/netZooR/example_datasets/dragon/dragon_layer1.csv','testthat/dragon_layer1.csv')
download.file('https://netzoo.s3.us-east-2.amazonaws.com/netZooR/example_datasets/dragon/dragon_layer2.csv','testthat/dragon_layer2.csv')
download.file('https://netzoo.s3.us-east-2.amazonaws.com/netZooR/example_datasets/dragon/dragon_python_cov.csv','testthat/dragon_python_cov.csv')
download.file('https://netzoo.s3.us-east-2.amazonaws.com/netZooR/example_datasets/dragon/dragon_python_prec.csv','testthat/dragon_python_prec.csv')
download.file('https://netzoo.s3.us-east-2.amazonaws.com/netZooR/example_datasets/dragon/dragon_python_parcor.csv','testthat/dragon_python_parcor.csv')
download.file('https://netzoo.s3.us-east-2.amazonaws.com/netZooR/example_datasets/dragon/risk_grid_netzoopy.csv','testthat/risk_grid_netzoopy.csv')

test_check("netZooR")
