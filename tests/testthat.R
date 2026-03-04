library(testthat)
library(netZooR)

# Download shared test data with error handling and caching
test_data_dir <- file.path("testthat")
dir.create(test_data_dir, showWarnings = FALSE, recursive = TRUE)

download_if_missing <- function(url, dest) {
  dest_path <- file.path(test_data_dir, dest)
  if (!file.exists(dest_path)) {
    tryCatch(
      download.file(url, dest_path, quiet = TRUE, mode = "wb"),
      error = function(e) {
        warning("Could not download ", dest, ": ", conditionMessage(e))
      }
    )
  }
}

if (!identical(Sys.getenv("NETZOOR_SKIP_DOWNLOADS"), "true")) {
  download_if_missing('https://netzoo.s3.us-east-2.amazonaws.com/netZooR/example_datasets/ppi_medium.txt', 'ppi_medium.txt')
  download_if_missing('https://netzoo.s3.us-east-2.amazonaws.com/netZooR/unittest_datasets/testDataset.RData', 'testDataset.RData')
  download_if_missing('https://netzoo.s3.us-east-2.amazonaws.com/netZooR/example_datasets/dragon/dragon_test_get_shrunken_covariance.csv', 'dragon_test_get_shrunken_covariance.csv')
  download_if_missing('https://netzoo.s3.us-east-2.amazonaws.com/netZooR/example_datasets/dragon/dragon_layer1.csv', 'dragon_layer1.csv')
  download_if_missing('https://netzoo.s3.us-east-2.amazonaws.com/netZooR/example_datasets/dragon/dragon_layer2.csv', 'dragon_layer2.csv')
  download_if_missing('https://netzoo.s3.us-east-2.amazonaws.com/netZooR/example_datasets/dragon/dragon_python_cov.csv', 'dragon_python_cov.csv')
  download_if_missing('https://netzoo.s3.us-east-2.amazonaws.com/netZooR/example_datasets/dragon/dragon_python_prec.csv', 'dragon_python_prec.csv')
  download_if_missing('https://netzoo.s3.us-east-2.amazonaws.com/netZooR/example_datasets/dragon/dragon_python_parcor.csv', 'dragon_python_parcor.csv')
  download_if_missing('https://netzoo.s3.us-east-2.amazonaws.com/netZooR/example_datasets/dragon/risk_grid_netzoopy.csv', 'risk_grid_netzoopy.csv')
  download_if_missing('https://netzoo.s3.us-east-2.amazonaws.com/netZooR/unittest_datasets/yarn/skin.rdata', 'skin.rdata')
}

test_check("netZooR")
