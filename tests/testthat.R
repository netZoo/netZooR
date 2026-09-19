library(testthat)
library(netZooR)

# Defines S3_BASE_URL and netzoo_download(); test_check() loads it again for the
# individual test files, so the bucket URL lives in exactly one place.
source(file.path("testthat", "helper-netzoo-data.R"))

# Download shared test data once; already-present files are reused.
test_data_dir <- file.path("testthat")
dir.create(test_data_dir, showWarnings = FALSE, recursive = TRUE)

download_if_missing <- function(relative_path) {
  dest <- file.path(test_data_dir, basename(relative_path))
  if (!file.exists(dest)) {
    netzoo_download(relative_path, dest_dir = test_data_dir)
  }
  invisible(dest)
}

if (!identical(Sys.getenv("NETZOOR_SKIP_DOWNLOADS"), "true")) {
  download_if_missing('example_datasets/ppi_medium.txt')
  download_if_missing('unittest_datasets/testDataset.RData')
  download_if_missing('example_datasets/dragon/dragon_test_get_shrunken_covariance.csv')
  download_if_missing('example_datasets/dragon/dragon_layer1.csv')
  download_if_missing('example_datasets/dragon/dragon_layer2.csv')
  download_if_missing('example_datasets/dragon/dragon_python_cov.csv')
  download_if_missing('example_datasets/dragon/dragon_python_prec.csv')
  download_if_missing('example_datasets/dragon/dragon_python_parcor.csv')
  download_if_missing('example_datasets/dragon/risk_grid_netzoopy.csv')
}

test_check("netZooR")
