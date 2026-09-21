# Test datasets are served from the netZoo S3 bucket. A failed download raises
# instead of silently leaving an HTTP error page on disk in place of the data.
S3_BASE_URL <- "https://netzoo-data.s3.us-east-2.amazonaws.com/netZooR"

netzoo_download <- function(relative_path, dest_dir = ".") {
  url <- paste(S3_BASE_URL, relative_path, sep = "/")
  dest <- file.path(dest_dir, basename(relative_path))
  fail <- function(reason) {
    stop("Could not download test data from ", url, ": ", reason, call. = FALSE)
  }
  status <- withCallingHandlers(
    tryCatch(
      utils::download.file(url, dest, quiet = TRUE, mode = "wb"),
      error = function(e) fail(conditionMessage(e))
    ),
    warning = function(w) fail(conditionMessage(w))
  )
  if (!identical(as.integer(status), 0L) || !file.exists(dest) ||
      file.size(dest) == 0) {
    fail("empty or incomplete response")
  }
  invisible(dest)
}
