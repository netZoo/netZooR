context("test create.diff.panda.style function")

test_that("test create.diff.panda.style function", {
  
  # lazy test, not rigorous enough.
  # can not connect to cytoscape
  expect_error(create.diff.panda.style(),"Failed to connect to localhost port 1234: Connection refused")
  
})
