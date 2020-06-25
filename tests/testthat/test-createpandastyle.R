context("test create.panda.style() function")

test_that("test create.panda.style() function", {
  
  # lazy test, not rigorous enough.
  # can not connect to cytoscape
  expect_error(create.panda.style(),"Failed to connect to localhost port 1234: Connection refused")
  
})

