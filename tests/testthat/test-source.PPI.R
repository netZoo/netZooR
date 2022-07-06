context("test source PPI")

test_that("sourcePPI works", {
  load("./testDataset.RData")
  tf <- data.frame(matrix(c("Rv0022c","Rv0023","Rv0042c","Rv0043c","Rv0047c","Rv0054"), 
                          nrow=6, byrow=T),stringsAsFactors=FALSE)
  # STRINGdb Version 10
  if(R.Version()$major=="3"){
    actual_PPI_V10 <- sourcePPI(tf,"10",83332)
    ppiV10$from = as.factor(ppiV10$from)
    ppiV10$to = as.factor(ppiV10$to)
    expect_equal(actual_PPI_V10, ppiV10)
  }
  # STRINGdb Version 11
  else if(R.Version()$major=="4"){
    actual_PPI_V11 <- sourcePPI(tf,"11",83332)
    expect_equal(actual_PPI_V11, ppiV11)
  }
}) 
