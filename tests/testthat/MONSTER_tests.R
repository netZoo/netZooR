library(MONSTER)
context("Data Checking")

test_that("Data is in an acceptable format", {
    df <- data.frame(a=rnorm(100),b=rnorm(100),c=rnorm(100),d=rnorm(100))
    mt <- matrix(rnorm(1000),ncol=10)
    eset <- ExpressionSet(assayData=mt)
    expect_equal(checkDataType(df), as.matrix(df))
    expect_equal(checkDataType(mt), mt)
    expect_equivalent(checkDataType(eset), mt)
    
})
