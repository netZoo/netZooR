# test dragon

context("test DRAGON helper functions")

test_that("[DRAGON] scale() function yields expected results", {
  x = seq(1:100)
  scale_unbiased = scale(x) # default is bias=F
  expect_equal(sd(scale_unbiased),1,tolerance = 1e-15)
  
  xbar = mean(x)
  xstd_bias = sqrt(1/length(x)*sum((x-xbar)^2))
  
  # calculate the unbiased standard deviation of the vector scaled
  # using the biased standard deviation. It will be slightly off from 1
  biased_sd_hand = sqrt(1/(length(x)-1)*sum(((x-xbar)/xstd_bias)^2))
  
  scale_biased = scale(x,bias=T)
  expect_equal(sd(scale_biased),biased_sd_hand, tolerance = 1e-15)
})

test_that("[DRAGON] VarS() function yields expected results", {}
)

test_that("[DRAGON] EsqS() function yields expected results", {}
)

