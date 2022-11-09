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

test_that("[DRAGON] VarS() function yields expected results", {
  myX = matrix(c(1,2,9,3,1,7,5,12,8),byrow=T,ncol=3)
  # gold standard calculated from DRAGON in netZooPy
  pyVarS = matrix(c(4,37,1,37,342.25,9.25,1,9.25,0.25),byrow=T,ncol=3)
  expect_equal(VarS(myX),pyVarS, tolerance=1e-15)
}
)

test_that("[DRAGON] EsqS() function yields expected results", {
  myX = matrix(c(1,2,9,3,1,7,5,12,8),byrow=T,ncol=3)
  # gold standard calculated from DRAGON in netZooPy
  pyEsqS = matrix(c(16,100,1,100,1369,0.25,1,0.25,1),byrow=T,ncol=3)
  expect_equal(EsqS(myX),pyEsqS, tolerance=1e-15)
}
)

test_that("[DRAGON] risk() function calculates correct risk",{
  Gamma1 = 1.1
  Gamma2 = 2
  T11 = 3
  T12 = 4
  T21 = 5
  T22 = 6
  T3 = 7
  T4 = 8
  # manual calc
  R_hand = (1-Gamma1^2)*T11 + (1-Gamma2^2)*T12 + 
    (1-Gamma1^2)^2*T21 +(1-Gamma2^2)^2*T22 +
    (1-Gamma1^2)*(1-Gamma2^2)*T3 + 
    Gamma1*Gamma2*T4
  R = risk(Gamma1,Gamma2,T11,T12,T21,T22,T3,T4)
  expect_equal(R,R_hand, tolerance = 1e-15)
}
)
