
otter <- function(W, P, C, lambda = 0.0035, gamma = 0.335, Iter = 300, eta = 0.00001, bexp = 1){
  #W: prior gene regulatory matrix
  #P: ppi matrix
  #C: correlation matrix
  #lambda: tuning parameter 
  #gamma: regularization parameter
  #Iter: number of total iterations
  #eta: learning rate
  #bexp: exponent influencing the learning rate
  
  #ADAM parameters
  b1 <- 0.9
  b2 <- 0.999
  eps <- 0.00000001
  b1t <- b1**bexp
  b2t <- b2**bexp
  
  dW <- dim(W)
  nTF <- dW[1]
  nGenes <- dW[2]
  P <- P*((1-lambda)/sum(diag(P))) + (1-lambda)*0.0013
  C <- C*(lambda/sum(diag(C)))
  W <- P%*%W
  W <- W/sqrt(sum(diag(W%*%t(W))))
  P = P - gamma*diag(nTF)
  m <- matrix(data = 0, nrow = nTF, ncol = nGenes)
  v <- m 
  for(i in 1:Iter){
    grad <- W%*%t(W)%*%W -P%*%W - W%*%C  
       m <- b1*m + (4*(1-b1))*grad
       v <- b2*v + (16*(1-b2))*grad^2 
      b1t <- b1t*b1
      b2t <- b2t*b2
      alpha <- sqrt(1-b2t)/(1-b1t)*eta
      epst <- eps*sqrt((1-b2t))
      #update of gene ragulatory matrix
       W <- W - alpha*(m/(epst+sqrt(v)))
  }
  return(W)
}