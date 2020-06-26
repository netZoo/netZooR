#' Run OTTER in R
#' 
#' Description:
#'               OTTER infers gene regulatory networks using TF DNA binding
#'               motif (W), TF PPI (P), and gene coexpression (C) through 
#'               minimzing the following objective:
#'                                  min f(W) 
#'               with f(W) = (1-lambda)*||WW' - P||^2 + lambda*||W'W - C||^2 + (gamma/2)*||W||^2
#'
#' Inputs:
#' @param W     : TF-gene regulatory network based on TF motifs as a
#'                       matrix of size (t,g), g=number of genes, t=number of TFs
#' @param P     : TF-TF protein interaction network as a matrix of size (t,t)
#' @param C     : gene coexpression as a matrix of size (g,g) 
#' @param lambda: tuning parameter in [0,1] (higher gives more weight to C)
#' @param gamma : regularization parameter
#' @param Iter  : number of iterations of the algorithm
#' @param eta   : learning rate
#' @param bexp  : exponent influencing learning rate (higher means smaller)
#'
#' Outputs:
#' @return W    : Predicted TF-gene complete regulatory network as an adjacency matrix of size (t,g).
#'
#' @examples
#'
#' W=matrix(rexp(100, rate=.1), ncol=10)
#' C=matrix(rexp(100, rate=.1), ncol=10)
#' P=matrix(rexp(100, rate=.1), ncol=10)
#'
#' # Run OTTER algorithm
#' W <- otter(W, P, C)
#'  
#' @export

otter <- function(W, P, C, lambda = 0.0035, gamma = 0.335, Iter = 32, eta = 0.00001, bexp = 1){
  #ADAM parameters
  b1 <- 0.9
  b2 <- 0.999
  eps <- 0.00000001
  b1t <- b1**bexp
  b2t <- b2**bexp
  
  dW <- dim(W)
  nTF <- dW[1]
  nGenes <- dW[2]
  P <- P+2.2
  P <- P/sum(diag(P))
  W <- P%*%W
  P <- P*(1-lambda)
  C <- C*(lambda/sum(diag(C)))
  W <- W/sum(diag(W%*%t(W)))
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
