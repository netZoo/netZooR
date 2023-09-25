#' Run COBRA in R
#' 
#' Description:
#'               COBRA decomposes a (partial) gene co-expression matrix as a 
#'               linear combination of covariate-specific components. 
#'               It can be applied for batch correction, differential co-expression 
#'               analysis controlling for variables, and to understand the impact of
#'               variables of interest to the observed co-expression. 
#'
#' Inputs:
#' @param X               : design matrix of size (n, q), n = number of samples, q = number of covariates
#' @param expressionData  : gene expression as a matrix of size (g, n), g = number of genes
#' @param standardize     : boolean flag to standardize the gene expression as a pre-processing step
#'
#' Outputs:
#' @return psi : impact of each covariate on the eigenvalues as a matrix of size (q, n)
#' @return Q   : eigenvectors corresponding to non-zero eigenvalues as a matrix of size (g, n)
#' @return D   : non-zero eigenvalues as a list of length n
#' @return G   : (standardized) gene expression as a matrix of size (g, n)
#'
#' @examples
#'
#' g <- 100 # number of genes
#' n <- 10 # number of samples
#' q <- 2 # number of covariates
#' X <- X <- cbind(rep(1, n), rbinom(n, 1, 0.5))
#' expressionData=matrix(rnorm(g*n, 1, 1), ncol = n, nrow = g)
#'
#' # Run COBRA algorithm
#' cobra_output <- cobra(X, expressionData)
#'
#' @export  

cobra <- function(X, expressionData, standardize=T){
  numSamples <- ncol(expressionData)
  N <- min(ncol(expressionData),nrow(expressionData))
  
  if (standardize){
    G_star <- expressionData-rowMeans(expressionData)
    G <- (G_star/sqrt(rowSums(G_star^2)))
    G <- as.matrix(G)
  } else {
    G <- expressionData
    G <- (G/sqrt(rowSums(G^2)))
    G <- as.matrix(G)
  }
  
  eigenG <- rARPACK::eigs_sym(tcrossprod(G),N)
  
  Q <- eigenG$vectors
  D <- diag(eigenG$values)
  
  hatmat <- ginv(crossprod(X))%*%t(X)
  Qinv <- ginv(Q) 
  QinvG <- Qinv%*%(G)
  
  est <- t(sapply(seq_len(nrow(hatmat)), function(hatmatRow){
    diag(QinvG%*%(numSamples*diag(hatmat[hatmatRow,]))%*%t(QinvG))
  }))
  
  list(psi=est, Q=Q, D=eigenG$values, G=G)
}