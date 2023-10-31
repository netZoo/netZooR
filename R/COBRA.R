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
#' @param method     : if pearson, the decomposition of the co-expression matrix is compouted. If pcorsh, COBRA decomposes the regularized partial co-expression
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

cobra <- function(X, expressionData, method = "pearson"){
  
  if(!(method %in% c("pearson", "pcorsh", "dragon"))){
    stop("Only Pearson and pcor methods are supported. Please make sure to provide a valid method argument.")
  }
  
  numSamples <- ncol(expressionData)
  if(method != "dragon"){
    N <- min(ncol(expressionData),nrow(expressionData)) 
  }
  C <- 0
  if(method == "pearson"){
    G_star <- expressionData-rowMeans(expressionData)
    expressionData <- (G_star/sqrt(rowSums(G_star^2)))
    expressionData <- as.matrix(expressionData)
    C <- tcrossprod(expressionData)
  }
  if(method == "dragon"){
    if(length(expressionData) != 2){
      stop("Dragon needs two layers, please provide a list of two matrices in expressionData or consider using the pcorsh argument")
    }
    if(dim(expressionData[[1]])[2] != dim(expressionData[[2]])[2]){
      stop("Dragon layers in expressionData list must have the same number of samples")
    }
    C <- dragon(t(expressionData[[1]]), layer2 = t(expressionData[[2]]), pval=FALSE)$cov
    N <- dim(expressionData[[1]])[2]
    numSamples <- N
    expressionData <- rbind(expressionData[[1]], expressionData[[2]])
    print(dim(expressionData))
  }
  if(method == "pcorsh"){
   C <- matrix(as.numeric(pcor.shrink(t(expressionData))), dim(expressionData)[1], dim(expressionData)[1])
  }
  
  eigenG <- rARPACK::eigs_sym(C,N)
  
  Q <- eigenG$vectors
  D <- diag(eigenG$values)
  
  hatmat <- ginv(crossprod(X))%*%t(X)
  Qinv <- ginv(Q) 
  QinvG <- Qinv%*%(expressionData)
  
  est <- t(sapply(seq_len(nrow(hatmat)), function(hatmatRow){
    diag(QinvG%*%(numSamples*diag(hatmat[hatmatRow,]))%*%t(QinvG))
  }))
  
  list(psi=est, Q=Q, D=eigenG$values, G=expressionData)
}
