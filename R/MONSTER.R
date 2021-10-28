monsterAnalysis <- setClass("monsterAnalysis", slots=c("tm","nullTM","numGenes","numSamples"))
setMethod("show","monsterAnalysis",function(object){monsterPrintMonsterAnalysis(object)})
setGeneric("get_tm", function(object) standardGeneric("get_tm"))
setMethod("get_tm", "monsterAnalysis", function(object){object@tm})

#' monsterPlotMonsterAnalysis
#'
#' plots the sum of squares of off diagonal mass (differential TF Involvement)
#'
#' @param x an object of class "monsterAnalysis"
#' @param ... further arguments passed to or from other methods.
#' @export
#' @return Plot of the dTFI for each TF against null distribution
#' @examples
#' data(yeast)
#' yeast$exp.cc[is.na(yeast$exp.cc)] <- mean(as.matrix(yeast$exp.cc),na.rm=TRUE)
#' design <- c(rep(1,25),rep(0,10),rep(NA,15))
#' monsterRes <- monster(yeast$exp.cc, design,
#' yeast$motif, nullPerms=10, numMaxCores=1)
#' monsterPlotMonsterAnalysis(monsterRes)
monsterPlotMonsterAnalysis <- function(x, ...){
  monsterdTFIPlot(x,...)
}
#' monsterPrintMonsterAnalysis
#'
#' summarizes the results of a MONSTER analysis
#'
#' @param x an object of class "monster"
#' @param ... further arguments passed to or from other methods.
#' @export
#' @return Description of transition matrices in object
#' @examples
#' data(yeast)
#' yeast$exp.cc[is.na(yeast$exp.cc)] <- mean(as.matrix(yeast$exp.cc),na.rm=TRUE)
#' design <- c(rep(1,25),rep(0,10),rep(NA,15))
#' monster(yeast$exp.cc,design,yeast$motif, nullPerms=10, numMaxCores=1)
monsterPrintMonsterAnalysis <- function(x, ...){
  cat("MONSTER object\n")
  cat(paste(x@numGenes, "genes\n"))
  cat(paste(x@numSamples[1],"baseline samples\n"))
  cat(paste(x@numSamples[2],"final samples\n"))
  cat(paste("Transition driven by", ncol(x@tm), "transcription factors\n"))
  cat(paste("Run with", length(x@nullTM), "randomized permutations.\n"))
}

#' MOdeling Network State Transitions from Expression and Regulatory data (MONSTER)
#'
#' This function runs the MONSTER algorithm.  Biological states are characterized by distinct patterns 
#' of gene expression that reflect each phenotype's active cellular processes. 
#' Driving these phenotypes are gene regulatory networks in which transcriptions factors control 
#' when and to what degree individual genes are expressed. Phenotypic transitions, such as those that 
#' occur when disease arises from healthy tissue, are associated with changes in these  networks. 
#' MONSTER is an approach to understanding these transitions. MONSTER models phenotypic-specific 
#' regulatory networks and then estimates a "transition matrix" that converts one state to another. 
#' By examining the properties of the transition matrix, we can gain insight into regulatory 
#' changes associated with phenotypic state transition.
#' Important note: the direct regulatory network observed from gene expression is currently
#' implemented as a regular correlation as opposed to the partial correlation described 
#' in the paper.
#' Citation: Schlauch, Daniel, et al. "Estimating drivers of cell state transitions using gene regulatory network models." 
#' BMC systems biology 11.1 (2017): 139. https://doi.org/10.1186/s12918-017-0517-y
#' @param expr Gene Expression dataset, can be matrix or data.frame of expression values or ExpressionSet. 
#' @param design Binary vector indicating case control partition. 1 for case and 0 for control.
#' @param motif Regulatory data.frame consisting of three columns.  For each row, a transcription factor (column 1) 
#' regulates a gene (column 2) with a defined strength (column 3), usually taken to be 0 or 1 
#' @param nullPerms number of random permutations to run (default 100).  Set to 0 to only 
#' calculate observed transition matrix. When mode is is 'buildNet' it randomly permutes the case and control expression
#' samples, if mode is 'regNet' it will randomly permute the case and control networks.
#' @param ni_method String to indicate algorithm method.  Must be one of "bere","pearson","cd","lda", or "wcd". Default is "bere"
#' @param ni.coefficient.cutoff numeric to specify a p-value cutoff at the network
#' inference step.  Default is NA, indicating inclusion of all coefficients.
#' @param numMaxCores requires doParallel, foreach.  Runs MONSTER in parallel computing 
#' environment.  Set to 1 to avoid parallelization, NA will take the default parallel pool in the computer.
#' @param outputDir character vector specifying a directory or path in which 
#' which to save MONSTER results, default is NA and results are not saved.
#' @param alphaw A weight parameter between 0 and 1 specifying proportion of weight 
#' to give to indirect compared to direct evidence. The default is 0.5 to give an equal weight to direct and indirect evidence.
#' @param mode A parameter telling whether to build the regulatory networks ('buildNet') or to use provided regulatory networks
#' ('regNet'). If set to 'regNet', then the parameters motif, ni_method, ni.coefficient.cutoff, and alphaw will be set to NA.
#' @export
#' @import doParallel
#' @import parallel
#' @import foreach
#' @importFrom methods new
#' @return An object of class "monsterAnalysis" containing results
#' 
#' @examples
#' # Example with the network reconstruction step
#' data(yeast)
#' design <- c(rep(0,20),rep(NA,10),rep(1,20))
#' yeast$exp.cc[is.na(yeast$exp.cc)] <- mean(as.matrix(yeast$exp.cc),na.rm=TRUE)
#' monsterRes <- monster(yeast$exp.cc[1:500,], design, yeast$motif, nullPerms=10, numMaxCores=1)
#' # Example with provided networks
#' \donttest{
#' pandaResult <- panda(pandaToyData$motif, pandaToyData$expression, pandaToyData$ppi)
#' case=get_regNet(pandaResult)
#' nelemReg=dim(get_regNet(pandaResult))[1]*dim(get_regNet(pandaResult))[2]
#' nGenes=length(colnames(get_regNet(pandaResult)))
#' control=matrix(rexp(nelemReg, rate=.1), ncol=nGenes)
#' colnames(control) = colnames(case)
#' rownames(control) = rownames(case) 
#' expr = as.data.frame(cbind(control,case))
#' design=c(rep(0,nGenes),rep(1, nGenes))
#' monsterRes <- monster(expr, design, motif=NA, nullPerms=10, numMaxCores=1, mode='regNet')
#' }

monster <- function(expr, 
                    design, 
                    motif, 
                    nullPerms=100,
                    ni_method="BERE",
                    ni.coefficient.cutoff = NA,
                    numMaxCores=1, 
                    outputDir=NA, alphaw=0.5, mode='buildNet'){
  
  if(mode=='regNet'){
    motif=NA
    alphaw=NA
    ni_method=NA
    ni.coefficient.cutoff=NA
    if(length(design == 1) != length(design == 0)){
      stop('case and control have a different number of genes')
    }
  }else{
    if(is.null(motif)){
      stop("motif may not be NULL")
    }
  }
  # Data type checking
  expr <- monsterCheckDataType(expr)
  # Parallelize
  # Initiate cluster
  if(!is.na(numMaxCores) && numMaxCores > 1){
    # Calculate the number of cores
    numCores <- detectCores() - 4
    numCores <- min(numCores, numMaxCores)
    cl <- makeCluster(numCores)
    registerDoParallel(cl)
    print("Running null permutations in parallel")
    print(paste(numCores,"cores used"))
  }
  
  iters <- nullPerms+1 # Two networks for each partition, plus observed partition
  print(paste(iters,"network transitions to be estimated"))
  
  #start time
  strt  <- Sys.time()
  #loop
  if(!is.na(outputDir)){
    dir.create(file.path(outputDir))  
    dir.create(file.path(outputDir,"tms"))        
  }
  
  # Remove unassigned data 
  expr <- expr[,design%in%c(0,1)]
  design <- design[design%in%c(0,1)]
  
  # Check column order
  if(mode == 'regNet'){
    numGenes = ncol(expr)/2
    if(any(colnames(expr[,design==1]) != colnames(expr[,design==0]))){
      stop('Please provide two regulatory networks with the same gene labels and 
           the same number of genes in the same order')
    }
  }else if(mode == 'buildNet'){
    numGenes = nrow(expr)
  }
  
  nullExpr <- expr
  if(numMaxCores == 1){
    transMatrices=list()
    for(i in seq_len(iters)){
      print(paste0("Running iteration ", i))
      if(i!=1){
        if(mode == 'regNet'){
          # Resample columns of provided network
          nullExpr[] <- expr[,sample(seq_along(colnames(expr)))]
        }else if(mode=='buildNet'){
          # Resample all entries in gene expression matrix then build null network
          nullExpr[] <- expr[sample(seq_along(c(expr)))]
      }
    }
      if(mode == 'buildNet'){
        nullExprCases <- nullExpr[,design==1]
        nullExprControls <- nullExpr[,design==0]
        
        tmpNetCases <- monsterMonsterNI(motif, nullExprCases, 
                                         method=ni_method, regularization="none",
                                         score="none", ni.coefficient.cutoff,
                                         verbose=TRUE, randomize = "none", cpp=FALSE,
                                         alphaw)
        tmpNetControls <- monsterMonsterNI(motif, nullExprControls, 
                                            method=ni_method, regularization="none",
                                            score="none", ni.coefficient.cutoff,
                                            verbose=TRUE, randomize = "none", cpp=FALSE,
                                            alphaw)
      }else if(mode == 'regNet'){
        tmpNetCases    = nullExpr[,design==1]
        tmpNetControls = nullExpr[,design==0]
      }
      transitionMatrix <- monsterTransformationMatrix(
        tmpNetControls, tmpNetCases, remove.diagonal=TRUE, method="ols")    
      print(paste("Finished running iteration", i))
      if (!is.na(outputDir)){
        saveRDS(transitionMatrix,file.path(outputDir,'tms',paste0('tm_',i,'.rds')))
      }
      transMatrices[[i]]=transitionMatrix
    }
    print(Sys.time()-strt)
  }else{
    transMatrices <- foreach(i=seq_len(iters),
                             .packages=c("netZooR","reshape2","penalized","MASS")) %dopar% {
                               print(paste0("Running iteration ", i))
                               if(i!=1){
                                 if(mode == 'regNet'){
                                   # Resample columns of provided network
                                   nullExpr[] <- expr[,sample(seq_along(colnames(expr)))]
                                 }else if(mode=='buildNet'){
                                   # Resample all entries in gene expression matrix then build null network
                                   nullExpr[] <- expr[sample(seq_along(c(expr)))]
                                 }
                               }
                               if(mode == 'buildNet'){
                                 nullExprCases <- nullExpr[,design==1]
                                 nullExprControls <- nullExpr[,design==0]
                                 
                                 tmpNetCases <- monsterMonsterNI(motif, nullExprCases, 
                                                                  method=ni_method, regularization="none",
                                                                  score="none", ni.coefficient.cutoff,
                                                                  verbose = FALSE, randomize = "none",
                                                                  alphaw)
                                 tmpNetControls <- monsterMonsterNI(motif, nullExprControls, 
                                                                     method=ni_method, regularization="none",
                                                                     score="none", ni.coefficient.cutoff,
                                                                     verbose = FALSE, randomize = "none",
                                                                     alphaw)
                               }else if(mode == 'regNet'){
                                 tmpNetCases    = nullExpr[,design==1]
                                 tmpNetControls = nullExpr[,design==0]
                               }
                               transitionMatrix <- monsterTransformationMatrix(
                                 tmpNetControls, tmpNetCases, remove.diagonal=TRUE, method="ols")    
                               print(paste("Finished running iteration", i))
                               if (!is.na(outputDir)){
                                 saveRDS(transitionMatrix,file.path(outputDir,'tms',paste0('tm_',i,'.rds')))
                               }
                               transitionMatrix
                             }
    
    print(Sys.time()-strt)
  }
  if(!is.na(numMaxCores)  && numMaxCores > 1){
    stopCluster(cl)
  }
  
  gc()
  return(
    monsterAnalysis(
      tm=transMatrices[[1]], 
      nullTM=transMatrices[-1], 
      numGenes=numGenes, 
      numSamples=c(sum(design==0), sum(design==1))))
}

#' Checks that data is something MONSTER can handle
#'
#' @param expr Gene Expression dataset
#' @return expr Gene Expression dataset in the proper form (may be the same as input)
#' @importFrom assertthat assert_that
#' @importFrom methods is
#' @export
#' @examples
#' expr.matrix <- matrix(rnorm(2000),ncol=20)
#' monsterCheckDataType(expr.matrix)
#' #TRUE
#' data(yeast)
#' class(yeast$exp.cc)
#' monsterCheckDataType(yeast$exp.cc)
#' #TRUE
monsterCheckDataType <- function(expr){
  assert_that(is.data.frame(expr)||is.matrix(expr)||is(expr,"ExpressionSet"))
  if("ExpressionSet" %in% class(expr)){
    if (requireNamespace("Biobase", quietly = TRUE)) {
      expr <- Biobase::exprs(expr)
    } 
  }
  if(is.data.frame(expr)){
    expr <- as.matrix(expr)
  }
  expr
}

globalVariables("i")

#' Bi-partite network analysis tools
#'
#' This function analyzes a bi-partite network.
#'
#' @param network.1 starting network, a genes by transcription factors data.frame with scores 
#' for the existence of edges between
#' @param network.2 final network, a genes by transcription factors data.frame with scores 
#' for the existence of edges between
#' @param by.tfs logical indicating a transcription factor based transformation.    If 
#' false, gives gene by gene transformation matrix
#' @param remove.diagonal logical for returning a result containing 0s across the diagonal
#' @param standardize logical indicating whether to standardize the rows and columns
#' @param method character specifying which algorithm to use, default='ols'
#' @return matrix object corresponding to transition matrix
#' @import MASS
#' @importFrom penalized optL1
#' @importFrom reshape2 melt
#' @export
#' @examples
#' data(yeast)
#' cc.net.1 <- monsterMonsterNI(yeast$motif,yeast$exp.cc[1:1000,1:20])
#' cc.net.2 <- monsterMonsterNI(yeast$motif,yeast$exp.cc[1:1000,31:50])
#' monsterTransformationMatrix(cc.net.1, cc.net.2)

monsterTransformationMatrix <- function(network.1, network.2, by.tfs=TRUE, standardize=FALSE, 
                                          remove.diagonal=TRUE, method="ols"){
  if(is.list(network.1)&&is.list(network.2)){
    if(by.tfs){
      net1 <- t(network.1$reg.net)
      net2 <- t(network.2$reg.net)
    } else {
      net1 <- network.1$reg.net
      net2 <- network.2$reg.net
    }
  } else if(is.matrix(network.1)&&is.matrix(network.2)){
    if(by.tfs){
      net1 <- t(network.1)
      net2 <- t(network.2)
    } else {
      net1 <- network.1
      net2 <- network.2
    }
  } else {
    stop("Networks must be lists or matrices")
  }
  
  if(!method%in%c("ols","kabsch","L1","orig")){
    stop("Invalid method.  Must be one of 'ols', 'kabsch', 'L1','orig'")
  }
  if (method == "kabsch"){
    tf.trans.matrix <- kabsch(net1,net2)
  }
  if (method == "orig"){
    svd.net2 <- svd(net2)
    tf.trans.matrix <- svd.net2$v %*% diag(1/svd.net2$d) %*% t(svd.net2$u) %*% net1
  }
  if (method == "ols"){
    net2.star <- vapply(seq_len(ncol(net1)), function(i,x,y){
      lm(y[,i]~x[,i])$resid
    }, net1, net2)
    tf.trans.matrix <- ginv(t(net1)%*%net1)%*%t(net1)%*%net2.star
    colnames(tf.trans.matrix) <- colnames(net1)
    rownames(tf.trans.matrix) <- colnames(net1)
    print("Using OLS method")
    
  }
  if (method == "L1"){
    net2.star <- vapply(seq_len(ncol(net1)), function(i,x,y){
      lm(y[,i]~x[,i])$resid
    }, net1, net2)
    tf.trans.matrix <- vapply(seq_len(ncol(net1)), function(i){
      z <- optL1(net2.star[,i], net1, fold=5, minlambda1=1, 
                 maxlambda1=2, model="linear", standardize=TRUE)
      coefficients(z$fullfit, "penalized")
    })
    colnames(tf.trans.matrix) <- rownames(tf.trans.matrix)
    print("Using L1 method")
    
  }
  if (standardize){
    tf.trans.matrix <- apply(tf.trans.matrix, 1, function(x){
      x/sum(abs(x))
    })
  }
  
  if (remove.diagonal){
    diag(tf.trans.matrix) <- 0
  }
  colnames(tf.trans.matrix) <- rownames(tf.trans.matrix)
  tf.trans.matrix
}

kabsch <- function(P,Q){
  
  P <- apply(P,2,function(x){
    x - mean(x)
  })
  Q <- apply(Q,2,function(x){
    x - mean(x)
  })
  covmat <- cov(P,Q)
  P.bar <- colMeans(P)
  Q.bar <- colMeans(Q)
  num.TFs <- ncol(P)        #n
  num.genes <- nrow(P)    #m
  
  #     covmat <- (t(P)%*%Q - P.bar%*%t(Q.bar)*(num.genes))
  
  svd.res <- svd(covmat-num.TFs*Q.bar%*%t(P.bar))
  
  # Note the scalar multiplier in the middle.
  # NOT A MISTAKE!
  c.k <- colSums(P %*% svd.res$v * Q %*% svd.res$u) - 
    num.genes*(P.bar%*%svd.res$v)*(Q.bar%*%svd.res$u)
  
  E <- diag(c(sign(c.k)))
  
  W <- svd.res$v %*% E %*% t(svd.res$u)
  rownames(W) <- colnames(P)
  colnames(W) <- colnames(P)
  W
}

#' Transformation matrix plot
#'
#' This function plots a hierachically clustered heatmap and 
#' corresponding dendrogram of a transaction matrix
#'
#' @param monsterObj monsterAnalysis Object
#' @param method distance metric for hierarchical clustering.    
#' Default is "Pearson correlation"
#' @export
#' @import ggplot2
#' @import grid
#' @rawNamespace import(stats, except= c(cov2cor,decompose,toeplitz,lowess,update,spectrum))
#' @return ggplot2 object for transition matrix heatmap
#' @examples
#' # data(yeast)
#' # design <- c(rep(0,20),rep(NA,10),rep(1,20))
#' # yeast$exp.cc[is.na(yeast$exp.cc)] <- mean(as.matrix(yeast$exp.cc),na.rm=TRUE)
#' # monsterRes <- monster(yeast$exp.cc, design, yeast$motif, nullPerms=10, numMaxCores=1)
#' data(monsterRes)
#' monsterHclHeatmapPlot(monsterRes)
monsterHclHeatmapPlot <- function(monsterObj, method="pearson"){
  x <- monsterObj@tm
  if(method=="pearson"){
    dist.func <- function(y) as.dist(cor(y))
  } else {
    dist.func <- dist
  }
  x <- scale(x)
  dd.col <- as.dendrogram(hclust(dist.func(x)))
  col.ord <- order.dendrogram(dd.col)
  
  dd.row <- as.dendrogram(hclust(dist.func(t(x))))
  row.ord <- order.dendrogram(dd.row)
  
  xx <- x[col.ord, row.ord]
  xx_names <- attr(xx, "dimnames")
  df <- as.data.frame(xx)
  colnames(df) <- xx_names[[2]]
  df$Var1 <- xx_names[[1]]
  df$Var1 <- with(df, factor(Var1, levels=Var1, ordered=TRUE))
  mdf <- melt(df)
  
  
  ddata_x <- dendro_data(dd.row)
  ddata_y <- dendro_data(dd.col)
  
  ### Set up a blank theme
  theme_none <- theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.title.x = element_text(colour=NA),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.line = element_blank()
  )
  ### Set up a blank theme
  theme_heatmap <- theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.title.x = element_text(colour=NA),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.line = element_blank()
  )
  ### Create plot components ###
  # Heatmap
  p1 <- ggplot(mdf, aes(x=variable, y=Var1)) +
    geom_tile(aes(fill=value)) + 
    scale_fill_gradient2() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  # Dendrogram 1
  p2 <- ggplot(segment(ddata_x)) +
    geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) +
    theme_none + theme(axis.title.x=element_blank())
  
  # Dendrogram 2
  p3 <- ggplot(segment(ddata_y)) +
    geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) +
    coord_flip() + theme_none
  
  ### Draw graphic ###
  
  grid.newpage()
  print(p1, vp=viewport(0.80, 0.8, x=0.400, y=0.40))
  print(p2, vp=viewport(0.73, 0.2, x=0.395, y=0.90))
  print(p3, vp=viewport(0.20, 0.8, x=0.910, y=0.43))
}

#' Principal Components plot of transformation matrix
#'
#' This function plots the first two principal components for a 
#' transaction matrix
#'
#' @param monsterObj a monsterAnalysis object resulting from a monster analysis
#' @param title The title of the plot
#' @param clusters A vector indicating the number of clusters to compute
#' @param alpha A vector indicating the level of transparency to be plotted
#' @return ggplot2 object for transition matrix PCA
#' @import ggdendro
#' @export
#' @examples
#' # data(yeast)
#' # design <- c(rep(0,20),rep(NA,10),rep(1,20))
#' # yeast$exp.cc[is.na(yeast$exp.cc)] <- mean(as.matrix(yeast$exp.cc),na.rm=TRUE)
#' # monsterRes <- monster(yeast$exp.cc, design, yeast$motif, nullPerms=100, numMaxCores=4)#' 
#' data(monsterRes)
#' # Color the nodes according to cluster membership
#' clusters <- kmeans(get_tm(monsterRes),3)$cluster 
#' monsterTransitionPCAPlot(monsterRes, 
#' title="PCA Plot of Transition - Cell Cycle vs Stress Response", 
#' clusters=clusters)
monsterTransitionPCAPlot <-    function(monsterObj, 
                                         title="PCA Plot of Transition", 
                                         clusters=1, alpha=1){
  tm.pca <- princomp(monsterObj@tm)
  odsm <- apply(monsterObj@tm,2,function(x){t(x)%*%x})
  odsm.scaled <- 2*(odsm-mean(odsm))/sd(odsm)+4
  scores.pca <- as.data.frame(tm.pca$scores)
  scores.pca <- cbind(scores.pca,'node.names'=rownames(scores.pca))
  ggplot(data = scores.pca, aes(x = Comp.1, y = Comp.2, label = node.names)) +
    geom_hline(yintercept = 0, colour = "gray65") +
    geom_vline(xintercept = 0, colour = "gray65") +
    geom_text(size = odsm.scaled, alpha=alpha, color=clusters) +
    ggtitle(title)
}

#' This function uses igraph to plot the transition matrix (directed graph) as a network.
#' The edges in the network should be read as A 'positively/negatively contributes to' the 
#' targeting of B in the target state.
#'
#' @param monsterObj monsterAnalysis Object
#' @param numEdges The number of edges to display
#' @param numTopTFs The number of TFs to display, only when rescale='significance'
#' @param rescale string to specify the order of edges. If set to 'significance', 
#' the TFs with the largest dTFI significance (smallest dTFI p-values) will be filtered first before
#' plotting the edges with the largest magnitude in the transition matrix. Otherwise
#' the filtering step will be skipped and the edges with the largest transitions will be plotted.
#' The plotted graph represents the top numEdges edges between the numTopTFs if rescale=='significance'
#' and top numEdges edges otherwise. The edge weight represents the observed transition edges standardized
#' by the null and the node size in the graph is proportional to the p-values of the dTFIs of each
#' TF. When rescale is set to 'significance', the results can be different between two MONSTER runs
#' if the number of permutations is not large enough to sample the null, that is why it is the seed should be set
#' prior to calling MONSTER to get reproducible results. If rescale is set to another value such as 'none', it will
#' produce deterministic results between two identical MONSTER runs.
#' @importFrom igraph graph.data.frame plot.igraph V E V<- E<-
#' @export
#' @return plot the transition matrix (directed graph) as a network.
#' @examples
#' # data(yeast)
#' # yeast$exp.cc[is.na(yeast$exp.cc)] <- mean(as.matrix(yeast$exp.cc),na.rm=TRUE)
#' # design <- c(rep(0,20),rep(NA,10),rep(1,20))
#' # monsterRes <- monster(yeast$exp.cc, design, yeast$motif, nullPerms=100, numMaxCores=4)#' 
#' data(monsterRes)
#' monstertransitionNetworkPlot(monsterRes, rescale='significance')
#' monstertransitionNetworkPlot(monsterRes, rescale='none')

monstertransitionNetworkPlot <- function(monsterObj, numEdges=100, numTopTFs=10, rescale='significance'){
  ## Calculate p-values for off-diagonals
  transitionSigmas <- function(tm.observed, tm.null){
    tm.null.mean <- apply(simplify2array(tm.null), seq_len(2), mean)
    tm.null.sd <- apply(simplify2array(tm.null), seq_len(2), sd)
    sigmas <- (tm.observed - tm.null.mean)/tm.null.sd
  }
  
  tm.sigmas <- transitionSigmas(monsterObj@tm, monsterObj@nullTM)
  diag(tm.sigmas) <- 0
  tm.sigmas.melt <- melt(tm.sigmas)
  
  adjMat <- monsterObj@tm
  diag(adjMat) <- 0
  adjMat.melt <- melt(adjMat)
  
  adj.combined <- merge(tm.sigmas.melt, adjMat.melt, by=c("Var1","Var2"))
  
  # adj.combined[,1] <- mappings[match(adj.combined[,1], mappings[,1]),2]
  # adj.combined[,2] <- mappings[match(adj.combined[,2], mappings[,1]),2]
  dTFI_pVals_All <- 1-2*abs(.5-monsterCalculateTmPvalues(monsterObj, 
                                                             method="z-score"))
  if(rescale=='significance'){
    topTFsIncluded <- names(sort(dTFI_pVals_All)[seq_len(numTopTFs)])
    topTFIndices <- 2>(is.na(match(adj.combined[,1],topTFsIncluded)) + 
                         is.na(match(adj.combined[,2],topTFsIncluded)))
    adj.combined <- adj.combined[topTFIndices,]
  }
  
  adj.combined <- adj.combined[
    abs(adj.combined[,4])>=sort(abs(adj.combined[,4]),decreasing=TRUE)[numEdges],]
  tfNet <- graph.data.frame(adj.combined, directed=TRUE)
  vSize <- -log(dTFI_pVals_All)
  vSize[vSize<0] <- 0
  vSize[vSize>3] <- 3
  
  V(tfNet)$size <- vSize[V(tfNet)$name]*5
  V(tfNet)$color <- "yellow"
  E(tfNet)$width <- (abs(E(tfNet)$value.x))*15/max(abs(E(tfNet)$value.x))
  E(tfNet)$color <-ifelse(E(tfNet)$value.x>0, "blue", "red")
  
  plot.igraph(tfNet, edge.arrow.size=2, vertex.label.cex= 1.5, vertex.label.color= "black",main="")
}

#' This function plots the Off diagonal mass of an 
#' observed Transition Matrix compared to a set of null TMs
#'
#' @param monsterObj monsterAnalysis Object
#' @param rescale string indicating whether to reorder transcription
#' factors according to their statistical significance and to 
#' rescale the values observed to be standardized by the null
#' distribution ('significance'), to reorder transcription
#' factors according to the largest dTFIs ('magnitude') with the TF x axis labels proportional to their significance
#' , or finally without ordering them ('none'). When rescale is set to 'significance', 
#' the results can be different between two MONSTER runs if the number of permutations is not large enough to sample 
#' the null, that is why it is the seed should be set prior to calling MONSTER to get reproducible results. 
#' If rescale is set to another value such as 'magnitude' or 'none', it will produce deterministic results 
#' between two identical MONSTER runs.
#' @param plot.title String specifying the plot title
#' @param highlight.tfs vector specifying a set of transcription 
#' factors to highlight in the plot
#' @param nTFs number of TFs to plot in x axis. -1 takes all TFs.
#' @return ggplot2 object for transition matrix comparing observed 
#' distribution to that estimated under the null 
#' @export
#' @examples
#' # data(yeast)
#' # yeast$exp.cc[is.na(yeast$exp.cc)] <- mean(as.matrix(yeast$exp.cc),na.rm=TRUE)
#' # design <- c(rep(0,20),rep(NA,10),rep(1,20))
#' # monsterRes <- monster(yeast$exp.cc, design, yeast$motif, nullPerms=100, numMaxCores=4)#' 
#' data(monsterRes)
#' monsterdTFIPlot(monsterRes)
monsterdTFIPlot <- function(monsterObj, rescale='none', plot.title=NA, highlight.tfs=NA,
                             nTFs=-1){
  if(is.na(plot.title)){
    plot.title <- "Differential TF Involvement"
  }
  num.iterations <- length(monsterObj@nullTM)
  # Calculate the off-diagonal squared mass for each transition matrix
  null.SSODM <- lapply(monsterObj@nullTM,function(x){
    apply(x,2,function(y){t(y)%*%y})
  })
  null.ssodm.matrix <- matrix(unlist(null.SSODM),ncol=num.iterations)
  null.ssodm.matrix <- t(apply(null.ssodm.matrix,1,sort))
  
  ssodm <- apply(monsterObj@tm,2,function(x){t(x)%*%x})
  
  p.values <- 1-pnorm(vapply(seq_along(ssodm),function(i){
    (ssodm[i]-mean(null.ssodm.matrix[i,]))/sd(null.ssodm.matrix[i,])
  }))
  t.values <- vapply(seq_along(ssodm),function(i){
    (ssodm[i]-mean(null.ssodm.matrix[i,]))/sd(null.ssodm.matrix[i,])
  })
  
  # Process the data for ggplot2
  combined.mat <- cbind(null.ssodm.matrix, ssodm)
  colnames(combined.mat) <- c(rep('Null',num.iterations),"Observed")
  
  
  if (rescale == 'significance'){
    combined.mat <- t(apply(combined.mat,1,function(x){
      (x-mean(x[-(num.iterations+1)]))/sd(x[-(num.iterations+1)])
    }))
    x.axis.order <- rownames(monsterObj@nullTM[[1]])[order(-t.values)]
    x.axis.size    <- 10 # pmin(15,7-log(p.values[order(p.values)]))
  } else if (rescale == 'none'){
    x.axis.order <- rownames(monsterObj@nullTM[[1]])
    x.axis.size    <- pmin(15,7-log(p.values))
  } else if (rescale == 'magnitude'){
    x.axis.order <- rownames(monsterObj@nullTM[[1]])[order(-combined.mat[, dim(combined.mat)[2]])]
    x.axis.size    <- pmin(15,7-log(p.values))
  }
  if(nTFs==-1){
    nTFs = length(x.axis.order)
  }
  null.SSODM.melt <- melt(combined.mat)[,-1][,c(2,1)]
  null.SSODM.melt$TF<-rep(rownames(monsterObj@nullTM[[1]]),num.iterations+1)
  
  ## Plot the data
  ggplot(null.SSODM.melt, aes(x=TF, y=value))+
    geom_point(aes(color=factor(Var2), alpha = .5*as.numeric(factor(Var2))), size=2) +
    scale_color_manual(values = c("blue", "red")) +
    scale_alpha(guide = "none") +
    scale_x_discrete(limits = x.axis.order[seq_len(nTFs)] ) +
    theme_classic() +
    theme(legend.title=element_blank(),
          axis.text.x = element_text(colour = 1+x.axis.order%in%highlight.tfs, 
                                     angle = 90, hjust = 1, 
                                     size=x.axis.size,face="bold")) +
    ylab("dTFI") +
    ggtitle(plot.title)
  
}

#' Calculate p-values for a tranformation matrix
#'
#' This function calculates the significance of an observed
#' transition matrix given a set of null transition matrices
#'
#' @param monsterObj monsterAnalysis Object
#' @param method one of 'z-score' or 'non-parametric'
#' @return vector of p-values for each transcription factor
#' @export
#' @examples
#' # data(yeast)
#' # design <- c(rep(0,20),rep(NA,10),rep(1,20))
#' # yeast$exp.cc[is.na(yeast$exp.cc)] <- mean(as.matrix(yeast$exp.cc),na.rm=TRUE)
#' # monsterRes <- monster(yeast$exp.cc, design, yeast$motif, nullPerms=100, numMaxCores=4)
#' data(monsterRes)
#' monsterCalculateTmPvalues(monsterRes)
monsterCalculateTmPvalues <- function(monsterObj, method="z-score"){
  num.iterations <- length(monsterObj@nullTM)
  # Calculate the off-diagonal squared mass for each transition matrix
  null.SSODM <- lapply(monsterObj@nullTM,function(x){
    apply(x,1,function(y){t(y)%*%y})
  })
  null.ssodm.matrix <- matrix(unlist(null.SSODM),ncol=num.iterations)
  null.ssodm.matrix <- t(apply(null.ssodm.matrix,1,sort))
  
  ssodm <- apply(monsterObj@tm,1,function(x){t(x)%*%x})
  
  # Get p-value (rank of observed within null ssodm)
  if(method=="non-parametric"){
    p.values <- vapply(seq_along(ssodm),function(i){
      1-findInterval(ssodm[i], null.ssodm.matrix[i,])/num.iterations
    })
  } else if (method=="z-score"){
    p.values <- pnorm(vapply(seq_along(ssodm),function(i){
      (ssodm[i]-mean(null.ssodm.matrix[i,]))/sd(null.ssodm.matrix[i,])
    }))
  } else {
    print('Undefined method')
  }
  p.values
}

globalVariables(c("Var1", "Var2","value","variable","xend","yend","y","Comp.1", "Comp.2","node.names","TF","i"))

#' Bipartite Edge Reconstruction from Expression data
#'
#' This function generates a complete bipartite network from 
#' gene expression data and sequence motif data
#' 
#' @param motif.data A motif dataset, a data.frame, matrix or exprSet containing 
#' 3 columns. Each row describes an motif associated with a transcription 
#' factor (column 1) a gene (column 2) and a score (column 3) for the motif.
#' @param expr.data An expression dataset, as a genes (rows) by samples (columns)
#' @param verbose logical to indicate printing of output for algorithm progress.
#' @param method String to indicate algorithm method.  Must be one of 
#' "bere","pearson","cd","lda", or "wcd". Default is "bere".
#' Important note: the direct regulatory network observed from gene expression is currently
#' implemented as a regular correlation as opposed to the partial correlation described 
#' in the paper (please see Schlauch et al., 2017,  https://doi.org/10.1186/s12918-017-0517-y)
#' @param ni.coefficient.cutoff numeric to specify a p-value cutoff at the network
#' inference step.  Default is NA, indicating inclusion of all coefficients.
#' @param randomize logical indicating randomization by genes, within genes or none
#' @param score String to indicate whether motif information will be 
#' readded upon completion of the algorithm
#' to give to indirect compared to direct evidence.  See documentation.
#' @param regularization String parameter indicating one of "none", "L1", "L2"
#' @param cpp logical use C++ for maximum speed, set to false if unable to run.
#' @export
#' @return matrix for inferred network between TFs and genes
#' @importFrom tidyr spread
#' @importFrom penalized penalized
#' @importFrom reshape2 dcast
#' @examples
#' data(yeast)
#' cc.net <- monsterMonsterNI(yeast$motif,yeast$exp.cc)

monsterMonsterNI <- function(motif.data, 
                              expr.data,
                              verbose=FALSE,
                              randomize="none",
                              method="bere",
                              ni.coefficient.cutoff=NA,
                              alphaw=1.0,
                              regularization="none",
                              score="motifincluded",
                              cpp=FALSE){
  if(verbose)
    print('Initializing and validating')
  # Create vectors for TF names and Gene names from Motif dataset
  tf.names   <- sort(unique(motif.data[,1]))
  num.TFs    <- length(tf.names)
  if (is.null(expr.data)){
    stop("Expression data null")
  } else {
    # Use the motif data AND the expr data (if provided) for the gene list
    gene.names <- sort(intersect(motif.data[,2],rownames(expr.data)))
    num.genes  <- length(gene.names)
    
    # Filter out the expr genes without motif data
    expr.data <- expr.data[rownames(expr.data) %in% gene.names,]
    
    # Keep everything sorted alphabetically
    expr.data      <- expr.data[order(rownames(expr.data)),]
    num.conditions <- ncol(expr.data);
    if (randomize=='within.gene'){
      expr.data <- t(apply(expr.data, 1, sample))
      if(verbose)
        print("Randomizing by reordering each gene's expression")
    } else if (randomize=='by.genes'){
      rownames(expr.data) <- sample(rownames(expr.data))
      expr.data           <- expr.data[order(rownames(expr.data)),]
      if(verbose)
        print("Randomizing by reordering each gene labels")
    }
  }
  
  # Bad data checking
  if (num.genes==0){
    stop("Validating data.  No matched genes.\n
            Please ensure that gene names in expression 
            file match gene names in motif file.")
  }
  
  strt<-Sys.time()
  if(num.conditions==0) {
    stop("Number of samples = 0")
    gene.coreg <- diag(num.genes)
  } else if(num.conditions<3) {
    stop('Not enough expression conditions detected to calculate correlation.')
  } else {
    if(verbose)
      print('Verified adequate samples, calculating correlation matrix')
    if(cpp){
      # C++ implementation
      gene.coreg <- rcpp_ccorr(t(apply(expr.data, 1, function(x)(x-mean(x))/(sd(x)))))
      rownames(gene.coreg)<- rownames(expr.data)
      colnames(gene.coreg)<- rownames(expr.data)
      
    } else if(!(method %in% c("BERE","pearson"))) {
      # Standard r correlation calculation
      gene.coreg <- cor(t(expr.data), method="pearson", use="pairwise.complete.obs")
    }
  }
  
  print(Sys.time()-strt)
  
  if(verbose)
    print('More data cleaning')
  # Convert 3 column format to matrix format
  colnames(motif.data) <- c('TF','GENE','value')
  if( method != "BERE"){
    regulatory.network <- spread(motif.data, GENE, value, fill=0)
    rownames(regulatory.network) <- regulatory.network[,1]
    # sort the TFs (rows), and remove redundant first column
    regulatory.network <- regulatory.network[order(rownames(regulatory.network)),-1]
    # sort the genes (columns)
    regulatory.network <- as.matrix(regulatory.network[,order(colnames(regulatory.network))])
    
    # Filter out any motifs that are not in expr dataset (if given)
    if (!is.null(expr.data)){
      regulatory.network <- regulatory.network[,colnames(regulatory.network) %in% gene.names]
    }
    
    # store initial motif network (alphabetized for rows and columns)
    #   starting.motifs <- regulatory.network
  }
  
  if(verbose)
    print('Main calculation')
  result <- NULL
  ########################################
  if (method=="BERE"){
    
    expr.data <- data.frame(expr.data)
    tfdcast <- dcast(motif.data,TF~GENE,fill=0)
    rownames(tfdcast) <- tfdcast[,1]
    tfdcast <- tfdcast[,-1]
    
    expr.data <- expr.data[sort(rownames(expr.data)),]
    tfdcast <- tfdcast[,sort(colnames(tfdcast)),]
    tfNames <- rownames(tfdcast)[rownames(tfdcast) %in% rownames(expr.data)]
    
    ## Filtering
    # filter out the TFs that are not in expression set
    tfdcast <- tfdcast[rownames(tfdcast)%in%tfNames,]
    
    # Filter out genes that aren't targetted by anything 7/28/15
    commonGenes <- intersect(colnames(tfdcast),rownames(expr.data))
    expr.data <- expr.data[commonGenes,]
    tfdcast <- tfdcast[,commonGenes]
    
    # check that IDs match
    if (prod(rownames(expr.data)==colnames(tfdcast))!=1){
      stop("ID mismatch")
    }
    
    ## Get direct evidence
    if ((1-alphaw)!=0){
      directCor <- t(cor(t(expr.data),t(expr.data[rownames(expr.data)%in%tfNames,]))^2)
    }else{
      directCor = matrix(0L, length(tfNames), length(commonGenes))
    }
    
    ## Get the indirect evidence    
    if(alphaw==0){
      result = matrix(0L, length(tfNames), length(commonGenes))
    }else{
      result <- t(apply(tfdcast, 1, function(x){
        cat(".")
        tfTargets <- as.numeric(x)
        z <- NULL
        if(regularization=="none"){
          z <- glm(tfTargets ~ ., data=expr.data, family="binomial")
          
          # 9/10/17
          # Adding argument to allow cutoffs based on p-values
          if(is.numeric(ni.coefficient.cutoff)){
            coefs <- coef(z)
            coefs[summary(z)$coef[,4]>ni.coefficient.cutoff] <- 0
            logit.res <- apply(expr.data,1,function(x){coefs[1] + sum(coefs[-1]*x)})
            return(exp(logit.res)/(1+exp(logit.res)))
            
          } else {
            return(predict(z, expr.data,type='response'))
          }
          
        } else {
          z <- penalized(tfTargets, expr.data,
                         lambda2=10, model="logistic", standardize=TRUE)
          # z <- optL1(tfTargets, expr.data, minlambda1=25, fold=5)
          
        }
        
        # Penalized Logistic Reg
        
        
        predict(z, expr.data)
      }))
    }
    
    ## Convert values to ranks
    if(alphaw<1 && alphaw>0){
      directCor <- matrix(rank(directCor), ncol=ncol(directCor))
      result <- matrix(rank(result), ncol=ncol(result))
    }
    consensus <- directCor*(1-alphaw) + result*alphaw
    rownames(consensus) <- rownames(tfdcast)
    colnames(consensus) <- rownames(expr.data)
    consensusRange <- max(consensus)- min(consensus)
    if(score=="motifincluded"){
      consensus <- as.matrix(consensus + consensusRange*regulatory.network)
    }
    result=consensus
  } else if (method=="pearson"){
    tfNames = levels(motif.data$TF)
    result <- t(cor(t(expr.data),t(expr.data[rownames(expr.data)%in%tfNames,]))^2)
    if(score=="motifincluded"){
      result <- as.matrix(consensus + consensusRange*regulatory.network)
    }
    result
  } else {
    strt<-Sys.time()
    # Remove NA correlations
    gene.coreg[is.na(gene.coreg)] <- 0
    correlation.dif <- sweep(regulatory.network,1,rowSums(regulatory.network),`/`)%*%
      gene.coreg - 
      sweep(1-regulatory.network,1,rowSums(1-regulatory.network),`/`)%*%
      gene.coreg
    result <- sweep(correlation.dif, 2, apply(correlation.dif, 2, sd),'/')
    #   regulatory.network <- ifelse(res>quantile(res,1-mean(regulatory.network)),1,0)
    
    print(Sys.time()-strt)
    ########################################
    if(score=="motifincluded"){
      result <- result + max(result)*regulatory.network
    }
    result
  }
  return(result)
}

#' Bipartite Edge Reconstruction from Expression data 
#' (composite method with direct/indirect)
#'
#' This function generates a complete bipartite network from 
#' gene expression data and sequence motif data. This NI method
#' serves as a default method for inferring bipartite networks
#' in MONSTER.  Running monsterBereFull can generate these networks
#' independently from the larger MONSTER method.
#'
#' @param motif.data A motif dataset, a data.frame, matrix or exprSet 
#' containing 3 columns. Each row describes an motif associated 
#' with a transcription factor (column 1) a gene (column 2) 
#' and a score (column 3) for the motif.
#' @param expr.data An expression dataset, as a genes (rows) by 
#' samples (columns) data.frame
#' @param alpha A weight parameter specifying proportion of weight 
#' to give to indirect compared to direct evidence.  See documentation.
#' @param lambda if using penalized, the lambda parameter in the penalized logistic regression
#' @param score String to indicate whether motif information will 
#' be readded upon completion of the algorithm
#' @importFrom reshape2 dcast
#' @importFrom penalized predict
#' @export
#' @return An matrix or data.frame
#' @examples
#' data(yeast)
#' monsterRes <- monsterBereFull(yeast$motif, yeast$exp.cc, alpha=.5)
monsterBereFull <- function(motif.data, 
                             expr.data, 
                             alpha=.5, 
                             lambda=10, 
                             score="motifincluded"){
  
  expr.data <- data.frame(expr.data)
  tfdcast <- dcast(motif.data,TF~GENE,fill=0)
  rownames(tfdcast) <- tfdcast[,1]
  tfdcast <- tfdcast[,-1]
  
  expr.data <- expr.data[sort(rownames(expr.data)),]
  tfdcast <- tfdcast[,sort(colnames(tfdcast)),]
  tfNames <- rownames(tfdcast)[rownames(tfdcast) %in% rownames(expr.data)]
  
  ## Filtering
  # filter out the TFs that are not in expression set
  tfdcast <- tfdcast[rownames(tfdcast)%in%tfNames,]
  
  # Filter out genes that aren't targetted by anything 7/28/15
  commonGenes <- intersect(colnames(tfdcast),rownames(expr.data))
  expr.data <- expr.data[commonGenes,]
  tfdcast <- tfdcast[,commonGenes]
  
  # check that IDs match
  if (prod(rownames(expr.data)==colnames(tfdcast))!=1){
    stop("ID mismatch")
  }
  ## Get direct evidence
  directCor <- t(cor(t(expr.data),t(expr.data[rownames(expr.data)%in%tfNames,]))^2)
  
  ## Get the indirect evidence    
  result <- t(apply(tfdcast, 1, function(x){
    cat(".")
    tfTargets <- as.numeric(x)
    
    # Ordinary Logistic Reg
    # z <- glm(tfTargets ~ ., data=expr.data, family="binomial")
    
    # Penalized Logistic Reg
    expr.data[is.na(expr.data)] <- 0
    z <- penalized(response=tfTargets, penalized=expr.data, unpenalized=~0,
                   lambda2=lambda, model="logistic", standardize=TRUE)
    #z <- optL1(tfTargets, expr.data, minlambda1=25, fold=5)
    
    
    predict(z, expr.data)
  }))
  
  ## Convert values to ranks
  directCor <- matrix(rank(directCor), ncol=ncol(directCor))
  result <- matrix(rank(result), ncol=ncol(result))
  
  consensus <- directCor*(1-alpha) + result*alpha
  rownames(consensus) <- rownames(tfdcast)
  colnames(consensus) <- rownames(expr.data)
  consensusRange <- max(consensus)- min(consensus)
  if(score=="motifincluded"){
    consensus <- as.matrix(consensus + consensusRange*tfdcast)
  }
  consensus
}

globalVariables(c("expr.data","lambda","rcpp_ccorr","GENE", "TF","value"))

#' MONSTER results from example cell-cycle yeast transition
#'
#'This data contains the MONSTER result from analysis of Yeast Cell cycle, included in data(yeast).  
#'This result arbitrarily takes the first 20 gene expression samples in yeast$cc to be the baseline condition, and the final 20 samples to be the final condition.
#'
#' @docType data
#' @keywords datasets
#' @name monsterRes
#' @usage data(monsterRes)
#' @format MONSTER obj
NULL

#' Toy data derived from three gene expression datasets and a mapping from transcription factors to genes.
#'
#'This data is a list containing gene expression data from three separate yeast studies along with data mapping yeast transcription factors with genes based on the presence of a sequence binding motif for each transcription factor in the vicinity of each gene. 
#' The motif data.frame, yeast$motif, describes a set of pairwise connections where a specific known sequence motif of a transcription factor was found upstream of the corresponding gene.   
#' The expression data, yeast$exp.ko, yeast$exp.cc, and yeast$exp.sr, are three gene expression datasets measured in conditions of gene knockout, cell cycle, and stress response, respectively. 
#' @docType data
#' @keywords datasets
#' @name yeast
#' @usage data(yeast)
#' @format A list containing 4 data.frames
#' @return A list of length 4
#' @references Glass K, Huttenhower C, Quackenbush J, Yuan GC. Passing Messages Between Biological Networks to Refine Predicted Interactions. PLoS One. 2013 May 31;8(5):e64832.
#' 
NULL
