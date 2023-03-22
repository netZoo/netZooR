#' Seeding PANDA Interactions to Derive Epigenetic Regulation
#'
#' This function runs the SPIDER algorithm
#'
#' @param motif A motif dataset, a data.frame, matrix or exprSet containing 3 columns.
#' Each row describes an motif associated with a transcription factor (column 1) a
#' gene (column 2) and a score (column 3) for the motif.
#' @param epifilter A binary matrix that is of the same size as motif that will be used as a mask to filter motif
#' for open chromatin region. Motif interactions that fall in open chromatin region will be kept and the others are removed.
#' @param expr An expression dataset, as a genes (rows) by samples (columns) data.frame
#' @param ppi A Protein-Protein interaction dataset, a data.frame containing 3 columns.
#' Each row describes a protein-protein interaction between transcription factor 1(column 1),
#' transcription factor 2 (column 2) and a score (column 3) for the interaction.
#' @param alpha value to be used for update variable, alpha (default=0.1)
#' @param hamming value at which to terminate the process based on hamming distance (default 10^-3)
#' @param iter sets the maximum number of iterations SPIDER can run before exiting.
#' @param progress Boolean to indicate printing of output for algorithm progress.
#' @param output a vector containing which networks to return.  Options include "regulatory",
#' "coregulatory", "cooperative".
#' @param zScale Boolean to indicate use of z-scores in output.  False will use [0,1] scale.
#' @param randomize method by which to randomize gene expression matrix.  Default "None".  Must
#' be one of "None", "within.gene", "by.genes".  "within.gene" randomization scrambles each row
#' of the gene expression matrix, "by.gene" scrambles gene labels.
#' @param cor.method Correlation method, default is "pearson".
#' @param scale.by.present Boolean to indicate scaling of correlations by percentage of positive samples.
#' @param remove.missing.ppi Boolean to indicate whether TFs in the PPI but not in the motif data should be
#' removed. Only when mode=='legacy'.
#' @param remove.missing.motif Boolean to indicate whether genes targeted in the motif data but not the
#' expression data should be removed. Only when mode=='legacy'.
#' @param remove.missing.genes Boolean to indicate whether genes in the expression data but lacking
#' information from the motif prior should be removed. Only when mode=='legacy'.
#' @param edgelist Boolean to indicate if edge lists instead of matrices should be returned. 
#' @param mode The data alignment mode. The mode 'union' takes the union of the genes in the expression matrix and the motif
#' and the union of TFs in the ppi and motif and fills the matrics with zeros for nonintersecting TFs and gens, 'intersection' 
#' takes the intersection of genes and TFs and removes nonintersecting sets, 'legacy' is the old behavior with PANDAR version 1.19.3.
#' #' Parameters remove.missing.ppi, remove.missingmotif, remove.missing.genes work only with mode=='legacy'.
#' @keywords keywords
#' @importFrom matrixStats rowSds
#' @importFrom matrixStats colSds
#' @importFrom Biobase assayData
#' @importFrom reshape melt.array
#' @export
#' @return An object of class "panda" containing matrices describing networks achieved by convergence
#' with SPIDER algorithm.\cr
#' "regNet" is the regulatory network\cr
#' "coregNet" is the coregulatory network\cr
#' "coopNet" is the cooperative network
#' @examples
#' data(pandaToyData)
#' pandaToyData$epifilter = pandaToyData$motif
#' nind=floor(runif(5000, min=1, max=dim(pandaToyData$epifilter)[1]))
#' pandaToyData$epifilter[nind,3] = 0
#' spiderRes <- spider(pandaToyData$motif,pandaToyData$expression,
#'                     pandaToyData$epifilter,pandaToyData$ppi,hamming=.1,progress=TRUE)
#' @references
#' Sonawane, Abhijeet Rajendra, et al. "Constructing gene regulatory networks using epigenetic data." npj Systems Biology and Applications 7.1 (2021): 1-13.
spider <- function(motif,expr=NULL,epifilter=NULL,ppi=NULL,alpha=0.1,hamming=0.001,
                  iter=NA,output=c('regulatory','coexpression','cooperative'),
                  zScale=TRUE,progress=FALSE,randomize=c("None", "within.gene", "by.gene"),cor.method="pearson",
                  scale.by.present=FALSE,edgelist=FALSE,remove.missing.ppi=FALSE,
                  remove.missing.motif=FALSE,remove.missing.genes=FALSE,mode="union"){
  
  randomize <- match.arg(randomize)  
  if(progress)
    print('Initializing and validating')
  
  if(any(epifilter[,c(1,2)] != motif[,c(1,2)])){
    stop('Chromatin accessibility data does not match motif data size and order.')
  }
  
  if(class(expr)=="ExpressionSet")
    expr <- assayData(expr)[["exprs"]]
  
  if (is.null(expr)){
    # Use only the motif data here for the gene list
    num.conditions <- 0
    if (randomize!="None"){
      warning("Randomization ignored because gene expression is not used.")
      randomize <- "None"
    }
  } else {
    if(mode=='legacy'){
      if(remove.missing.genes){
        # remove genes from expression data that are not in the motif data
        n <- nrow(expr)
        expr <- expr[which(rownames(expr)%in%motif[,2]),]
        message(sprintf("%s genes removed that were not present in motif", n-nrow(expr)))
      }
      if(remove.missing.motif){
        # remove genes from motif data that are not in the expression data
        n <- nrow(motif)
        motif <- motif[which(motif[,2]%in%rownames(expr)),]
        epifilter <- epifilter[which(motif[,2]%in%rownames(expr)),]
        message(sprintf("%s motif edges removed that targeted genes missing in expression data", n-nrow(motif)))
      }
      # Use the motif data AND the expr data (if provided) for the gene list
      # Keep everything sorted alphabetically
      expr <- expr[order(rownames(expr)),]
    }else if(mode=='union'){
      gene.names=unique(union(rownames(expr),unique(motif[,2])))
      tf.names  =unique(union(unique(ppi[,1]),unique(motif[,1])))
      num.TFs    <- length(tf.names)
      num.genes  <- length(gene.names)
      # gene expression matrix
      expr1=as.data.frame(matrix(0,num.genes,ncol(expr)))
      rownames(expr1)=gene.names
      expr1[which(gene.names%in%rownames(expr)),]=expr[]
      expr=expr1
      #PPI matrix
      tfCoopNetwork <- matrix(0,num.TFs,num.TFs)
      colnames(tfCoopNetwork)=tf.names
      rownames(tfCoopNetwork)=tf.names
      Idx1 <- match(ppi[,1], tf.names);
      Idx2 <- match(ppi[,2], tf.names);
      Idx <- (Idx2-1)*num.TFs+Idx1;
      tfCoopNetwork[Idx] <- ppi[,3];
      Idx <- (Idx1-1)*num.TFs+Idx2;
      tfCoopNetwork[Idx] <- ppi[,3];
      #Motif matrix
      regulatoryNetwork=matrix(0,num.TFs,num.genes)
      colnames(regulatoryNetwork)=gene.names
      rownames(regulatoryNetwork)=tf.names
      Idx1=match(motif[,1], tf.names);
      Idx2=match(motif[,2], gene.names);
      Idx=(Idx2-1)*num.TFs+Idx1;
      regulatoryNetwork[Idx]=motif[,3]*epifilter[,3]
    }else if(mode=='intersection'){
      gene.names=unique(intersect(rownames(expr),unique(motif[,2])))
      tf.names  =unique(intersect(unique(ppi[,1]),unique(motif[,1])))
      num.TFs    <- length(tf.names)
      num.genes  <- length(gene.names)
      # gene expression matrix
      expr1=as.data.frame(matrix(0,num.genes,ncol(expr)))
      rownames(expr1)=gene.names
      interGeneNames=gene.names[which(gene.names%in%rownames(expr))]
      expr1[interGeneNames,]=expr[interGeneNames,]
      expr=expr1
      #PPI matrix
      tfCoopNetwork <- matrix(0,num.TFs,num.TFs)
      colnames(tfCoopNetwork)=tf.names
      rownames(tfCoopNetwork)=tf.names
      Idx1 <- match(ppi[,1], tf.names);
      Idx2 <- match(ppi[,2], tf.names);
      Idx <- (Idx2-1)*num.TFs+Idx1;
      indIdx=!is.na(Idx)
      Idx=Idx[indIdx] #remove missing TFs
      tfCoopNetwork[Idx] <- ppi[indIdx,3];
      Idx <- (Idx1-1)*num.TFs+Idx2;
      indIdx=!is.na(Idx)
      Idx=Idx[indIdx] #remove missing TFs
      tfCoopNetwork[Idx] <- ppi[indIdx,3];
      #Motif matrix
      regulatoryNetwork=matrix(0,num.TFs,num.genes)
      colnames(regulatoryNetwork)=gene.names
      rownames(regulatoryNetwork)=tf.names
      Idx1=match(motif[,1], tf.names);
      Idx2=match(motif[,2], gene.names);
      Idx=(Idx2-1)*num.TFs+Idx1;
      indIdx=!is.na(Idx)
      Idx=Idx[indIdx] #remove missing genes
      regulatoryNetwork[Idx]=motif[indIdx,3]*epifilter[indIdx,3];          
    }
    num.conditions <- ncol(expr)
    if (randomize=='within.gene'){
      expr <- t(apply(expr, 1, sample))
      if(progress)
        print("Randomizing by reordering each gene's expression")
    } else if (randomize=='by.gene'){
      rownames(expr) <- sample(rownames(expr))
      expr           <- expr[order(rownames(expr)),]
      if(progress)
        print("Randomizing by reordering each gene labels")
    }
  }
  
  if (mode=='legacy'){
    # Create vectors for TF names and Gene names from motif dataset
    tf.names   <- sort(unique(motif[,1]))
    gene.names <- sort(unique(rownames(expr)))
    num.TFs    <- length(tf.names)
    num.genes  <- length(gene.names)
  }
  
  # Bad data checking
  if (num.genes==0){
    stop("Error validating data.  No matched genes.\n  Please ensure that gene names in expression data match gene names in motif data")
  }
  
  if(num.conditions==0) {
    warning('No expression data given.  SPIDER will run based on an identity co-regulation matrix')
    geneCoreg <- diag(num.genes)
  } else if(num.conditions<3) {
    warning('Not enough expression conditions detected to calculate correlation. Co-regulation network will be initialized to an identity matrix.')
    geneCoreg <- diag(num.genes)
  } else {
    
    if(scale.by.present){
      num.positive=(expr>0)%*%t((expr>0))
      geneCoreg <- cor(t(expr), method=cor.method, use="pairwise.complete.obs")*(num.positive/num.conditions)
    } else {
      geneCoreg <- cor(t(expr), method=cor.method, use="pairwise.complete.obs")
    }
    if(progress)
      print('Verified sufficient samples')
  }
  if (any(is.na(geneCoreg))){ #check for NA and replace them by zero
    diag(geneCoreg)=1
    geneCoreg[is.na(geneCoreg)]=0
  }
  
  if (any(duplicated(motif))) {
    warning("Duplicate edges have been found in the motif data. Weights will be summed.")
    motif <- aggregate(motif[,3], by=list(motif[,1], motif[,2]), FUN=sum)
  }
  
  # Prior Regulatory Network
  if(mode=='legacy'){
    Idx1=match(motif[,1], tf.names);
    Idx2=match(motif[,2], gene.names);
    Idx=(Idx2-1)*num.TFs+Idx1;
    regulatoryNetwork=matrix(data=0, num.TFs, num.genes);
    regulatoryNetwork[Idx]=motif[,3]
    colnames(regulatoryNetwork) <- gene.names
    rownames(regulatoryNetwork) <- tf.names
    # PPI data
    # If no ppi data is given, we use the identity matrix
    tfCoopNetwork <- diag(num.TFs)
    # Else we convert our two-column data.frame to a matrix
    if (!is.null(ppi)){
      if(any(duplicated(ppi))){
        warning("Duplicate edges have been found in the PPI data. Weights will be summed.")
        ppi <- aggregate(ppi[,3], by=list(ppi[,1], ppi[,2]), FUN=sum)
      }
      if(remove.missing.ppi){
        # remove edges in the PPI data that target TFs not in the motif
        n <- nrow(ppi)
        ppi <- ppi[which(ppi[,1]%in%tf.names & ppi[,2]%in%tf.names),]
        message(sprintf("%s PPI edges removed that were not present in motif", n-nrow(ppi)))
      }
      Idx1 <- match(ppi[,1], tf.names);
      Idx2 <- match(ppi[,2], tf.names);
      Idx <- (Idx2-1)*num.TFs+Idx1;
      tfCoopNetwork[Idx] <- ppi[,3];
      Idx <- (Idx1-1)*num.TFs+Idx2;
      tfCoopNetwork[Idx] <- ppi[,3];
    }
    colnames(tfCoopNetwork) <- tf.names
    rownames(tfCoopNetwork) <- tf.names
  }
  
  ## Run SPIDER ##
  tic=proc.time()[3]
  
  # adjusting degree distribution
  regulatoryNetwork = degreeAdjust(regulatoryNetwork)
  
  if(progress)
    print('Normalizing networks...')
  regulatoryNetwork = normalizeNetwork(regulatoryNetwork)
  tfCoopNetwork     = normalizeNetwork(tfCoopNetwork)
  geneCoreg         = normalizeNetwork(geneCoreg)
  
  if(progress)
    print('Learning Network...')
  
  minusAlpha = 1-alpha
  step=0
  hamming_cur=1
  if(progress)
    print("Using tanimoto similarity")
  while(hamming_cur>hamming){
    if ((!is.na(iter))&&step>=iter){
      print(paste("Reached maximum iterations, iter =",iter),sep="")
      break
    }
    Responsibility=tanimoto(tfCoopNetwork, regulatoryNetwork)
    Availability=tanimoto(regulatoryNetwork, geneCoreg)
    RA = 0.5*(Responsibility+Availability)
    
    hamming_cur=sum(abs(regulatoryNetwork-RA))/(num.TFs*num.genes)
    regulatoryNetwork=minusAlpha*regulatoryNetwork + alpha*RA
    
    ppi=tanimoto(regulatoryNetwork, t(regulatoryNetwork))
    ppi=update.diagonal(ppi, num.TFs, alpha, step)
    tfCoopNetwork=minusAlpha*tfCoopNetwork + alpha*ppi
    
    CoReg2=tanimoto(t(regulatoryNetwork), regulatoryNetwork)
    CoReg2=update.diagonal(CoReg2, num.genes, alpha, step)
    geneCoreg=minusAlpha*geneCoreg + alpha*CoReg2
    
    if(progress)
      message("Iteration", step,": hamming distance =", round(hamming_cur,5))
    step=step+1
  }
  
  toc=proc.time()[3] - tic
  if(progress)
    message("Successfully ran SPIDER on ", num.genes, " Genes and ", num.TFs, " TFs.\nTime elapsed:", round(toc,2), "seconds.")
  prepResult(zScale, output, regulatoryNetwork, geneCoreg, tfCoopNetwork, edgelist, motif)
}

#' Function to adjust the degree so that the hub nodes are not penalized in z-score transformation
#'
#' @param A Input adjacency matrix
degreeAdjust <- function(A){
  k1 <- colSums(A)/dim(A)[1]
  k2 <- rowSums(A)/dim(A)[2]
  B <- (matrix(replicate(dim(A)[1],k1),nrow=dim(A)[1]))^2
  B <- B + (matrix(t(replicate(dim(A)[2],k2)),nrow=dim(A)[1]))^2
  A <- A * sqrt(B);
}

normalizeNetwork<-function(X){
  X <- as.matrix(X)
  
  nr = nrow(X)
  nc = ncol(X)
  dm = c(nr,nc)
  
  # overall values
  mu0=mean(X)
  std0=sd(X)*sqrt((nr*nc-1)/(nr*nc))
  
  # operations on rows
  mu1=rowMeans(X) # operations on rows
  std1=rowSds(X)*sqrt((nc-1)/nc)
  
  mu1=rep(mu1, nc)
  dim(mu1) = dm
  std1=rep(std1,nc)
  dim(std1)= dm
  
  Z1=(X-mu1)/std1
  
  # operations on columns
  mu2=colMeans(X) # operations on columns
  std2=colSds(X)*sqrt((nr-1)/nr)
  
  mu2 = rep(mu2, each=nr)
  dim(mu2) = dm
  std2= rep(std2, each=nr)
  dim(std2) = dm
  
  Z2=(X-mu2)/std2
  
  # combine and return
  normMat=Z1/sqrt(2)+Z2/sqrt(2)
  
  # checks and defaults for missing data
  Z0=(X-mu0)/std0;
  f1=is.na(Z1); f2=is.na(Z2);
  normMat[f1]=Z2[f1]/sqrt(2)+Z0[f1]/sqrt(2);
  normMat[f2]=Z1[f2]/sqrt(2)+Z0[f2]/sqrt(2);
  normMat[f1 & f2]=2*Z0[f1 & f2]/sqrt(2);
  
  normMat
}

tanimoto<-function(X,Y){
  
  nc = ncol(Y)
  nr = nrow(X)
  dm = c(nr,nc)
  
  Amat=(X %*% Y)
  Bmat=colSums(Y*Y)
  
  Bmat = rep(Bmat,each=nr)
  dim(Bmat) = dm
  #Bmat=matrix(rep(Bmat, each=nr), dm)
  
  Cmat=rowSums(X*X)
  Cmat=rep(Cmat,nc)
  dim(Cmat) = dm
  #Cmat=matrix(rep(Cmat, nc), dm)
  
  den = (Bmat+Cmat-abs(Amat))
  Amat=Amat/sqrt(den)
  
  return(Amat)
}

update.diagonal<-function(diagMat, num, alpha, step){
  seqs = seq(1, num*num, num+1)
  diagMat[seqs]=NaN;
  diagstd=rowSds(diagMat,na.rm=TRUE)
  diagstd[is.na(diagstd)]=0#replace NA with zeros
  diagstd=diagstd*sqrt( (num-2)/(num-1) );
  diagMat[seqs]=diagstd*num*exp(2*alpha*step);
  return(diagMat);
}

prepResult <- function(zScale, output, regulatoryNetwork, geneCoreg, tfCoopNetwork, edgelist, motif){
  resList <- list()
  numGenes = dim(geneCoreg)[1]
  numTFs   = dim(tfCoopNetwork)[1]
  numEdges = sum(regulatoryNetwork!=0)
  if (!zScale){
    regulatoryNetwork <- pnorm(regulatoryNetwork)
    geneCoreg         <- pnorm(geneCoreg)
    tfCoopNetwork     <- pnorm(tfCoopNetwork)
  }
  if("regulatory"%in%output){
    if(edgelist){
      regulatoryNetwork <- melt.array(regulatoryNetwork)
      colnames(regulatoryNetwork) <- c("TF", "Gene", "Score")
      regulatoryNetwork$Motif <- as.numeric(with(regulatoryNetwork, paste0(TF, Gene))%in%paste0(motif[,1],motif[,2]))
    }
    resList$regNet <- regulatoryNetwork
  }
  if("coexpression"%in%output){
    if(edgelist){
      geneCoreg <- melt.array(geneCoreg)
      colnames(geneCoreg) <- c("Gene.x", "Gene.y", "Score")
    }
    resList$coregNet <- geneCoreg
  }
  if("cooperative"%in%output){
    if(edgelist){
      tfCoopNetwork <- melt.array(tfCoopNetwork)
      colnames(tfCoopNetwork) <- c("TF.x", "TF.y", "Score")
    }
    resList$coopNet <- tfCoopNetwork
  }
  pandaObj(regNet=regulatoryNetwork, coregNet=geneCoreg, coopNet=tfCoopNetwork, numGenes=numGenes, numTFs=numTFs, numEdges=numEdges)
}

pandaObj <- setClass("panda", slots=c("regNet","coregNet","coopNet","numGenes","numTFs","numEdges"))
setMethod("show","panda",function(object){print.panda(object)})