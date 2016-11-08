#' Passing Messages between Biological Networks to Refine Predicted Interactions
#'
#' This function runs the PANDA algorithm
#'
#' @param motif A motif dataset, a data.frame, matrix or exprSet containing 3 columns.
#' Each row describes an motif associated with a transcription factor (column 1) a
#' gene (column 2) and a score (column 3) for the motif.
#' @param expr An expression dataset, as a genes (rows) by samples (columns) data.frame
#' @param ppi A Protein-Protein interaction dataset, a data.frame containing 3 columns.
#' Each row describes a protein-protein interaction between transcription factor 1(column 1),
#' transcription factor 2 (column 2) and a score (column 3) for the interaction.
#' @param alpha value to be used for update variable, alpha (default=0.1)
#' @param hamming value at which to terminate the process based on hamming distance (default 10^-3)
#' @param iter sets the maximum number of iterations PANDA can run before exiting.
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
#' removed.
#' @param remove.missing.motif Boolean to indicate whether genes targeted in the motif data but not the
#' expression data should be removed.
#' @param remove.missing.genes Boolean to indicate whether genes in the expression data but lacking
#' information from the motif prior should be removed.
#' @param edgelist Boolean to indicate if edge lists instead of matrices should be returned. 
#' @keywords keywords
#' @importFrom matrixStats rowSds
#' @importFrom matrixStats colSds
#' @importFrom Biobase assayData
#' @importFrom reshape melt.array
#' @export
#' @return An object of class "panda" containing matrices describing networks achieved by convergence
#' with PANDA algorithm.\cr
#' "regNet" is the regulatory network\cr
#' "coregNet" is the coregulatory network\cr
#' "coopNet" is the cooperative network
#' @examples
#' data(pandaToyData)
#' pandaRes <- panda(pandaToyData$motif,
#'            pandaToyData$expression,pandaToyData$ppi,hamming=.1,progress=TRUE)
#' @references
#' Glass K, Huttenhower C, Quackenbush J, Yuan GC. Passing Messages Between Biological Networks
#' to Refine Predicted Interactions. PLoS One. 2013 May 318(5):e64832.
panda <- function(motif,expr=NULL,ppi=NULL,alpha=0.1,hamming=0.001,
    iter=NA,output=c('regulatory','coexpression','cooperative'),
    zScale=TRUE,progress=FALSE,randomize=c("None", "within.gene", "by.gene"),cor.method="pearson",
    scale.by.present=FALSE,edgelist=FALSE,remove.missing.ppi=FALSE,
    remove.missing.motif=FALSE,remove.missing.genes=FALSE){

  randomize <- match.arg(randomize)  
  if(progress)
        print('Initializing and validating')

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
          message(sprintf("%s motif edges removed that targeted genes missing in expression data", n-nrow(motif)))
        }
        # Use the motif data AND the expr data (if provided) for the gene list
        # Keep everything sorted alphabetically
        expr <- expr[order(rownames(expr)),]
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
    
    # Create vectors for TF names and Gene names from motif dataset
    tf.names   <- sort(unique(motif[,1]))
    gene.names <- sort(unique(rownames(expr)))
    num.TFs    <- length(tf.names)
    num.genes  <- length(gene.names)

    # Bad data checking
    if (num.genes==0){
        stop("Error validating data.  No matched genes.\n  Please ensure that gene names in expression data match gene names in motif data")
    }

    if(num.conditions==0) {
        warning('No expression data given.  PANDA will run based on an identity co-regulation matrix')
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
    
    if (any(duplicated(motif))) {
      warning("Duplicate edges have been found in the motif data. Weights will be summed.")
      motif <- aggregate(motif[,3], by=list(motif[,1], motif[,2]), FUN=sum)
    }

    # Prior Regulatory Network
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

    ## Run PANDA ##
    tic=proc.time()[3]

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
        message("Successfully ran PANDA on ", num.genes, " Genes and ", num.TFs, " TFs.\nTime elapsed:", round(toc,2), "seconds.")
    prepResult(zScale, output, regulatoryNetwork, geneCoreg, tfCoopNetwork, edgelist, motif)
}

prepResult <- function(zScale, output, regulatoryNetwork, geneCoreg, tfCoopNetwork, edgelist, motif){
    resList <- list()
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
    if("coregulatory"%in%output){
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
    pandaObj(regNet=regulatoryNetwork, coregNet=geneCoreg, coopNet=tfCoopNetwork)
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

dFunction<-function(X,Y){
    A <- cor(X,Y)
    A[A<0]<-0
    A
}

update.diagonal<-function(diagMat, num, alpha, step){
    seqs = seq(1, num*num, num+1)
    diagMat[seqs]=NaN;
    diagstd=rowSds(diagMat,na.rm=TRUE)*sqrt( (num-2)/(num-1) );
    diagMat[seqs]=diagstd*num*exp(2*alpha*step);
    return(diagMat);
}

spreadNet <- function(df){
    df[,3]<- as.numeric(df[,3])
    row_names <- unique(df[,1])
    col_names <- unique(df[,2])
    spread.df <- data.frame(matrix(0,nrow=length(row_names),ncol=length(col_names)),row.names=row_names)
    colnames(spread.df) <- col_names
    for(i in 1:nrow(df)){
        spread.df[as.character(df[i,1]),as.character(df[i,2])] <- df[i,3]
    }
    spread.df
}

#' Top edges
#'
#' topedges gets a network from a panda obj with a specified cutoff based on magnitude of edgeweight.
#'
#' @param x an object of class "panda"
#' @param count an optional integer indicating number of top edges to be included in regulatory network.
#' @param cutoff an optional numeric indicating the z-score edge weight cutoff to be used to identify edges. Default is 2.0.  Not used if count is not NA.
#' @param networks an optional vector specifying which networks to be included in output.  May be any combination of c("coregulation","cooperation","regulatory").
#' @keywords keywords
#' @export
#' @return An object of class "panda" containing binary matrices indicating the existence of an edge between two nodes.  For regulatory network the matrix indicates an edge between a transcription factor (row) and a gene (column)
#' @examples
#' \donttest{
#' data(pandaToyData)
#' pandaRes <- panda(pandaToyData$motif,
#'            pandaToyData$expression,pandaToyData$ppi,hamming=.001,progress=TRUE)
#' topPandaRes <- topedges(pandaRes,1000)
#' }
#' data(pandaResult)
#' topPandaRes <- topedges(pandaResult,1000)
topedges <- function(x, count=NA, cutoff=2.0, networks=c("coregulation","cooperation","regulatory")){
    if(class(x)!="panda"){
        warning(paste(sep="","Cannot run topedges on object of class '",class(x),"'.  Must be of class 'panda'"))
      stop
    }
    if (!is.na(count)){
        cutoff <- sort(x@regNet)[length(x@regNet)-(count-1)]
    }
    regulatoryNetwork <- apply(x@regNet>cutoff, 2,as.numeric)
    rownames(regulatoryNetwork)<-rownames(x@regNet)
    geneCoreg <- apply(x@coregNet>cutoff, 2,as.numeric)
    rownames(geneCoreg)<-rownames(x@coregNet)
    tfCoopNetwork <- apply(x@coopNet>cutoff, 2,as.numeric)
    rownames(tfCoopNetwork)<-rownames(x@coopNet)

    res <- pandaObj(regNet=regulatoryNetwork, coregNet=geneCoreg, coopNet=tfCoopNetwork)
    res
}

#' Subnetwork
#'
#' subnetwork gets a bipartite network containing only the transcription factors or genes and their respective connections
#'
#' @param x an object of class "panda"
#' @param nodes character vector containing the transcription factor or gene labels to subset
#' @param subTf an optional logical indicating whether to subset by transcription factor.  Default is TRUE.
#' @keywords keywords
#' @export
#' @return An matrix describing the subsetted bipartite network.
#' @examples
#' \donttest{
#' data(pandaToyData)
#' pandaRes <- panda(pandaToyData$motif,
#'            pandaToyData$expression,pandaToyData$ppi,hamming=.001,progress=TRUE)
#' topPandaRes <- topedges(pandaRes,1000)
#' subnet.pandaRes <- subnetwork(topPandaRes,c("AR","ARID3A","ELK1"))
#' }
#' data(pandaResult)
#' topPandaRes <- topedges(pandaResult,1000)
#' subnetwork(topPandaRes,c("AR","ARID3A","ELK1"))
subnetwork <- function(x, nodes, subTf=TRUE){
    if(class(x)!="panda"){
        warning(paste(sep="","Cannot run subnetwork on object of class '",class(x),"'.    Must be of class 'panda'"))
        stop
    }
    if (subTf){
        subnet <- x@regNet[nodes,]
        edgeexists <- colSums(subnet)>0
        subnet <- subnet[,edgeexists]
    } else {
        subnet <- x@regNet[,nodes]
        edgeexists <- rowSums(subnet)>0
        subnet <- subnet[edgeexists,]
    }
    subnet
}

#' targetedGenes
#'
#' Gets a set of genes targeted by a specified transcription factor.  This function can be applied to a graph that
#' is not complete, subsetting the edges which have non-zero edge weight.  See function topEdges for dichotomizing edgeweights.
#'
#' @param x an object of class "panda"
#' @param tfs transcription factors to query
#' @keywords keywords
#' @export
#' @return A vector of targeted genes
#' @examples
#' \donttest{
#' data(pandaToyData)
#' pandaRes <- panda(pandaToyData$motif,
#'            pandaToyData$expression,pandaToyData$ppi,hamming=.001)
#' topPandaRes <- topedges(pandaRes,1000)
#' targetedGenes(topPandaRes,c("AR","ELK1"))
#' }
#' data(pandaResult)
#' topPandaRes <- topedges(pandaResult,1000)
targetedGenes <- function(x, tfs){
    if(class(x)!="panda"){
        warning(paste(sep="","Cannot run subnetwork on object of class '",class(x),"'.  Must be of class 'panda'"))
        stop
    }
    subnet <- x@regNet[tfs,,drop=FALSE]
    edgeexists <- colSums(subnet)>0
    targeted <- colnames(x@regNet)[edgeexists]
    targeted
}

#' Plot graph
#'
#' plotGraph plots a bipartite graph
#'
#' @param x an object of class "panda"
#' @keywords keywords
#' @importFrom igraph graph.incidence
#' @importFrom igraph layout.bipartite
#' @export
#' @return An matrix describing the subsetted bipartite network.
#' @examples
#' \donttest{
#' data(pandaToyData)
#' pandaRes <- panda(pandaToyData$motif,
#'            pandaToyData$expression,pandaToyData$ppi,hamming=.001,progress=TRUE)
#' topPandaRes <- topedges(pandaRes,1000)
#' subnet.pandaRes <- subnetwork(topPandaRes,c("AR","ARID3A","ELK1"))
#' plotGraph(subnet.pandaRes)
#' }
#' data(pandaResult)
#' topPandaRes <- topedges(pandaResult, 1000)
#' subnet.pandaRes <- subnetwork(topPandaRes,c("AR","ARID3A","ELK1"))
#' plotGraph(subnet.pandaRes)
plotGraph <- function(x){
    plot(graph.incidence(x), layout=layout.bipartite)
}