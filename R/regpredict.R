#' Bipartite Edge Reconstruction from Expression data
#'
#' This function generates a complete bipartite network from 
#' gene expression data and sequence motif data 
#' 
#' @param motif A motif dataset, a data.frame, matrix or exprSet containing 
#' 3 columns. Each row describes an motif associated with a transcription 
#' factor (column 1) a gene (column 2) and a score (column 3) for the motif.
#' @param expr An expression dataset, as a genes (rows) by samples (columns)
#'  data.frame
#' @param verbose logical to indicate printing of output for algorithm progress.
#' @param method String to indicate algorithm method.  Must be one of 
#' "bere","pearson","panda","cd","lda", or "wcd". Default is "bere"
#' @param score String to indicate whether motif information will be 
#' readded upon completion of the algorithm
#' @param alphaw A weight parameter specifying proportion of weight 
#' to give to indirect compared to direct evidence.  See documentation.
#' @param verbose logical to indicate printing of output for 
#' @param cpp logical use C++ for maximum speed, set to false if unable to run.
#' @keywords keywords
#' @export
#' @return matrix for inferred network between TFs and genes
#' @examples
#' data(yeast)
monsterNI <- function(motif.data, 
                    expr.data,
                    verbose=FALSE,
                    randomize="none",
                    method="bere",
                    alphaw=.5,
                    score="motifincluded",
                    cpp=FALSE){
    if(verbose)
        print('Initializing and validating')
    # Create vectors for TF names and Gene names from Motif dataset
    tf.names   <- sort(unique(motif.data[,1]))
    num.TFs    <- length(tf.names)
    if (is.null(expr.data)){
        stop("Error: Expression data null")
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
        stop("Error validating data.  No matched genes.\n
            Please ensure that gene names in expression 
            file match gene names in motif file.")
    }
  
    strt<-Sys.time()
    if(num.conditions==0) {
        stop("Error: Number of samples = 0")
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

        } else {
            # Standard r correlation calculation
            gene.coreg <- cor(t(expr.data), method="pearson", use="pairwise.complete.obs")
        }
    }
  
    print(Sys.time()-strt)
  
    if(verbose)
        print('More data cleaning')
    # Convert 3 column format to matrix format
    colnames(motif.data) <- c('TF','GENE','value')
    regulatory.network <- tidyr::spread(motif.data, GENE, value, fill=0)
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
    
    
    if(verbose)
        print('Main calculation')
    result <- NULL
    ########################################
    if (method=="BERE"){
        ## Get direct evidence
        directCor <- t(cor(t(exprData),t(exprData[rownames(exprData)%in%tfNames,]))^2)
        
        ## Get the indirect evidence    
        result <- t(apply(tfdcast, 1, function(x){
            cat(".")
            tfTargets <- as.numeric(x)
            
            # Ordinary Logistic Reg
            #         z <- glm(tfTargets ~ ., data=exprData, family="binomial")
            
            # Penalized Logistic Reg
            z <- penalized(tfTargets, exprData, 
                           lambda2=lambda, model="logistic", standardize=TRUE)
            #         z <- optL1(tfTargets, exprData, minlambda1=25, fold=5)
            
            
            predict(z, exprData)
        }))
        
        ## Convert values to ranks
        directCor <- matrix(rank(directCor), ncol=ncol(directCor))
        result <- matrix(rank(result), ncol=ncol(result))
        
        consensus <- directCor*(1-alphaw) + result*alphaw
        rownames(consensus) <- rownames(tfdcast)
        colnames(consensus) <- rownames(exprData)
        consensusRange <- max(consensus)- min(consensus)
        if(score=="motifincluded"){
            consensus <- as.matrix(consensus + consensusRange*tfdcast)
        }
        consensus
    } else if (method=="pearson"){
        result <- t(cor(t(exprData),t(exprData[rownames(exprData)%in%tfNames,]))^2)
        if(score=="motifincluded"){
            result <- as.matrix(consensus + consensusRange*tfdcast)
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

#' Bipartite Edge Reconstruction from Expression data (LDA method)
#'
#' This function generates a complete bipartite network from 
#' gene expression data and sequence motif data 
#' 
#' @param motif A motif dataset, a data.frame, matrix or exprSet 
#' containing 3 columns. Each row describes an motif associated 
#' with a transcription factor (column 1) a gene (column 2) and 
#' a score (column 3) for the motif.
#' @param expr An expression dataset, as a genes (rows) by 
#' samples (columns) data.frame
#' @param verbose logical to indicate printing of output for 
#' algorithm progress.
#' @param method String to indicate algorithm method.  Must 
#' be one of "cd","lda", or "wcd". Default is correlation 
#' difference "cd".
#' @param score String to indicate whether motif information 
#' will be readded upon completion of the algorithm
#' @keywords keywords
#' @export
#' @return TBD, An object of class "bere" (currently matrix or 
#' data.frame) 
#' @examples
#' 1+1
ldaBERE <- function(motifs, expData, score="motifincluded"){
    require(MASS)
    expData <- data.frame(expData)
    tfdcast <- dcast(motifs,V1~V2,fill=0)
    rownames(tfdcast) <- tfdcast[,1]
    tfdcast <- tfdcast[,-1]
    
    expData <- expData[sort(rownames(expData)),]
    tfdcast <- tfdcast[,sort(colnames(tfdcast)),]
    # check that IDs match
    if (prod(rownames(expData)==colnames(tfdcast))!=1){
        stop("ID mismatch")
    }
    result <- t(apply(tfdcast, 1, function(x){
        cat(".")
        tfTargets <- as.numeric(x)
        z <- lda(tfTargets ~ ., expData)
        
        predict(z, expData)$posterior[,2]
    }))
  
    if(score=="motifincluded"){
        result <- as.matrix(result + tfdcast)
    }
    result
}

#' Bipartite Edge Reconstruction from Expression data 
#' (composite method with direct/indirect)
#'
#' This function generates a complete bipartite network from 
#' gene expression data and sequence motif data 
#' 
#' @param motif A motif dataset, a data.frame, matrix or exprSet 
#' containing 3 columns. Each row describes an motif associated 
#' with a transcription factor (column 1) a gene (column 2) 
#' and a score (column 3) for the motif.
#' @param expr An expression dataset, as a genes (rows) by 
#' samples (columns) data.frame
#' @param alpha A weight parameter specifying proportion of weight 
#' to give to indirect compared to direct evidence.  See documentation.
#' @param verbose logical to indicate printing of output for 
#' algorithm progress.
#' @param method String to indicate algorithm method.  Must be 
#' one of "cd","lda", or "wcd". Default is correlation difference 
#' "cd".
#' @param score String to indicate whether motif information will 
#' be readded upon completion of the algorithm
#' @keywords keywords
#' @importFrom reshape2 dcast
#' @export
#' @return TBD, An object of class "bere" (currently matrix or 
#' data.frame) 
#' @examples
#' 1+1
bereFull <- function(motifs, 
                    exprData, 
                    alpha=.5, 
                    penalized=TRUE, 
                    lambda=10, 
                    score="motifincluded"){
    require(MASS)
    exprData <- data.frame(exprData)
    tfdcast <- dcast(motifs,V1~V2,fill=0)
    rownames(tfdcast) <- tfdcast[,1]
    tfdcast <- tfdcast[,-1]
    
    exprData <- exprData[sort(rownames(exprData)),]
    tfdcast <- tfdcast[,sort(colnames(tfdcast)),]
    tfNames <- rownames(tfdcast)[rownames(tfdcast) %in% rownames(exprData)]
    
    ## Filtering
    # filter out the TFs that are not in expression set
    tfdcast <- tfdcast[rownames(tfdcast)%in%tfNames,]
    
    # Filter out genes that aren't targetted by anything 7/28/15
    commonGenes <- intersect(colnames(tfdcast),rownames(exprData))
    exprData <- exprData[commonGenes,]
    tfdcast <- tfdcast[,commonGenes]
    
    # check that IDs match
    if (prod(rownames(exprData)==colnames(tfdcast))!=1){
        stop("ID mismatch")
    }

    ## Get direct evidence
    directCor <- t(cor(t(exprData),t(exprData[rownames(exprData)%in%tfNames,]))^2)
    
    ## Get the indirect evidence    
    result <- t(apply(tfdcast, 1, function(x){
        cat(".")
        tfTargets <- as.numeric(x)
        
        # Ordinary Logistic Reg
#         z <- glm(tfTargets ~ ., data=exprData, family="binomial")
        
        # Penalized Logistic Reg
        z <- penalized(tfTargets, exprData, 
                        lambda2=lambda, model="logistic", standardize=TRUE)
#         z <- optL1(tfTargets, exprData, minlambda1=25, fold=5)
        
        
        predict(z, exprData)
    }))
    
    ## Convert values to ranks
    directCor <- matrix(rank(directCor), ncol=ncol(directCor))
    result <- matrix(rank(result), ncol=ncol(result))
    
    consensus <- directCor*(1-alpha) + result*alpha
    rownames(consensus) <- rownames(tfdcast)
    colnames(consensus) <- rownames(exprData)
    consensusRange <- max(consensus)- min(consensus)
    if(score=="motifincluded"){
        consensus <- as.matrix(consensus + consensusRange*tfdcast)
    }
    consensus
}

#' function for building a network based only on node degree
#'
#' This function creates an unsophisticated graph based solely on 
#' the node degrees
#' 
#' @param motif A motif dataset, a data.frame, matrix or exprSet 
#' containing 3 columns. Each row describes an motif associated 
#' with a transcription factor (column 1) a gene (column 2) and 
#' a score (column 3) for the motif.
#' @keywords keywords
#' @export
#' @return network 
#' @examples
#' 1+1
degreeApproach <- function(motifs){
    tfDegree <- table(motifs[,c(1,3)])
    geneDegree <- table(motifs[,c(2,3)])
    
    tfMatrix <- matrix(rep(tfDegree[,2],nrow(geneDegree)), 
                        ncol=nrow(geneDegree))
    geneMatrix <- t(matrix(rep(geneDegree[,2],nrow(tfDegree)), 
                        nrow=nrow(geneDegree)))
    
    result <- tfMatrix+geneMatrix
    rownames(result) <- rownames(tfDegree)
    colnames(result) <- rownames(geneDegree)
    result
}