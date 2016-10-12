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
#'
#' @param expr Gene Expression dataset, can be matrix or data.frame of expression values or ExpressionSet
#' @param design Binary vector indicating case control partition
#' @param motif Regulatory data.frame consisting of three columns.  For each row, a transcription factor (column 1) 
#' regulates a gene (column 2) with a defined strength (column 3), usually taken to be 0 or 1 
#' @param nullPerms number of random permutations to run (default 100).  Set to 0 to only 
#' calculate observed transition matrix
#' @param numMaxCores requires doParallel, foreach.  Runs MONSTER in parallel computing 
#' environment.  Set to 1 to avoid parallelization.
#' @param outputDir character vector specifying a directory or path in which 
#' which to save MONSTER results, default is NA and results are not saved.
#' @export
#' @import doParallel
#' @import parallel
#' @import foreach
#' @importFrom methods new
#' @return An object of class "monsterAnalysis" containing results
#' @seealso \code{\link{monsterAnalysis-class}}
#' @examples
#' data(yeast)
#' design <- c(rep(0,20),rep(NA,10),rep(1,20))
#' monsterRes <- monster(yeast$exp.cc[1:500,], design, yeast$motif, nullPerms=10, numMaxCores=4)
#' plot(monsterRes)
monster <- function(expr, 
                    design, 
                    motif=NULL, 
                    nullPerms=100, 
                    numMaxCores=1, 
                    outputDir=NA){
    
    # Data type checking
    expr <- checkDataType(expr)

    # Parallelize
    # Initiate cluster
    if(!is.na(numMaxCores)){
        # Calculate the number of cores
        numCores <- detectCores() - 4
        numCores <- min(numCores, numMaxCores)
        
        cl <- makeCluster(numCores)
        registerDoParallel(cl)
        iters <- nullPerms+1 # Two networks for each partition, plus observed partition
        print("Running null permutations in parallel")
        print(paste(numCores,"cores used"))
        print(paste(iters,"network transitions to be estimated"))
    }
    
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
    
    nullExpr <- expr
    transMatrices <- foreach(i=1:iters,
                             .packages=c("MONSTER","reshape2","penalized","MASS")) %dopar% {
        print(paste0("Running iteration ", i))
        if(i!=1){
            nullExpr[] <- expr[sample(1:length(c(expr)))]
        }
        nullExprCases <- nullExpr[,design==1]
        nullExprControls <- nullExpr[,design==0]

        tmpNetCases <- monsterNI(motif, nullExprCases)
        tmpNetControls <- monsterNI(motif, nullExprControls)
        transitionMatrix <- transformation.matrix(
            tmpNetControls, tmpNetCases, remove.diagonal=TRUE, method="ols")    
        print(paste("Finished running iteration", i))
        if (!is.na(outputDir)){
            saveRDS(transitionMatrix,file.path(outputDir,'tms',paste0('tm_',i,'.rds')))
        }
        transitionMatrix
    }
    
    print(Sys.time()-strt)
    if(!is.na(numMaxCores)){
        stopCluster(cl)
    }
    
    gc()
    return(
        monsterAnalysis(
            tm=transMatrices[[1]], 
            nullTM=transMatrices[-1], 
            numGenes=nrow(expr), 
            numSamples=c(sum(design==0), sum(design==1))))
}

#' Checks that data is something MONSTER can handle
#'
#' @param expr Gene Expression dataset
#' @return expr Gene Expression dataset in the proper form (may be the same as input)
#' @importFrom assertthat assert_that
#' @keywords internal
#' @examples
#' expr.matrix <- matrix(rnorm(2000),ncol=20)
#' checkDataType(expr.matrix)
#' #TRUE
#' data(yeast)
#' class(yeast$exp.cc)
#' checkDataType(yeast$exp.cc)
#' #TRUE
#' 
checkDataType <- function(expr){
    assert_that(is.data.frame(expr)||is.matrix(expr)||class(expr)=="ExpressionSet")
    if(class(expr)=="ExpressionSet"){
        expr <- exprs(expr)
    }
    if(is.data.frame(expr)){
        expr <- as.matrix(expr)
    }
    expr
}