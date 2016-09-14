#' MOdeling Network State Transitions from Expression and Regulatory data (MONSTER)
#'
#' This function runs the MONSTER algorithm
#'
#' @param expr Gene Expression dataset
#' @param design Binary vector indicating case control partition
#' @param motif Regulatory data.frame
#' @param nullPerms number of random permutations to run (default 100).  Set to 0 to only 
#' calculate observed transition matrix
#' @param numMaxCores requires doParallel, foreach.  Runs MONSTER in parallel computing 
#' environment.  Set to 1 to avoid parallelization.
#' @param outputDir save MONSTER results in a directory
#' @keywords keywords
#' @export
#' @return An object of class "monster" containing results
#' @examples
#' data(yeast)
#' design <- c(rep(0,ncol(yeast$exp.cc)), rep(1,ncol(yeast$exp.sr)))
#' monsterRes <- monster(cbind(yeast$exp.cc, yeast$exp.sr), 
#'     design, yeast$motif, nullPerms=10, numMaxCores=4)
#' plot(monsterRes)
#' 
#' 
monster <- function(expr, 
                    design, 
                    motif=NULL, 
                    nullPerms=100, 
                    numMaxCores=1, 
                    outputDir=NA, 
                    ...){
    require(doParallel)
    require(foreach)
    
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
        monsterObj(
            tm=transMatrices[[1]], 
            nullTM=transMatrices[-1], 
            numGenes=nrow(expr), 
            numSamples=c(sum(design==0), sum(design==1))))
    
}

#' Checks that data is something MONSTER can handle
#'
#' This function runs the MONSTER algorithm
#'
#' @param expr Gene Expression dataset
#' @return expr Gene Expression dataset in the proper form (may be the same as input)
#' @importFrom assertthat assert_that
#' 
checkDataType <- function(expr){
    requireNamespace("assertthat")
    require(assertthat)
    assert_that(is.data.frame(expr)||is.matrix(expr)||class(expr)=="ExpressionSet")
    if(class(expr)=="ExpressionSet"){
        expr <- exprs(expr)
    }
    if(is.data.frame(expr)){
        expr <- as.matrix(expr)
    }
    expr
}