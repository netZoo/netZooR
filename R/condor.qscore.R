#'Calculate Qscore for all nodes 
#'
#'Qscore is designed to calculate the fraction of the modularity 
#'contributed by each node to its community's modularity 
#' @param condor.object output of \code{\link{condor.cluster}} or 
#' \code{\link{condor.modularity.max}}
#' @return condor.object list has \code{condor.object$qscores} added to it.
#' this includes two data.frames, \code{blue.qscore} and \code{red.qscore}
#' which have the qscore for each red and blue node.
#' 
#' @examples
#' r = c(1,1,1,2,2,2,3,3,3,4,4);
#' b = c(1,2,3,1,2,4,2,3,4,3,4);
#' reds <- c("Alice","Sue","Janine","Mary")
#' blues <- c("Bob","John","Ed","Hank")
#' elist <- data.frame(red=reds[r],blue=blues[b])
#' condor.object <- create.condor.object(elist)
#' condor.object <- condor.cluster(condor.object)
#' condor.object <- condor.qscore(condor.object)   
#' 
#' @import igraph
#' @import Matrix
#' @export
#' 
condor.qscore = function(condor.object){
    
        if(is.null(condor.object$red.memb) | is.null(condor.object$blue.memb)){
        stop("Community Memberships missing. Run condor.cluster or condor.modularity.max first!")
        }
    bo <- condor.object
    bo$blue.memb <- bo$blue.memb[order(bo$blue.memb[,"blue.names"]),]
    bo$red.memb <- bo$red.memb[order(bo$red.memb[,"red.names"]),]
    bo$Qcoms <- bo$Qcoms[order(bo$Qcoms[,"community"]),]
    condor.object <- bo
    
    #Convert the edgelist to a sparseMatrix object
    esub <- condor.object$edges
    reds = as.integer(factor(esub[,1]))
    blues = as.integer(factor(esub[,2]))
    
    #extract weights, if any
    if(ncol(esub) > 2){
      weights <- esub[, 3]
    }
    else {
      weights <- rep(1, nrow(esub))
    }
    
    #Sparse matrix with the upper right block of the true Adjacency matrix. notices the dimension is reds x blues
    A = sparseMatrix(i=reds,j=blues,x=weights,dims=c(length(unique(reds)),length(unique(blues))),index1=TRUE);
    if(nrow(A) < ncol(A)){A <- t(A)}
    #if(max(blues) > max(reds)){blues <- reds;}
    p = nrow(A)
    q = ncol(A)
    N = p+q
    ki = rowSums(A)
    dj = colSums(A)
    m = sum(ki) # m = sum(dj) too
    R1 = condor.object$red.memb
    T1 = condor.object$blue.memb
    r1 = cbind(as.numeric(factor(R1[,1])),R1[,2])
    t1 = cbind(as.numeric(factor(T1[,1])),T1[,2]) 
    Rtrans = sparseMatrix(i=r1[,2],j=r1[,1],x=1,dims=c(max(r1[,2]),length(unique(r1[,1]))),index1=TRUE);    
    T2 = sparseMatrix(i=t1[,1],j=t1[,2],x=1,dims=c(max(t1[,1]),max(t1[,2])),index1=TRUE);
    Qcoms <- condor.object$Qcoms
    Qjk = vector(length=q)
    for(j in 1:max(t1[,1])){
        if(j %% 1000 == 0){print(paste(j,t1[j,]))}
        Bj = A[,j] - (ki*dj[j])/m;
        Qjk[j] = ((Rtrans[t1[j,2],] %*% Bj)/(2*m))*(1/Qcoms[t1[j,2],1])
    }  
    Qik = vector(length=p)
    for(i in 1:max(r1[,1])){
        if(i %% 1000 == 0){print(i)}
        Bi = A[i,] - (ki[i]*dj)/m;
        Qik[i] = ((Bi %*% T2[,r1[i,2]])/(2*m))*(1/Qcoms[r1[i,2],1])  
    }    
    condor.object$qscores = list(blue.qscore=data.frame(T1,Qjk),red.qscore=data.frame(R1,Qik))
    return(condor.object)
}