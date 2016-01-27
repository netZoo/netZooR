#' Iteratively maximize bipartite modularity.
#' 
#' This function is based on the bipartite modularity
#' as defined in "Modularity and community detection in bipartite networks"
#' by Michael J. Barber, Phys. Rev. E 76, 066102 (2007)
#' This function uses a slightly different implementation from the paper. It 
#' does not use the "adaptive BRIM" method for identifying the number of 
#' modules. Rather, it simply continues to iterate until the difference in 
#' modularity between iterations is less that 10^-4. Starting from a random 
#' initial condition, this could take some time. Use 
#' \code{\link{condor.cluster}} for quicker runtimes and likely better 
#' clustering, it initializes the blue 
#' node memberships by projecting the blue nodes into a unipartite "blue" 
#' network and then identify communities in that network using a standard 
#' unipartite community detection algorithm run on the projected network.
#' See \code{\link{condor.cluster}} for more details on that.
#' This function loads the entire adjacency matrix in memory, so if your
#' network has more than ~50,000 nodes, you may want to use
#' \code{\link{condor.modularity.max}}, which is slower, but does not store
#' the matrices in memory. Or, of course, you could move to a larger machine.
#' @param condor.object is a list created by 
#' \code{\link{create.condor.object}}. \code{condor.object$edges} must 
#' contain the edges in the giant connected component of a bipartite network 
#' @param T0 is a two column data.frame with the initial community 
#' assignment for each "blue" node, assuming there are more reds than blues, 
#' though this is not strictly necessary. The first column contains the 
#' node name, the second column the community assignment.
#' @param weights edgeweights for each edge in \code{edgelist}.j
#' @return Qcoms data.frame with modularity of each community. 
#' @return modularity modularity value after each iteration.
#' @return red.memb community membership of the red nodes
#' @return blue.memb community membership of the blue.nodes
#' @examples 
#' r = c(1,1,1,2,2,2,3,3,3,4,4);
#' b = c(1,2,3,1,2,4,2,3,4,3,4);
#' reds <- c("Alice","Sue","Janine","Mary")
#' blues <- c("Bob","John","Ed","Hank")
#' elist <- data.frame(red=reds[r],blue=blues[b])
#' condor.object <- create.condor.object(elist)
#' #randomly assign blues to their own community
#' T0 <- data.frame(nodes=blues,coms=1:4)
#' condor.object <- condor.matrix.modularity(condor.object,T0=T0)
#' @import nnet
#' @export
#' 
condor.matrix.modularity = function(condor.object,T0=cbind(1:q,rep(1,q)),weights=1){
    #define a local function
    machine.which.max = function(X){
        thresh <- .Machine$double.eps
        X[abs(X) <= thresh] <- 0
        true_max <- which.max(X)
        return(true_max)
    }
    machine.max = function(X){
        thresh <- .Machine$double.eps
        X[abs(X) <= thresh] <- 0
        true_max <- max(X)
        return(true_max)
    }
    machine.which.is.max = function(X){
        thresh <- .Machine$double.eps
        X[abs(X) <= thresh] <- 0
        true_max <- which.is.max(X)
        return(true_max)
    }
    #Convert the edgelist to a sparseMatrix object
    esub <- condor.object$edges
    #make sure there's only one connected component
    g.component.test <- graph.data.frame(esub,directed=FALSE)
    if(!is.connected(g.component.test)){
        stop("More than one connected component,
             method requires only one connected component")
    }
    reds <- as.integer(factor(esub[,1]))
    red.names <- levels(factor(esub[,1]))
    blues <- as.integer(factor(esub[,2]))
    blue.names <- levels(factor(esub[,2]))
    #ensure that nrows > ncols
    if(length(red.names) < length(blue.names)){
        stop("Adjacency matrix dimension error: This code requires nrows > ncols")
    }
    
    #The upper right block of the true Adjacency matrix. notices the dimension is reds x blues
    #A = sparseMatrix(i=reds,j=blues,x=weights,dims=c(length(unique(reds)),length(unique(blues))),index1=TRUE);
    A = matrix(0,nrow=length(unique(reds)),ncol=length(unique(blues)))
    edges = cbind(reds,blues)
    A[edges] <- weights
    rownames(A) <- red.names
    colnames(A) <- blue.names
    
    p = nrow(A)
    q = ncol(A)
    N = p+q
    
    #initialize community assignments for blue nodes.
    T0[,1] = as.integer(factor(T0[,1]))
    T0 = T0[order(T0[,1]),]
    Tind <- data.matrix(T0)
    Rind = data.matrix(cbind(1:p,rep(0,length=p)))
    cs = sort(unique(Tind[,2]))
    
    #variables to track modularity changes after each loop
    Qhist <- vector();
    Qnow  <- 0
    deltaQ <- 1

    
    #Sparse matrix with the upper right block of the true Adjacency matrix. Notice the dimension is reds x blues
    ki = rowSums(A)
    dj = colSums(A)
    m = sum(ki) # m = sum(dj) too
    
    #Make B tilde, Btilde = Aij - (ki*dj)/m
    Btilde = -(ki %o% dj)/m ###
    Btilde[edges] = weights + Btilde[edges]
    #initialize Tm
    #Tm = sparseMatrix(i=Tind[,1],j=Tind[,2],x=1,index1=TRUE)
    Tm = matrix(0,nrow=q,ncol=max(Tind[,2]))
    Tm[Tind] <- 1 
    #Begin iterations?
    
    #______________________________________
    iter=1
    while(round(deltaQ,digits=4) > 0){
    
    ### Step 1, assign red nodes
    # T tilde = Btilde %*% T
    Ttilde = Btilde %*% Tm
    
    #Find first max for all rows, update R
    Rind[,2] <- unname(apply(Ttilde, 1, machine.which.max))
    
    #Special condition if all nodes are stuck in one large community,
        #if TRUE, randomly assign two nodes to new communities.
    if(iter > 1 && length(unique(c(Rind[,2],Tind[,2]))) == 1){
        random_nodes <- sample(1:nrow(Rind),2)
        Rind[random_nodes,2] <- max(c(Rind[,2],Tind[,2])) + 1:2
    }
    
    #Check to see if new communities should be made
    negative_contribution <- unname(apply(Ttilde, 1,machine.max)) < 0
    if( any( negative_contribution )){
        #add new communities
        num_new_com <- sum(negative_contribution)
        cs <- length(unique(c(Tind[,2],Rind[,2])))
        Rind[negative_contribution,2] <- (cs + 1):(cs + num_new_com)
    }
    
    #Rm = sparseMatrix(i=Rind[,1],j=Rind[,2],x=1,index1=TRUE)
    Rm = matrix(0,nrow=p,ncol=max(Rind[,2]))
    Rm[data.matrix(Rind)] <- 1
    ### Step 2, assign blue nodes
    
    # R tilde = transpose(Btilde) %*% R
    Rtilde = crossprod(Btilde,Rm)
    
    #Find first max for all rows, update T
    Tind[,2] <- unname(apply(Rtilde, 1, machine.which.max))
    
    #Check to see if new communities should be made
    negative_contribution <- unname(apply(Rtilde, 1,machine.max)) < 0
    if( any( negative_contribution )){
        #add new communities
        num_new_com <- sum(negative_contribution)
        cs <- length(unique(c(Tind[,2],Rind[,2])))
        Tind[negative_contribution,2] <- (cs + 1):(cs + num_new_com)
    }
    #Tm dimensions, note the extra empty community.
    #Tm = sparseMatrix(i=Tind[,1],j=Tind[,2],x=1,index1=TRUE) 
    Tm = matrix(0,nrow=q,ncol=max(Tind[,2]))
    Tm[data.matrix(Tind)] <- 1 
    
    Qthen <- Qnow
    #replace with diag(crossprod(T,BTR))/m
    Qcom <- diag(crossprod(Rm,Btilde %*% Tm))/m
    Qnow <- sum(Qcom)
    if(abs(Qnow) < .Machine$double.eps){Qnow <- 0}
    Qhist = c(Qhist,Qnow)
    
    print(paste("Q =",Qnow,sep=" "))
        if(round(Qnow,digits=4) != 0 && round(Qnow,digits=4) != 0){
            deltaQ = Qnow - Qthen
        }
        iter=iter+1
    }
    
    #__________end_while_____________
    
    cs <- 1:ncol(Tm)
    if(any(sort(unique(Tind[,2])) != sort(unique(Rind[,2])))){stop("number of red and blue communities unequal")}
    #drop empty communities
    qcom_temp <- cbind(Qcom,cs)
    qcom_out <- qcom_temp[abs(Qcom) > .Machine$double.eps,]
    
    #if communities were dropped, relabel so community labels can function
    #as row/column indices in the modularity matrix, B_ij.
    if(nrow(qcom_out) < nrow(qcom_temp)){
        qcom_out[,2] <- as.integer(factor(qcom_out[,2]))
        Rind[,2] <- as.integer(factor(Rind[,2]))
        Tind[,2] <- as.integer(factor(Tind[,2]))
    }
    colnames(qcom_out) <- c("Qcom","community")
    condor.object$Qcoms = qcom_out
    condor.object$modularity=Qhist
    condor.object$red.memb=data.frame(red.names=red.names[Rind[,1]],com=Rind[,2])
    condor.object$blue.memb=data.frame(blue.names=blue.names[Tind[,1]],com=Tind[,2])
    
    return(condor.object)
}
