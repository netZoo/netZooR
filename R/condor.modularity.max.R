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
#' See \code{\link{condor.cluster}} for more details that.
#' @param condor.object is a list created by 
#' \code{\link{create.condor.object}}. \code{condor.object$edges} must 
#' contain the edges in the giant connected component of a bipartite network 
#' @param T0 is a two column data.frame with the initial community 
#' assignment for each "blue" node, assuming there are more reds than blues, 
#' though this is not strictly necessary. The first column contains the 
#' node name, the second column the community assignment.
#' @param weights edgeweights for each edge in \code{edgelist}.
#' @param deltaQmin convergence parameter determining the minimum required increase
#' in the modularity for each iteration. Default is min(10^-4,1/(number of edges)),
#' with number of edges determined by \code{nrow(condor.object$edges)}. User can
#' set this parameter by passing a numeric value to deltaQmin.
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
#' T0 <- data.frame(nodes=blues,coms=1)
#' condor.object <- condor.modularity.max(condor.object,T0=T0)
#' @import Matrix
#' @import nnet
#' @export
#' 
condor.modularity.max = function(condor.object,T0=cbind(1:q,rep(1,q)),weights=1,deltaQmin="default"){

    #assign convergence parameter
    if(deltaQmin == "default"){
        #number of edges
        m = nrow(condor.object$edges)
        deltaQmin <- min(10^-4,1/m)
    }
    if(deltaQmin != "default" & !is.numeric(deltaQmin)){
        stop("deltaQmin must be either 'default' or a numeric value")
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
    
    #Sparse matrix with the upper right block of the true Adjacency matrix. notices the dimension is reds x blues
    A = sparseMatrix(i=reds,j=blues,x=weights,dims=c(length(unique(reds)),length(unique(blues))),index1=TRUE);
    rownames(A) <- red.names
    colnames(A) <- blue.names
    
    p = nrow(A)
    q = ncol(A)
    N = p+q
    
    #Sparse matrix with the upper right block of the true Adjacency matrix. Notice the dimension is reds x blues
    ki = rowSums(A)
    dj = colSums(A)
    m = sum(ki) # m = sum(dj) too
    
    #initialize community assignments for red and blue nodes.
    T0[,1] = as.integer(factor(T0[,1]))
    Tmat <- T0
    R = cbind(1:p,rep(0,length=p))
    cs = sort(unique(Tmat[,2]))
    #variables to track modularity changes after each loop
    Qhist <- vector();
    Qnow  <- 0
    deltaQ <- 1
    #______________________________________
    iter=1
    while(deltaQ > deltaQmin){
        btr <- BTR <- bt <- BT <- vector();
        #calculate T tilde
        for(i in 1:p){
            if(i %% 2500 == 0){print(sprintf("%s%% through iteration %s",round(i/p*100,digits=1),iter))}
            #find the optimal community for node i
            bt <- rep(0,length(cs))
            for(k in cs){
                ind <- Tmat[,2] == k
                if(any(ind)){
                    #bt[k] = sum(A[i,ind] - (ki[i]*dj[ind])/m)
                    bt[k] = sum((A[i,] - (ki[i]*dj)/m)[ind])
                }
            }
            #note that which.max returns the FIRST max if more than one max
            h = which.max(bt)
            if(bt[h] < 0){
                print("making new comm")
                R[i,2] <- max(Tmat[,2])+1
                cs <- sort(c(cs,R[i,2]))
            }
            #if(length(h) > 1){h <- sample(h,1)}
            if(bt[h] >= 0){
                R[i,2] <- h # assign blue vertex i to comm k such that Q is maximized
                bt[-h] <- 0 # BTR is zero if i is not in k (see definition of Q)
            }
            
            #BT <- rbind(BT,bt)
        }
        #calculate R tilde, i.e., B_transpose * R
        for(j in 1:q)
        {
            #initialize jth row of (B_transpose) * R
            btr = rep(0,length(cs))
            #calculate Q for assigning j to different communities
            for(k in cs)
            {
                #if node j is in community k, else BTR[j,k] = 0
                ind <- R[,2] == k
                if(any(ind)){
                    #btr[k] = sum(A[ind,j]-(ki[ind]*dj[j])/m)
                    btr[k] = sum((A[,j]-(ki*dj[j])/m)[ind])
                }
            }
            
            g = which.max(btr)
            #if there is no comm. assignment to increase modularity, make
            #a new community with that node only
            if(btr[g] < 0)
            {
                print("making new comm")
                Tmat[j,2] <- max(R[,2])+1
                btr <- rep(0,length(cs)+1)
                cs <- c(cs,Tmat[j,2])
                btr[length(cs)] <- sum(A[,j]-(ki*dj[j])/m)
                print(btr)
            }
            if(btr[g] >= 0)
            {
                Tmat[j,2] <- g
                btr[-g] <- 0
            }
            
            #add another column to BTR if a new community is added
            if( !is.vector(BTR) && dim(BTR)[2] < length(btr)){BTR <- cbind(BTR,0) }
            BTR <- rbind(BTR,btr)
        }
        
        Tt =  t(sparseMatrix(i=Tmat[,1],j=Tmat[,2],x=1,dims=c(q,length(cs)),index1=TRUE))
        Qthen <- Qnow
        Qcom <- diag(Tt %*% BTR)/m
        Qnow <- sum(Qcom)
        Qhist = c(Qhist,Qnow)
        
        print(paste("Q =",Qnow,sep=" "))
        if(Qnow != 0){
            deltaQ = Qnow - Qthen
        }
        iter=iter+1
    }
    #__________end_while_____________
    qcom_temp <- cbind(Qcom,sort(unique(cs)))
    #drop empty communities
    qcom_out <- qcom_temp[qcom_temp[,1] > 0,]
    #if communities were dropped, relabel so community labels can function
    #as row/column indices in the modularity matrix, B_ij.
    if(nrow(qcom_out) < nrow(qcom_temp)){
        qcom_out[,2] <- as.integer(factor(qcom_out[,2]))
        R[,2] <- as.integer(factor(R[,2]))
        Tmat[,2] <- as.integer(factor(Tmat[,2]))
    }
    colnames(qcom_out) <- c("Qcom","community")
    condor.object$Qcoms = qcom_out
    condor.object$modularity=Qhist
    condor.object$red.memb=data.frame(red.names=red.names[R[,1]],com=R[,2])
    condor.object$blue.memb=data.frame(blue.names=blue.names[Tmat[,1]],com=Tmat[,2])
    
    return(condor.object)
    }
