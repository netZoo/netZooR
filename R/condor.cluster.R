#' Main clustering function for condor.
#' 
#' This function performs community structure clustering using
#' the bipartite modularity described in
#' \code{\link{condor.modularity.max}}. This function uses a standard
#' (non-bipartite) community structure clustering of the uni-partite,
#' weighted projection of the original bipartite graph as an initial
#' guess for the bipartite modularity.
#' @param condor.object Output of make.condor.object. This function uses 
#' \code{condor.object$edges}
#' @param cs.method is a string to specify which unipartite community 
#' structure algorithm should be used for the seed clustering. 
#' Options are \code{LCS} (\code{\link[igraph]{multilevel.community}}), 
#' \code{LEC} (\code{\link[igraph]{leading.eigenvector.community}}), 
#' \code{FG} (\code{\link[igraph]{fastgreedy.community}}).
#' @param project Provides options for initial seeding of the bipartite 
#' modularity maximization.
#' If TRUE, the nodes in the first column of \code{condor.object$edges}
#' are projected and clustered using \code{cs.method}. If FALSE, the 
#' complete bipartite network is clustered using the unipartite clustering 
#' methods listed in \code{cs.method}.
#' @param low.memory If TRUE, uses \code{\link{condor.modularity.max}}
#' instead of \code{\link{condor.matrix.modularity}}. This is a slower
#' implementation of the modularity maximization, which does not store any
#' matrices in memory. Useful on a machine with low RAM. However, runtimes
#' are (much) longer.
#' @param deltaQmin convergence parameter determining the minimum required increase
#' in the modularity for each iteration. Default is min(10^{-4},1/(number of edges)),
#' with number of edges determined by \code{nrow(condor.object$edges)}. User can
#' set this parameter by passing a numeric value to deltaQmin.
#' @return \code{condor.object} with \code{\link{condor.modularity.max}} output
#' included.
#' @examples 
#' r = c(1,1,1,2,2,2,3,3,3,4,4);
#' b = c(1,2,3,1,2,4,2,3,4,3,4);
#' reds <- c("Alice","Sue","Janine","Mary")
#' blues <- c("Bob","John","Ed","Hank")
#' elist <- data.frame(red=reds[r],blue=blues[b])
#' condor.object <- create.condor.object(elist)
#' condor.object <- condor.cluster(condor.object)
#' @import igraph
#' @import Matrix
#' @export
#' 
condor.cluster <- function(condor.object,cs.method="LCS",project=TRUE,low.memory=FALSE,deltaQmin="default"){
    
    elist <- condor.object$edges
    
    #throw exception if subgraphs are disjoint
    if (max(table(elist[, 1])) == 1) {
      stop("Unable to cluster. Subgraphs are disjoint.")
    }
    
    #extract weights, if any
    if(ncol(elist) > 2){
      weights <- elist[, 3]
    }
    else {
      weights <- rep(1, nrow(elist))
    }
    
    #make sure there's only one connected component
    g.component.test <- graph.data.frame(elist,directed=FALSE)
    if(!is.connected(g.component.test)){
        stop("More than one connected component detected,
             method requires only one connected component")
    }  
    
    G <- graph.data.frame(elist,directed=FALSE)
    
    project.weights <- weights
    #Use unipartite community structure method for first pass
    #project network into gene space to obtain initial community assignment for genes
    if(project){
        #SNP row indices
        reds = as.integer(factor(elist[,1]))
        red.names = levels(factor(elist[,1]))
        #Gene column indices
        blues = as.integer(factor(elist[,2]))
        blue.names = levels(factor(elist[,2]))
        N = max(blues)
        
        #Sparse matrix with the upper right block of the true Adjacency matrix. Notice the dimension is reds x blues
        sM = sparseMatrix(i=reds,j=blues,x=weights,dims=c(length(unique(reds)),length(unique(blues))),index1=T);
        #Project into gene space, projected adjacency matrix has dim = genes x genes
        gM = t(sM) %*% sM;
        rm(sM)
        gc()
        colnames(gM) <- blue.names
        rownames(gM) <- blue.names
        G1 = graph.adjacency(gM,mode="undirected",weighted=TRUE,diag=FALSE);
        #if(clusters(G1)$no > 1){print("Warning more than one component! May cause indexing error")}
        #V(G1)$name <- sort(unique(as.vector(esub[,2])))
        #remove loops and multiple edges
        gcc.initialize = simplify(max.component(G1))
        project.weights <- edge.attributes(gcc.initialize)$weight
    }
    
    #option to treat the bipartite network as if it is unipartite
    #for community initialization only
    if(!project) {
        gcc.initialize <- G 
        blue.names = levels(factor(elist[,2]))
        #blue.indx <- V(G)$name %in% blue.names
    }
    
    if(cs.method=="LCS"){cs0 = multilevel.community(gcc.initialize, weights=project.weights)}
    if(cs.method=="LEC"){cs0 = leading.eigenvector.community(gcc.initialize, weights=project.weights)}
    if(cs.method=="FG"){cs0 = fastgreedy.community(gcc.initialize, weights=project.weights)}
    print(paste("modularity of projected graph",max(cs0$modularity)))
    
    #initial condition for blues' community membership
    if(project){
        T0 <- data.frame(as.integer(factor(blue.names)),as.vector(membership(cs0)))
    }
    else {
        blue.membs <- membership(cs0)[blue.names]
        T0 <- data.frame(as.integer(factor(blue.names)),blue.membs)
    }
    #run bipartite modularity maximization using initial assignments T0
    if(low.memory){
    condor.object <- condor.modularity.max(condor.object,T0=T0,weights=weights,deltaQmin=deltaQmin)
    }
    if(!low.memory){
        condor.object <- condor.matrix.modularity(condor.object,T0=T0,weights=weights,deltaQmin=deltaQmin)
    }
    
    
    return(condor.object)
    }


