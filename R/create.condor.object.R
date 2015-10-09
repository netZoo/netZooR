#' Create list amenable to analysis using \code{condor} package.
#' 
#' Converts an edge list into a \code{list} which is then an input for 
#' other functions in the \code{condor} package.
#' @param edgelist a data.frame with 'red' nodes in the first column and
#' 'blue' nodes in the second column, representing links from the node in
#' the first column to the node in the second column. There must be more
#' unique 'red' nodes than 'blue' nodes. Optionally, a third column may be
#' provided to create a weighted network.
#' @param return.gcc if TRUE, returns the giant connected component
#' @return G is an igraph graph object with a 'color' attribute
#' based on the colnames of edgelist. This can be accessed via
#' V(g)$color, which returns a vector indicating red/blue. Use V(g)$name
#' with V(g)$color to identify red/blue node names
#' @return edges corresponding to graph G. If return.gcc=TRUE, includes only
#' those edges in the giant connected component.
#' @return Qcoms output from \code{\link{condor.cluster}} or 
#' \code{\link{condor.modularity.max}}
#' @return modularity \code{NULL} output from \code{\link{condor.cluster}} 
#' or \code{\link{condor.modularity.max}}
#' @return red.memb \code{NULL} output from \code{\link{condor.cluster}} 
#' or \code{\link{condor.modularity.max}}
#' @return blue.memb \code{NULL} output from \code{\link{condor.cluster}} 
#' or \code{\link{condor.modularity.max}}
#' @return qscores \code{NULL} output from \code{\link{condor.qscore}} 
#' @examples 
#' r = c(1,1,1,2,2,2,3,3,3,4,4);
#' b = c(1,2,3,1,2,4,2,3,4,3,4);
#' reds <- c("Alice","Sue","Janine","Mary")
#' blues <- c("Bob","John","Ed","Hank")
#' elist <- data.frame(red=reds[r],blue=blues[b])
#' condor.object <- create.condor.object(elist)
#' 
#' @import igraph
#' @export 
#' 
create.condor.object <- function(edgelist,return.gcc=TRUE){
    
    # make sure first to columns are of class character
    edgelist[, 1] <- as.character(edgelist[, 1])
    edgelist[, 2] <- as.character(edgelist[, 2])
    if(any(is.na(edgelist))) {
        stop("NA's detected. Remove these from edgelist")
    }
    if(any(edgelist[, 1] == "") | any(edgelist[, 2] == "")) {
        stop("Empty strings detected. Remove these from edgelist")
    }
    if(sum(edgelist[, 1] %in% edgelist[, 2]) > 0) {
        stop("edgelist contains one or more nodes that appear in both red and blue columns.
        Check to make sure network is truly bipartite and nodes of each type appear in the
             same column of 'edgelist'.")
    }
    
    g <- graph.data.frame(edgelist,directed=FALSE)
    blue.indx <- V(g)$name %in% unique(edgelist[, 2])
    V(g)$color <- "red"
    V(g)$color[blue.indx] <- "blue"
    
    if(ncol(edgelist) > 2) {
      message("Weights detected. Building condor object with weighted edges.")
      E(g)$weight <- as.numeric(edgelist[, 3])
      # omit 0 edges
      edgelist <- edgelist[edgelist[, 3] != 0, ]
    }
    
    if(!return.gcc){ g.out <- g}
    if(return.gcc){ gcc <- max.component(g); g.out <- gcc }
    blue.names <- V(g.out)$name[V(g.out)$name %in% unique(edgelist[, 2])]
    red.names <- V(g.out)$name[V(g.out)$name %in% unique(edgelist[, 1])]
    edges <- edgelist[edgelist[, 2] %in% blue.names,]
    
    return(list(G=g.out,edges=edges,Qcoms=NULL,modularity=NULL,red.memb=NULL,
                blue.memb=NULL,qscores=NULL))
}


max.component = function(g){
    # return largest connected component of the iGraph graph object g
    g.clust = clusters(g);
    maxclust.id = which(g.clust$csize == max(g.clust$csize))[1];
    h = induced.subgraph(g, which(g.clust$membership == maxclust.id)); # 1-indexed here
    return(h);
}
