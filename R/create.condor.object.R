#' Create list amenable to analysis using \code{CONDOR} package.
#' 
#' Converts an edge list into a \code{list} which is then an input for 
#' other functions in the \code{\link{CONDOR}} package.
#' @param edgelist a two column data.frame with colnames 'red' and 'blue'
#' representing links from the node in the first column to the node in the 
#' second column.
#' @param return.gcc if TRUE, returns the giant connected component
#' @return G is an igraph graph object with a 'color' attribute
#' based on the colnames of edgelist. This can be accessed via
#' V(g)$color, which returns a vector indicating red/blue. Use V(g)$name
#' with V(g)$color to identify red/blue node names
#' @return edges corresponding to graph G. If return.gcc=TRUE, includes only
#' those edges in the giant connected component.
#' @return Qcom output from \code{\link{condor.cluster}} or 
#' \code{\link{BRIM}}
#' @return modularity \code{NULL} output from \code{\link{condor.cluster}} 
#' or \code{\link{BRIM}}
#' @return red.memb \code{NULL} output from \code{\link{condor.cluster}} 
#' or \code{\link{BRIM}}
#' @return blue.memb \code{NULL} output from \code{\link{condor.cluster}} 
#' or \code{\link{BRIM}}
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
    
    if(sum(colnames(edgelist) %in% c("red","blue")) != 2)
    {
        stop("edgelist colnames must be labeled 'red' and 'blue'")
    }
    if(sum(is.na(edgelist$red) + is.na(edgelist$blue)) > 0)
    {
        stop("NA's detected. Remove these from edgelist")
    }
    if(sum(edgelist$red == "NA") + sum(edgelist$blue == "NA") > 0)
    {
        stop("NA's detected. Remove these from edgelist")
    }
    if(sum(edgelist$red == "") + sum(edgelist$blue == "") > 0)
    {
        stop("Empty strings detected. Remove these from edgelist")
    }
    if(sum(edgelist$red %in% edgelist$blue) > 0)
    {
        stop("edgelist contains one or more nodes that appear in both red and blue columns.
        Check to make sure network is truly bipartite and nodes of each type appear in the
             same column of 'edgelist'.")
    }
    
    g <- graph.data.frame(edgelist,directed=FALSE)
    blue.indx <- V(g)$name %in% unique(edgelist$blue)
    V(g)$color <- "red"
    V(g)$color[blue.indx] <- "blue"
    
    if(!return.gcc){ g.out <- g}
    if(return.gcc){ gcc <- max.component(g); g.out <- gcc }
    blue.names <- V(g.out)$name[V(g.out)$name %in% unique(edgelist$blue)]
    red.names <- V(g.out)$name[V(g.out)$name %in% unique(edgelist$red)]
    edges <- edgelist[edgelist$blue %in% blue.names,]
    
    return(list(G=g.out,edges=edges,Qcom=NULL,modularity=NULL,red.memb=NULL,
                blue.memb=NULL,qscores=NULL))
}


max.component = function(g){
    # return largest connected component of the iGraph graph object g
    g.clust = clusters(g);
    maxclust.id = which(g.clust$csize == max(g.clust$csize))[1];
    h = induced.subgraph(g, which(g.clust$membership == maxclust.id)); # 1-indexed here
    return(h);
}
