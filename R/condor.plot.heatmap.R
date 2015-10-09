#' Plot weighted adjacency matrix with links grouped by community
#' 
#' This function will generate the network link 'heatmap' for a weighted network
#' @param condor.object output of either \code{\link{condor.cluster}} or 
#' \code{\link{condor.modularity.max}}
#' @param main plot title
#' @param xlab x axis label
#' @param ylab y axis label
#' @return produces a \code{\link[graphics]{plot}} output.
#' @examples
#' data(small1976)
#' condor.object <- create.condor.object(small1976)
#' condor.object <- condor.cluster(condor.object, project=FALSE)
#' condor.plot.heatmap(condor.object)
#' @import gplots
#' @export
#'  
condor.plot.heatmap = function(condor.object, main="", xlab="blues", ylab="reds"){
  bo <- condor.object
  attach(bo)
  # convert edge lists to adjacency matrices (n reds x m blues)
  adj = get.adjacency(G, attr="weight", sparse=FALSE)
  # reorder reds according to community membership
  reds = as.character(red.memb[order(red.memb[,2]),1])
  adj = adj[reds,]
  # reorder blues according to community membership
  blues = as.character(blue.memb[order(blue.memb[,2]),1])
  adj = adj[,blues]
  rowsep = cumsum(as.vector(table(red.memb[,2])))
  colsep = cumsum(as.vector(table(blue.memb[,2])))
  labCol <- as.character(sort(blue.memb[,2]))
  labCol[duplicated(labCol)] <- ""
  labRow <- as.character(sort(red.memb[,2]))
  labRow[duplicated(labRow)] <- ""
  heatmap.2(adj, Rowv=FALSE, Colv=FALSE, dendrogram="none", keysize=1.25,
            col=colorpanel(10, "white", "black"), scale="none",
            key=TRUE, symkey=FALSE, density.info="none", trace="none",
            main=main, sepcol="#DDDDDD", colsep=colsep, rowsep=rowsep,
            sepwidth = c(0.025, 0.025), ylab=ylab, xlab=xlab, margins=c(3,3),
            labCol=labCol, labRow=labRow, offsetRow=0, offsetCol=0,
            breaks=sort(c(0.1,seq(0, max(adj),length.out=10))))
  detach(bo)
}
