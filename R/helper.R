#' Panda Matlab importer
#'
#' Imports the files from the \code{exportPanda.m} file.
#'
#' @param dir Working directory to search for the numeric files.
#' @param celldata Name of the 'celldata.dat' file.
#' @return Two column vector of "regulator" and "target"
#' @importFrom plyr ddply
#' @export
#' @examples
#' 
#' \donttest{
#' # determine gene degree
#'  pandaFiles = importPandaMatlab()
#'  indegree <- ddply(pandaFiles[,2:ncol(pandaFiles)], .(targer), numcolwise(sum))
#'  row.names(indegree) <- indegree[,1]
#'  indegree <- indegree[,-1]
#'  # to export the file
#'  networkfiles = list.files(pattern="numeric")
#'  write.table(indegree,paste("indegree_",networkfiles,sep=""),
#'              sep="\t",quote=F,row.names=T,col.names=T)
#' }
#'
importPandaMatlab<-function(dir=getwd(),celldata="celldata.dat"){
  # read in network output from "exportPANDA.m"
  networkfile <- list.files(dir,pattern="numeric")
  pandaFiles <- read.delim(networkfile, header=F)
  # this is output from "exportPANDA.m"
  edgelist <- read.delim(celldata, sep=" ", header=F)
  pandaFiles <- cbind(edgelist, pandaFiles)
  colnames(pandaFiles)[1:2] <- c("regulator", "target") # regulator and target columns
  pandaFiles
}

#' Network exporter
#'
#' exportNetwork writes the panda network (RegNet) to foreign file formats.
#'
#' @param x An object of class "panda" or an incidence network matrix.
#' @param file A connection or a string giving the file name to output the network.
#' @param format Character string giving the file format supported by igraph.
#' @keywords keywords
#' @importFrom igraph graph.incidence
#' @importFrom igraph write.graph
#' @export
#' @examples
#' \donttest{
#' data(pandaToyData)
#' pandaRes <- panda(pandaToyData$motif, pandaToyData$expression, pandaToyData$ppi)
#' topPandaRes <- topedges(pandaRes, 1000)
#' exportNetwork(topPandaRes, 'topPandaRes.txt')
#' }
exportNetwork <- function(x, file, format='ncol') {
    if (class(x) == "panda") {
        x <- x@regNet
    }
    g <- graph.incidence(x)
    write.graph(g, file, format)
}
