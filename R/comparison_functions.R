#' Calculate regulatory network degree
#'
#' Calculates the transcription factor out-degree or
#' gene in-degree for the estimated panda regulatory network.
#' 
#' @param x An object of class "panda" or matrix
#' @param type Character string - 'tf' or 'gene'
#' @param filter Boolean to force negative degrees to zero
#' @param trim Boolean to trim using topedges or not at a cutoff (weights become binary 1,0)
#' @param ... Options to be passed to topedges function
#' @export
#' @examples
#' \donttest{
#' data(pandaToyData)
#' pandaRes <- panda(pandaToyData$motif,
#'            pandaToyData$expression,pandaToyData$ppi,hamming=.001,progress=TRUE)
#' calcDegree(pandaRes)
#' calcDegree(pandaRes,trim=TRUE,cutoff=1.5)
#' }
#' data(pandaResult)
#' calcDegree(pandaResult,type="tf",trim=TRUE,1000)
#' calcDegree(pandaResult,type="gene",trim=TRUE,1000)
#'
calcDegree <- function(x,type=c("tf","gene"),filter=FALSE,trim=FALSE,...){
	type = match.arg(type)
    if( !(class(x)%in%c("panda","matrix")) ){
      stop(paste(sep="","Cannot run calcDegree on object of class '",class(x),"'.  Must be of class 'panda' or 'matrix'."))
    }
    if(class(x)=="panda"){
    	if(trim==TRUE) x=topedges(x,...)
    	x = x@regNet
    }
    if(type=="tf"){
    	res = rowSums(x)
	} else {
		res = colSums(x)
	}
	if(filter==TRUE){
		if(length(which(res<0))>0) res[res<0] = 0
	}
	return(res)
}

#' Calculate difference in degrees
#'
#' Calculates the transcription factor out-degree or
#' gene in-degree for two different panda regulatory networks.
#' This is useful in comparing networks from two phenotypes.
#' 
#' @param x An object of class "panda" or matrix
#' @param y A second object of class "panda" or matrix
#' @param filter Boolean to force negative degrees to zero
#' @param type Character string - 'tf' or 'gene'
#' @param trim Boolean to trim using topedges or not at a cutoff (weights become binary 1,0)
#' @param ... Options to be passed to topedges function
#' @export
#' @examples
#' \donttest{
#' data(pandaToyData)
#' pandaRes <- panda(pandaToyData$motif,
#'            pandaToyData$expression,pandaToyData$ppi,hamming=.001,progress=TRUE)
#' pandaRes2 <- panda(pandaToyData$motif,
#'            pandaToyData$expression,pandaToyData$ppi,hamming=.1,progress=TRUE)
#' calcDegreeDifference(pandaRes,pandaRes2)
#' calcDegreeDifference(pandaRes,pandaRes2,trim=TRUE,cutoff=1.5)
#' }
#'
calcDegreeDifference <- function(x,y,type=c("tf","gene"),filter=FALSE,trim=FALSE,...){
	xdegree = calcDegree(x,type,trim,...)
	ydegree = calcDegree(y,type,trim,...)
	xdegree-ydegree
}