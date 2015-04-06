#' Summary.panda
#'
#' summarizes the results of a PANDA analysis
#'
#' @param object an object of class "panda"
#' @param ... further arguments passed to or from other methods.
#' @keywords keywords
#' @export
#' @return Summary description of panda S4 object
#' @examples
#' \donttest{
#' data(pandaToyData)
#' panda.res <- panda(pandaToyData$motif,
#'            pandaToyData$expression,pandaToyData$ppi,hamming=.001,progress=TRUE)
#' summary(panda.res)}
summary.panda <- function(object, ...){
    l <- list(coregNet=dim(object@coregNet),regNet=dim(object@regNet),coopNet=dim(object@coopNet))
    message("PANDA network for", nrow(object@coregNet),"genes and",nrow(object@coopNet)," transcription factors.")
}
#' print.panda
#'
#' summarizes the results of a PANDA analysis
#'
#' @param x an object of class "panda"
#' @param ... further arguments passed to or from other methods.
#' @keywords keywords
#' @export
#' @return Summary description of panda S4 object
#' @examples
#' \donttest{
#' data(pandaToyData)
#' panda.res <- panda(pandaToyData$motif,
#'            pandaToyData$expression,pandaToyData$ppi,hamming=.001,progress=TRUE)
#' print(panda.res)}
print.panda <- function(x, ...){
    l <- list(coregNet=dim(x@coregNet),regNet=dim(x@regNet),coopNet=dim(x@coopNet))
    message("PANDA network for", nrow(x@coregNet),"genes and",nrow(x@coopNet),"transcription factors.")
    message("\nSlots:")
    message(slotNames(x)[1],"\t: Regulatory network of",nrow(x@coopNet)," transcription factors to", nrow(x@coregNet),"genes.")
    message(slotNames(x)[2],": Co-regulation network of", nrow(x@coregNet),"genes.")
    message(slotNames(x)[3],"\t: Cooperative network of", nrow(x@coopNet),"transcription factors.\n")
    numEdges <- sum(x@regNet!=0)
    message("Regulatory graph contains ",numEdges,"edges.")
    if (numEdges==nrow(x@regNet)*ncol(x@regNet)){
        message("Regulatory graph is complete.")
    } else {
        message("Regulatory graph is not complete.")
    }
}
#' Plot.panda
#'
#' summarizes the results of a PANDA analysis
#'
#' @param x an object of class "panda"
#' @param ... further arguments passed to or from other methods.
#' @keywords keywords
#' @export
#' @return Plot of the distribution of edge weights in the regulatory network.
#' @examples
#' \donttest{
#' data(pandaToyData)
#' panda.res <- panda(pandaToyData$motif,
#'            pandaToyData$expression,pandaToyData$ppi,hamming=.001,progress=TRUE)
#' plot(panda.res)
#' }
plot.panda <- function(x, ...){
    message("PANDA network for", nrow(x@coregNet),"genes and",nrow(x@coopNet)," transcription factors.")
    message("Mean edge weight = ", mean(x@regNet))
    message("Min edge weight = ", min(x@regNet))
    message("Max edge weight = ", max(x@regNet))
    hist(x@regNet, main="Distribution of edge weights")
}
