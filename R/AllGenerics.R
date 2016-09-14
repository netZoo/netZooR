#' plot.monster
#'
#' plots the sum of squares of off diagonal mass (differential TF Involvement)
#'
#' @param x an object of class "monster"
#' @param ... further arguments passed to or from other methods.
#' @keywords keywords
#' @export
#' @return Plot of the dTFI for each TF against null distribution
#' @examples
#' \donttest{
#' data(yeast)
#' monsterRes <- monster(yeast$exp.ko,c(rep(1,42),rep(0,49),rep(NA,15)),yeast$motif, nullPerms=10, numMaxCores=4)
#' plot(monsterRes)
#' }
plot.monster <- function(x, ...){
    dTFIPlot(x,...)
}
#' print.panda
#'
#' summarizes the results of a MONSTER analysis
#'
#' @param x an object of class "monster"
#' @param ... further arguments passed to or from other methods.
#' @keywords keywords
#' @export
#' @return Description of transition matrices in object
#' @examples
#' \donttest{
#' data(yeast)
#' monster(yeast$exp.ko,c(rep(1,42),rep(0,49),rep(NA,15)),yeast$motif, nullPerms=10, numMaxCores=4)
#' }
print.monster <- function(x, ...){
    cat("MONSTER object\n")
    cat(paste(x@numGenes, "genes\n"))
    cat(paste(x@numSamples[1],"baseline samples\n"))
    cat(paste(x@numSamples[2],"final samples\n"))
    cat(paste("Transition driven by", ncol(x@tm), "transcription factors\n"))
    cat(paste("Run with", length(x@nullTM), "randomized permutations.\n"))
}
