#' Normalize gene mutation scores by gene length.
#' @param x Mutation data, in the format of a matrix, including the number of mutations for samples (rows) and genes (columns).
#' @param cagenes A vector of genes, for example of cancer-associated genes. This will be used to subset the gene-level mutation data to.
#' @param exonsize A vector of gene lengths. This will be used to normalize the gene mutation scores.
#' @return Mutation rate-adjusted gene mutation scores.
#' @export
#
corgenelength <- function(x, cagenes, exonsize, ...){

	# subset mutation data to cancer-associated genes
		x <- x[,which(colnames(x) %in% cagenes)]

	# match with mutation data
		x <- x[,order(colnames(x))]
		x <- x[,which(colnames(x) %in% names(exonsize))]
		exonsize <- exonsize[which(names(exonsize) %in% colnames(x))]

	# normalize the number of mutations by gene length
		y <- sweep(x,2,exonsize,'/')

	return(y) # y is mutrate
}
