#' Convert .gmt files into a binary matrix.
#' @param signature A file containing gene sets (signatures) in .gmt format. These gene sets will be used to de-sparsify the gene-level mutation scores.
#' @param cagenes A vector of genes, for example of cancer-associated genes. This will be used to subset the gene-level mutation data to.
#' @return A matrix containing gene set mutation scores.
#' @export
#
# OBS! cagenes should be optional
# dependencies: utils
convertgmt <- function(signature, cagenes, ...){

	# determine the maximum number of genes a signature can have in the .gmt file
		ncols <- scan(signature, what="character")
		ncols <- length(unique(ncols))

	# read in the signature dataset
		sign <- utils::read.table(signature, header = FALSE, sep = "\t", col.names = paste0("V",seq_len(ncols)), fill = TRUE)
	  row.names(sign) <- sign[,1]
	  sign <- sign[,3:ncol(sign)]
	  sign <- as.matrix(sign)

	# convert the signature dataset to a binary matrix with pathways and genes
		allgenes <- unique(unlist(c(sign)))
		allgenes <- allgenes[order(allgenes)]
		allgenes <- allgenes[which(allgenes!="")] # remove "empty" genes
		signmat <- matrix(0, nrow(sign), length(allgenes))
		row.names(signmat) <- row.names(sign)
		colnames(signmat) <- allgenes
		for(i in 1:nrow(sign)){
			signmat[i,which(allgenes %in% sign[i,])] <- 1
		}

	# subset pathways to cancer-associated genes
		signmat <- signmat[,which(colnames(signmat) %in% cagenes)]

	# remove signatures without mutations
		signmat <- signmat[which(apply(signmat,1,sum)>0),]
		signmat <- signmat[,which(apply(signmat,2,sum)>0)]

	return(signmat)
}
