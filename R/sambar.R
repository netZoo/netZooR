#' Main SAMBAR function.
#' @param mutdata Mutation data in matrix format. The number of mutations should be listed for samples (rows) and genes (columns).
#' @param signature A file containing gene sets (signatures) in .gmt format. These gene sets will be used to de-sparsify the gene-level mutation scores.
#' @param cagenes A vector of genes, for example of cancer-associated genes. This will be used to subset the gene-level mutation data to.
#' @param kmin The minimum number of subtypes the user wants to assess. Defaults to 2.
#' @param kmax The maximum number of subtypes the user wants to assess. Defaults to 4.
#' @return A list of samples and the subtypes to which these samples are assigned, for each k.
#' @export
#
# dependencies: vegan, stats
sambar <- function(mutdata=data(mut.ucec), esize=data(exon.size), signatureset=system.file("extdata", "c2.cp.v5.0.symbols.gmt", package = "SAMBAR", mustWork = TRUE), cangenes=data(gene), kmin=2, kmax=4, ...){

	# convert gmt file to binary matrix, subset to cancer-associated genes
		edg <- convertgmt(signature=signatureset, cagenes=genes)

	# correct number of mutations for gene length (returns gene mutation scores)
		mutrate <- corgenelength(x=mut.ucec, cagenes=genes, exonsize=exon.size)

	# transform mutrate ### OBS! probably want to start with patient ids in columns
		mutrate <- t(mutrate)

  # correct gene scores for the number of pathways each gene belongs to
	  mutrate <- mutrate[which(row.names(mutrate) %in% colnames(edg)),]
	  genefreq <- apply(edg,2,sum)
	  genefreq <- genefreq[which(names(genefreq) %in% row.names(mutrate))]
	  mutratecor <- mutrate/genefreq

	# summarize gene mutation scores into pathway mutation scores
		signpat <- desparsify(edgx=edg, mutratecorx=mutratecor)

	# calculate binomial distance between samples
		distance <- vegan::vegdist(t(signpat), method="binomial")

	# cluster tree
		cluster <- stats::hclust(distance, method="complete") # hclust from "stats"

	# cut cluster based on k
		groups <- list()
		cnt <- 1
		for(k in kmin:kmax){
			groups[[cnt]] <- stats::cutree(cluster, k=k)
			cnt <- cnt+1
		}
		names(groups) <- kmin:kmax

	return(groups)
}
