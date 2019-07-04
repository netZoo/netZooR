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
sambar <- function(mutdata=mut.ucec, esize=exon.size, signatureset=system.file("extdata", "h.all.v6.1.symbols.gmt", package = "SAMBAR", mustWork = TRUE), cangenes=genes, kmin=2, kmax=4, ...){

	# convert gmt file to binary matrix, subset to cancer-associated genes
		edg <- convertgmt(signature=signatureset, cagenes=cangenes)

	# correct number of mutations for gene length (returns gene mutation scores)
		mutlength <- corgenelength(x=mutdata, cagenes=cangenes, exonsize=esize)

	# transform mutlength
		mutlength <- t(mutlength)

	# correct for patient-specific mutation rate
		# calculate mutation rate
			patmutrate <- apply(mutlength, 2, sum)
			patmut0 <- which(patmutrate==0)
		# remove patients with mutationrate==0
			if(length(patmut0)>0){
				mutlength <- mutlength[,-patmut0,drop=F]
				patmutrate <- patmutrate[-patmut0]
			}
			
		# correct for mutation rate		
			mutrate <- mutlength
		  for (p in 1:ncol(mutlength)){
		    mutrate[,p] <- mutlength[,p]/patmutrate[p]
		  }

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
