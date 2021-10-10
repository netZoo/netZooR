#' Convert .gmt files into a binary matrix.
#' @param signature A file containing gene sets (signatures) in .gmt format. These gene sets will be used to de-sparsify the gene-level mutation scores.
#' @param cagenes A vector of genes, for example of cancer-associated genes. This will be used to subset the gene-level mutation data to.
#' @return A matrix containing gene set mutation scores.
#' @export
#
# OBS! cagenes should be optional
# dependencies: utils
sambar.convertgmt <- function(signature, cagenes){
  
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

#' Normalize gene mutation scores by gene length.
#' @param x Mutation data, in the format of a matrix, including the number of mutations for samples (rows) and genes (columns).
#' @param cagenes A vector of genes, for example of cancer-associated genes. This will be used to subset the gene-level mutation data to.
#' @param exonsize A vector of gene lengths. This will be used to normalize the gene mutation scores.
#' @return Mutation rate-adjusted gene mutation scores.
#' @export
#
sambar.corgenelength <- function(x, cagenes, exonsize){
  
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

#' De-sparsify gene-level mutation scores into gene set-level mutation scores.
#' @param edgx A binary matrix containing information on which genes belong to which gene sets. Output from the sambar.convertgmt function.
#' @param mutratecorx Gene-level mutation scores corrected for the number of gene sets each gene belongs to (from sambar function).
#' @return De-sparsified mutation data.
#' @export
#
sambar.desparsify <- function(edgx, mutratecorx){ # edgx=edg, mutratecorx=mutratecor
  
  # de-sparsify the data
  despar <- matrix(, nrow=nrow(edgx), ncol=ncol(mutratecorx))
  row.names(despar) <- row.names(edgx)
  colnames(despar) <- colnames(mutratecorx)
  for (p in 1:ncol(despar)){
    for (s in 1:nrow(despar)){
      junk <- edgx[s,]
      junk <- names(junk[which(junk==1)]) # check which genes are in signature
      pjunk <- mutratecorx[,p]
      pjunk <- sum(pjunk[which(names(pjunk) %in% junk)]) # sum of mutation scores in signature/pathway s in patient/sample p
      despar[s,p] <- pjunk/length(junk) # correct for pathway length
    }
    # if(p%%10==0){ cat(sprintf("%.0f of %.0f loops\n", p/10, ncol(x)/10)) }
  }
  
  # remove rows and columns with only 0's
  idx1 <- which(apply(despar,1,sum)!=0)
  idx2 <- which(apply(despar,2,sum)!=0)
  if(length(idx1)==0 | length(idx2)==0){
    print("Not enough samples or signatures with mutation scores >0\n")
  } else {
    despar <- despar[idx1,idx2]
  }
  
  return(despar)
}

#' Gene length
#'
#' A vector of gene lengths. This will be used to normalize the gene mutation scores by the gene's length. This example is based on hg19 gene symbols. The gene length is based on the number of non-overlapping exons.
#' Data were downloaded and pre-processed as described in
#' \href{https://doi.org/10.1101/228031}{Kuijjer et al.}
#'
#' @docType data
#' @keywords datasets
#' @name exon.size 
#' @usage data(exon.size)
#' @format A integer vector of size 23459, with gene symbols as names
NULL

#' Example of a gene list
#'
#' List of cancer-associated genes to subset the mutation data to, as described in
#' \href{https://doi.org/10.1101/228031}{Kuijjer et al.}
#'
#' @docType data
#' @keywords datasets
#' @name genes
#' @usage data(genes)
#' @format A character vector of length 2352
NULL

#' Example of mutation data
#'
#' Somatic mutations of Uterine Corpus Endometrial Carcinoma from The Cancer Genome Atlas.
#' Data were downloaded and pre-processed as described in
#' \href{https://doi.org/10.1101/228031}{Kuijjer et al.}
#'
#' @docType data
#' @keywords datasets
#' @name mut.ucec
#' @usage data(mut.ucec)
#' @format A table with 248 rows and 19754 columns
NULL

#' Main SAMBAR function.
#' @param mutdata Mutation data in matrix format. The number of mutations should be listed for samples (rows) and genes (columns).
#' @param esize A integer vector of gene lengths, with gene symbols as names.
#' @param signatureset A file containing gene sets (signatures) in .gmt format. These gene sets will be used to de-sparsify the gene-level mutation scores.
#' @param cangenes A vector of genes, for example of cancer-associated genes. This will be used to subset the gene-level mutation data to.
#' @param kmin The minimum number of subtypes the user wants to assess. Defaults to 2.
#' @param kmax The maximum number of subtypes the user wants to assess. Defaults to 4.
#' @rawNamespace import(vegan, except=diversity)
#' @rawNamespace import(stats, except= c(cov2cor,decompose,toeplitz,lowess,update,spectrum))
#' @return A list of samples and the subtypes to which these samples are assigned, for each k.
#' @export
sambar <- function(mutdata=mut.ucec, esize=exon.size, signatureset=system.file("extdata", "h.all.v6.1.symbols.gmt", package = "netZooR", mustWork = TRUE), cangenes=genes, kmin=2, kmax=4){
  
  # convert gmt file to binary matrix, subset to cancer-associated genes
  edg <- sambar.convertgmt(signature=signatureset, cagenes=cangenes)
  
  # correct number of mutations for gene length (returns gene mutation scores)
  mutlength <- sambar.corgenelength(x=mutdata, cagenes=cangenes, exonsize=esize)
  
  # transform mutlength
  mutlength <- t(mutlength)
  
  # correct for patient-specific mutation rate
  # calculate mutation rate
  patmutrate <- apply(mutlength, 2, sum)
  patmut0 <- which(patmutrate==0)
  # remove patients with mutationrate==0
  if(length(patmut0)>0){
    mutlength <- mutlength[,-patmut0,drop=FALSE]
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
  signpat <- sambar.desparsify(edgx=edg, mutratecorx=mutratecor)
  
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

globalVariables(c("exon.size", "genes", "mut.ucec"))
