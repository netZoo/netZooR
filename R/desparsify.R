#' De-sparsify gene-level mutation scores into gene set-level mutation scores.
#' @param edgx A binary matrix containing information on which genes belong to which gene sets. Output from the convertgmt function.
#' @param mutratecorx Gene-level mutation scores corrected for the number of gene sets each gene belongs to (from sambar function).
#' @return De-sparsified mutation data.
#' @export
#
desparsify <- function(edgx, mutratecorx, ...){ # edgx=edg, mutratecorx=mutratecor

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
