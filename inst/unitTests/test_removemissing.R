test_removemissing <- function(){
  data(pandaToyData)
  motif <- pandaToyData$motif
  ppi <- pandaToyData$ppi
  
  # add motif edges that target genes missing in the expression data
  motif <- rbind(motif, c("AHR", "MISSING_GENE", 1))
  motif[,3] <- as.numeric(motif[,3])
  # add PPI edges between TFs not present in the motif
  ppi <- rbind(ppi, c("AHR", "MISSING_TF1", 1))
  ppi <- rbind(ppi, c("MISSING_TF2", "AHR", 1))
  ppi <- rbind(ppi, c("MISSING_TF1", "MISSING_TF3", 1))
  ppi[,3] <- as.numeric(ppi[,3])
  panda.r <- panda(motif, pandaToyData$expression, ppi,
                   edgelist=TRUE,
                   remove.missing.ppi=TRUE,
                   remove.missing.motif=TRUE,
                   remove.missing.genes=TRUE)
  reg.r <- panda.r@regNet
  checkTrue(class(reg.r)=="data.frame")
  checkTrue(all(reg.r$TF%in%pandaToyData$motif[,1]))
  checkTrue(all(!reg.r$TF%in%c("MISSING_TF1","MISSING_TF2","MISSING_TF3")))
  checkTrue(all(reg.r$Gene%in%pandaToyData$motif[,2]))
}
