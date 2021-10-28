#' Run EGRET in R
#' 
#' Description:
#'               NOTE: Beta version. EGRET infers individual-specific gene regulatory networks using inidividual level data - 
#'               a genotype vcf file (v) and QBiC binding predictions (q) -  as well as population/reference level data - 
#'                eQTLs (b), a motif-gene prior (m), PPI network (p), and gene expression (e). An annotation file g is also used to
#'                map TF names to their corresponding ensemble ids.
#'
#' Inputs:
#' @param v     : Data frame of VCF file containing SNPs of the individual in question
#' @param b     : Data frame of eQTL data, each row containing an eQTL which exist within motif regions adjacent to the eGene, with columns
#'  TF, gene, variant position,variant chromosome, eQTL beta value.
#' @param q     : Data frame of QBiC predictions of the effect of eQTL variants on TF binding. Each row represents an 
#' eQTL variant with a predicted negative (disruptive) effect on the binding of the TF corresponding to the motif in which the eQTL variant 
#' resides. Colums are: eQTL variant as chr[chrNum]_position, TF, adjacent eGene, QBiC binding effect size and QBiC binding effect (should be negative)  
#' @param m : Motif prior data frame. Each row represents an edge in the bipartite motif prior, with columns TF, gene and edge weight. 
#' The edge weight should be 1 or 0 based on the presence/absence of the TF motif in the promoter region of the gene.
#' @param p : PPI network data frame. Each row represents an edgem with columns TF, TF and interaction weight.
#' @param e  : Gene expression data frame in which each row represents a gene and each column represents the expression of that gene in a sample.
#'  The first column should contain gene IDs.
#' @param g   : Data frame mapping gene names to gene ids, with columns containing the gene ID the corresponding gene name.
#' @param t   : A string containing a name for the EGRET run. Output files will be labelled with this tag.
#'
#' Outputs:
#' @return EGRET    : Predicted genotye-specific gene regulatory network saved as tag_egret.RData
#' @return BASELINE : A Baseline (PANDA) genotype-agnostic gene regulatory network saved as tag_panda.RData
#'
#' @examples
#' 
#' # Run EGRET algorithm
#' toy_qbic_path <- system.file("extdata", "toy_qbic.txt", package = "netZooR", mustWork = TRUE)
#' toy_genotype_path <- system.file("extdata", "toy_genotype.vcf", package = "netZooR", mustWork = TRUE)
#' toy_motif_path <- system.file("extdata", "toy_motif_prior.txt", package = "netZooR", mustWork = TRUE)
#' toy_expr_path <- system.file("extdata", "toy_expr.txt", package = "netZooR", mustWork = TRUE)
#' toy_ppi_path <- system.file("extdata", "toy_ppi_prior.txt", package = "netZooR", mustWork = TRUE)
#' toy_eqtl_path <- system.file("extdata", "toy_eQTL.txt", package = "netZooR", mustWork = TRUE)
#' toy_map_path <- system.file("extdata", "toy_map.txt", package = "netZooR", mustWork = TRUE)
#' qbic <- read.table(file = toy_qbic_path, header = FALSE)
#' vcf <- read.table(toy_genotype_path, header = FALSE, sep = "\t", stringsAsFactors = FALSE, colClasses = c("character", "numeric", "character", "character", "character", "character", "character", "character", "character", "character"))
#' motif <- read.table(toy_motif_path, sep = "\t", header = FALSE)
#' expr <- read.table(toy_expr_path, header = FALSE, sep = "\t", row.names = 1)
#' ppi <- read.table(toy_ppi_path, header = FALSE, sep = "\t")
#' qtl <- read.table(toy_eqtl_path, header = FALSE)
#' nameGeneMap <- read.table(toy_map_path, header = FALSE)
#' tag <- "my_toy_egret_run"
#' \donttest{
#' runEgret(qtl,vcf,qbic,motif,expr,ppi,nameGeneMap,tag)
#' }
#' @export

runEgret <- function(b,v,q,m,e,p,g,t){
  
  # set the inputs
  qbic <- q
  vcf <- v
  motif <- m
  expr <- e
  ppi <- p
  qtl <- b
  tag <- t
  
  # prepare and merge data
  colnames(motif) <- c("tf"  ,    "gene"  ,  "edgeP")
  nameGeneMap <- g
  colnames(nameGeneMap) <- c("gene","name")
  expressedGeneNames <- nameGeneMap$name[which(nameGeneMap$gene %in% rownames(expr))]
  ppiFiltered <- ppi[((ppi$V1 %in% expressedGeneNames) & (ppi$V2 %in% expressedGeneNames)),]
  colnames(qtl) <- c("tf",	"gene", "snpPos",	"chr",	"effect")
  qtl$snpID <- paste0(qtl$chr,"_",qtl$snpPos)
  colnames(qbic) <- c("snpID",	"tf",	"gene",	"qbicEffectSize")
  maxabs <- function(x){
    location <- which(abs(x) == max(abs(x)))
    return(x[location])
  }
  qbic_uniq <- dplyr::distinct(qbic)
  qbic_uniq$catid <- paste0(qbic_uniq$snpID,qbic_uniq$tf,qbic_uniq$gene,qbic_uniq$qbicEffectSize)
  qbic_ag <- aggregate(qbic_uniq$qbicEffectSize, by = list(qbic_uniq$snpID,qbic_uniq$tf,qbic_uniq$gene), FUN = maxabs)
  colnames(qbic_ag) <- c("snpID","tf","gene","qbicEffectSize")
  qbic_ag$catid <- paste0(qbic_ag$snpID,qbic_ag$tf,qbic_ag$gene,qbic_ag$qbicEffectSize)
  qbic_ag$qbicEffect <- qbic_uniq$qbicEffect[match(qbic_ag$catid,qbic_uniq$catid)]
  qbic_ag$absQbicEffectSize <- abs(qbic_ag$qbicEffectSize)
  colnames(vcf) <- c("CHROM",  "POS"  ,   "ID"  ,    "REF"  ,   "ALT"   ,  "QUAL"   , "FILTER" , "INFO"   , "FORMAT", "NA12878")
  snp_ids <- paste0(vcf$CHROM,"_",vcf$POS)
  rownames(vcf) <- snp_ids
  vcf <- tidyr::separate(vcf, NA12878, c("allele1", "allele2"), "\\|", remove = TRUE)
  vcf$alt_allele_count <- as.numeric(vcf$allele1) + as.numeric(vcf$allele2)
  vcf$snp_id <- snp_ids
  qbic_ag$alt_allele_count <- vcf$alt_allele_count[match(qbic_ag$snpID, vcf$snp_id)]
  qtl$alt_allele_count <- vcf$alt_allele_count[match(qtl$snpID, vcf$snp_id)]
  QTL_tf_gene_pairs <- dplyr::distinct(qtl[,c(1:7)])
  QTL_tf_gene_pairs$edgeE <- rep(1,nrow(QTL_tf_gene_pairs))
  QTL_tf_gene_pairs$alt_allele_count[is.na(QTL_tf_gene_pairs$alt_allele_count)] <- 0
  QTL_tf_gene_pairs$qtlTF <- paste0(QTL_tf_gene_pairs$tf,QTL_tf_gene_pairs$snpID)
  
  # make prior
  qbic_ag$snpTF <- paste0(qbic_ag$tf,qbic_ag$snpID)
  QTL_tf_gene_pairs$qbicEffectSize <- qbic_ag$qbicEffectSize[match(QTL_tf_gene_pairs$qtlTF, qbic_ag$snpTF)]
  QTL_tf_gene_pairs$qbicEffect <- qbic_ag$qbicEffect[match(QTL_tf_gene_pairs$qtlTF, qbic_ag$snpTF)]
  QTL_tf_gene_pairs$absQbicEffectSize <- qbic_ag$absQbicEffectSize[match(QTL_tf_gene_pairs$qtlTF, qbic_ag$snpTF)]
  QTL_tf_gene_pairs[is.na(QTL_tf_gene_pairs)] <- 0
  QTL_tf_gene_pairs$absQtlQffect_altAlleleCount_qbicEffect <- abs(QTL_tf_gene_pairs$effect) * QTL_tf_gene_pairs$alt_allele_count * QTL_tf_gene_pairs$qbicEffect
  QTL_tf_gene_pairs$edgeID <- paste0(QTL_tf_gene_pairs$tf,QTL_tf_gene_pairs$gene)
  mods <- unique(QTL_tf_gene_pairs[,c(1,2)])
  mods$edgeID <- paste0(mods$tf,mods$gene)
  absQtlQffect_altAlleleCount_qbicEffect <- aggregate(absQtlQffect_altAlleleCount_qbicEffect ~ edgeID, data=QTL_tf_gene_pairs, FUN=sum)
  mods$absQtlQffect_altAlleleCount_qbicEffect <- absQtlQffect_altAlleleCount_qbicEffect$absQtlQffect_altAlleleCount_qbicEffect[match(mods$edgeID, absQtlQffect_altAlleleCount_qbicEffect$edgeID)]
  
  combined <- merge(as.data.frame(motif),mods, all.x=TRUE, by.x = c(1,2),by.y = c(1,2))
  combined[is.na(combined)] <- 0
  combined$egretPrior <- combined$edgeP - combined$absQtlQffect_altAlleleCount_qbicEffect
  egretPrior <- combined[,c(1,2,6)]
  priorfile <- paste0("priors_",tag,".txt")
  write.table(combined, file = priorfile, col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
  
  # Run message passing
  resultsP <- panda(as.data.frame(motif), expr=expr, ppi=ppiFiltered,progress=TRUE, remove.missing.ppi = TRUE, remove.missing.motif = TRUE, remove.missing.genes = TRUE, mode='legacy')
  filenameP <- paste0(tag,"_panda.RData")
  regnetP <- resultsP@regNet
  save(regnetP, file = filenameP)
  resultsE <- panda(as.data.frame(egretPrior), expr=expr, ppi=ppiFiltered,progress=TRUE, remove.missing.ppi = TRUE, remove.missing.motif = TRUE, remove.missing.genes = TRUE, mode='legacy')
  filenameE <- paste0(tag,"_egret.RData")
  regnetE <- resultsE@regNet
  save(regnetE, file = filenameE)
  
}


