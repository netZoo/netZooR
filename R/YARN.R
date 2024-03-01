#' Annotate your Expression Set with biomaRt
#'
#' @param obj ExpressionSet object.
#' @param genes Genes or rownames of the ExpressionSet.
#' @param filters getBM filter value, see getBM help file.
#' @param attributes getBM attributes value, see getBM help file.
#' @param biomart BioMart database name you want to connect to. Possible database names can be retrieved with teh function listMarts.
#' @param dataset Dataset you want to use. To see the different datasets available within a biomaRt you can e.g. do: mart = useMart('ensembl'), followed by listDatasets(mart).
#' @param ... Values for useMart, see useMart help file.
#'
#' @return ExpressionSet object with a fuller featureData.
#' @export
#'
#' @importFrom biomaRt useMart
#' @importFrom biomaRt getBM
#' @importFrom Biobase featureNames
#' @importFrom Biobase fData
#' @importFrom Biobase fData<-
#' @importFrom Biobase ExpressionSet
#' @importClassesFrom Biobase ExpressionSet
#'
#' @examples
#' download.file('https://netzoo.s3.us-east-2.amazonaws.com/netZooR/unittest_datasets/yarn/bladder.rdata',destfile='netZooR/data/bladder.rdata')
#' download.file('https://netzoo.s3.us-east-2.amazonaws.com/netZooR/unittest_datasets/yarn/skin',destfile='netZooR/data/skin.rdata')
#' data(skin)
#' # subsetting and changing column name just for a silly example
#' skin <- skin[1:10,]
#' colnames(fData(skin)) = paste("names",1:6)
#' biomart<-"ENSEMBL_MART_ENSEMBL";
#' genes <- sapply(strsplit(rownames(skin),split="\\."),function(i)i[1])
#' newskin <-annotateFromBiomart(skin,genes=genes,biomar=biomart)
#' head(fData(newskin)[,7:11])
#'
annotateFromBiomart <- function(obj,genes=featureNames(obj),filters="ensembl_gene_id",
                                attributes=c("ensembl_gene_id","hgnc_symbol","chromosome_name","start_position","end_position"),
                                biomart="ensembl",dataset="hsapiens_gene_ensembl",...){
  mart <- useMart(biomart=biomart,dataset=dataset,...)
  anno <- getBM(attributes = attributes, filters = filters,
                values = genes, mart = mart)
  if(nrow(anno)<length(genes)){
    warning("getBM returned fewer rows than genes queried.")
  }
  if(nrow(anno)>length(genes)){
    warning(sprintf("getBM returned more rows than genes queried. Using first call of %s.",colnames(anno)[1]))
    throw = which(duplicated(anno[,1]))
    anno  = anno[-throw,]
  }
  anno = anno[match(genes,anno[,"ensembl_gene_id"]),]
  if(!is.null(fData(obj))) anno = cbind(fData(obj),anno)
  fData(obj) = anno
  obj
}

#' Check for wrong annotation of a sample using classical MDS and control genes.
#'
#' @param obj ExpressionSet object.
#' @param phenotype phenotype column name in the phenoData slot to check.
#' @param controlGenes Name of controlGenes, ie. 'Y' chromosome. Can specify 'all'.
#' @param columnID Column name where controlGenes is defined in the featureData slot if other than 'all'.
#' @param plotFlag TRUE/FALSE Whether to plot or not
#' @param legendPosition Location for the legend.
#' @param ... Extra parameters for \code{\link{plotCMDS}} function.
#'
#' @importFrom graphics legend
#'
#' @return Plots a classical multi-dimensional scaling of the 'controlGenes'. Optionally returns co-ordinates.
#' @export
#'
#' @examples
#' data(bladder)
#' checkMisAnnotation(bladder,'GENDER',controlGenes='Y',legendPosition='topleft')
#'
checkMisAnnotation <- function(obj, phenotype, controlGenes = "all",
                               columnID = "chromosome_name", plotFlag = TRUE, legendPosition = NULL,
                               ...) {
  if (tolower(controlGenes) != "all") {
    obj <- filterGenes(obj, labels = controlGenes, featureName = columnID,
                       keepOnly = TRUE)
  }
  if (length(phenotype) == 1) {
    phenotype <- factor(pData(obj)[, phenotype])
  }
  res <- plotCMDS(obj, pch = 21, bg = phenotype, plotFlag = plotFlag,
                  ...)
  if (!is.null(legendPosition))
    legend(legendPosition, legend = levels(phenotype), fill = 1:length(levels(phenotype)))
  invisible(res)
}

#' Check tissues to merge based on gene expression profile
#'
#' @param obj ExpressionSet object.
#' @param majorGroups Column name in the phenoData slot that describes the general body region or site of the sample.
#' @param minorGroups Column name in the phenoData slot that describes the specific body region or site of the sample.
#' @param filterFun Filter group specific genes that might disrupt PCoA analysis.
#' @param plotFlag TRUE/FALSE whether to plot or not
#' @param ... Parameters that can go to \code{\link[yarn]{checkMisAnnotation}}
#'
#' @return CMDS Plots of the majorGroupss colored by the minorGroupss. Optional matrix of CMDS loadings for each comparison.
#' @export
#'
#' @seealso checkTissuesToMerge
#'
#' @examples
#' data(skin)
#' checkTissuesToMerge(skin,'SMTS','SMTSD')
#'
checkTissuesToMerge <- function(obj, majorGroups, minorGroups,
                                filterFun = NULL, plotFlag = TRUE, ...) {
  if (length(majorGroups) == 1) {
    region <- factor(pData(obj)[, majorGroups])
  } else {
    region <- factor(majorGroups)
  }
  if (!is.null(filterFun)) {
    obj <- filterFun(obj)
  }
  result <- lapply(levels(region), function(i) {
    keepSamples <- which(region == i)
    objSubset <- obj[, keepSamples]
    objSubset <- objSubset[which(rowSums(exprs(objSubset)) > 0), ]
    res <- checkMisAnnotation(objSubset, phenotype = minorGroups,
                              controlGenes = "all", plotFlag = plotFlag, main = i,
                              ...)
    res
  })
  invisible(result)
}

#' Download GTEx files and turn them into ExpressionSet object
#'
#' Downloads the V6 GTEx release and turns it into an ExpressionSet object.
#'
#' @param type Type of counts to download - default genes.
#' @param file File path and name to automatically save the downloaded GTEx expression set. Saves as a RDS file.
#' @param ... Does nothing currently.
#'
#' @return Organized ExpressionSet set.
#' @export
#'
#' @importFrom downloader download
#' @importFrom readr read_tsv
#' @importFrom readr problems
#' @importFrom Biobase AnnotatedDataFrame
#' @importFrom Biobase phenoData<-
#' @importFrom Biobase pData<-
#'
#' @examples
#' # obj <- downloadGTEx(type='genes',file='~/Desktop/gtex.rds')
downloadGTEx <- function(type = "genes", file = NULL, ...) {
  phenoFile <- "http://www.gtexportal.org/static/datasets/gtex_analysis_v6/annotations/GTEx_Data_V6_Annotations_SampleAttributesDS.txt"
  pheno2File <- "http://www.gtexportal.org/static/datasets/gtex_analysis_v6/annotations/GTEx_Data_V6_Annotations_SubjectPhenotypesDS.txt"
  geneFile <- "http://www.gtexportal.org/static/datasets/gtex_analysis_v6/rna_seq_data/GTEx_Analysis_v6_RNA-seq_RNA-SeQCv1.1.8_gene_reads.gct.gz"
  
  message("Downloading and reading files")
  pdFile <- tempfile("phenodat", fileext = ".txt")
  download(phenoFile, destfile = pdFile)
  pd <- read_tsv(pdFile)
  pd <- as.matrix(pd)
  rownames(pd) <- pd[, "SAMPID"]
  ids <- sapply(strsplit(pd[, "SAMPID"], "-"), function(i) paste(i[1:2],
                                                                 collapse = "-"))
  
  pd2File <- tempfile("phenodat2", fileext = ".txt")
  download(pheno2File, destfile = pd2File)
  pd2 <- read_tsv(pd2File)
  pd2 <- as.matrix(pd2)
  rownames(pd2) <- pd2[, "SUBJID"]
  pd2 <- pd2[which(rownames(pd2) %in% unique(ids)), ]
  pd2 <- pd2[match(ids, rownames(pd2)), ]
  rownames(pd2) <- colnames(counts)
  
  pdfinal <- AnnotatedDataFrame(data.frame(cbind(pd, pd2)))
  
  if (type == "genes") {
    countsFile <- tempfile("counts", fileext = ".gz")
    download(geneFile, destfile = countsFile)
    cnts <- suppressWarnings(read_tsv(geneFile, skip = 2))
    genes <- unlist(cnts[, 1])
    geneNames <- unlist(cnts[, 2])
    counts <- cnts[, -c(1:2)]
    counts <- as.matrix(counts)
    rownames(counts) <- genes
    for (i in 1:nrow(problems(cnts))) {
      counts[problems(cnts)$row[i], problems(cnts)$col[i]] <- 1e+05
    }
    throwAway <- which(rowSums(counts) == 0)
    counts <- counts[-throwAway, ]
    genes <- sub("\\..*", "", rownames(counts))
    
    host <- "dec2013.archive.ensembl.org"
    biomart <- "ENSEMBL_MART_ENSEMBL"
    dataset <- "hsapiens_gene_ensembl"
    attributes <- c("ensembl_gene_id", "hgnc_symbol", "chromosome_name",
                    "start_position", "end_position", "gene_biotype")
  }
  
  message("Creating ExpressionSet")
  pdfinal <- pdfinal[match(colnames(counts), rownames(pdfinal)),
  ]
  es <- ExpressionSet(as.matrix(counts))
  phenoData(es) <- pdfinal
  pData(es)["GTEX-YF7O-2326-101833-SM-5CVN9", "SMTS"] <- "Skin"
  pData(es)["GTEX-YEC3-1426-101806-SM-5PNXX", "SMTS"] <- "Stomach"
  
  message("Annotating from biomaRt")
  es <- annotateFromBiomart(obj = es, genes = genes, host = host,
                            biomart = biomart, dataset = dataset, attributes = attributes)
  
  message("Cleaning up files")
  unlink(pdFile)
  unlink(pd2File)
  unlink(countsFile)
  
  if (!is.null(file))
    saveRDS(es, file = file)
  return(es)
}

#' Extract the appropriate matrix
#'
#' This returns the raw counts, log2-transformed raw counts, or normalized expression.
#' If normalized = TRUE then the log paramater is ignored.
#'
#' @param obj ExpressionSet object or objrix.
#' @param normalized TRUE / FALSE, use the normalized matrix or raw counts
#' @param log TRUE/FALSE log2-transform.
#'
#' @importFrom Biobase assayData
#' @return matrix
#' @examples
#'
#' data(skin)
#' head(netZooR:::extractMatrix(skin,normalized=FALSE,log=TRUE))
#' head(netZooR:::extractMatrix(skin,normalized=FALSE,log=FALSE))
#'
extractMatrix <- function(obj, normalized = FALSE, log = TRUE) {
  if (class(obj) == "ExpressionSet") {
    if (!normalized) {
      obj <- exprs(obj)
    } else {
      if (!"normalizedMatrix" %in% names(assayData(obj)))
        stop("normalizedMatrix assayData missing")
      obj <- assayData(obj)[["normalizedMatrix"]]
      if (log & normalized)
        message("normalizedMatrix is assumed to already be log-transformed")
      log <- FALSE
    }
  }
  if (log == TRUE) {
    obj <- log2(obj + 1)
  }
  obj
}

#' Filter specific genes
#'
#' The main use case for this function is the removal of sex-chromosome genes.
#' Alternatively, filter genes that are not protein-coding.
#'
#' @param obj ExpressionSet object.
#' @param labels Labels of genes to filter or keep, eg. X, Y, and MT
#' @param featureName FeatureData column name, eg. chr
#' @param keepOnly Filter or keep only the genes with those labels
#'
#' @return Filtered ExpressionSet object
#' @export
#'
#' @importFrom Biobase exprs
#' @importFrom Biobase fData
#'
#' @examples
#' data(skin)
#' filterGenes(skin,labels = c('X','Y','MT'),featureName='chromosome_name')
#' filterGenes(skin,labels = 'protein_coding',featureName='gene_biotype',keepOnly=TRUE)
#'
filterGenes <- function(obj, labels = c("X", "Y", "MT"), featureName = "chromosome_name",
                        keepOnly = FALSE) {
  features <- fData(obj)[, featureName]
  if (keepOnly == FALSE) {
    throwAwayGenes <- which(features %in% labels)
  } else {
    throwAwayGenes <- which(!features %in% labels)
  }
  obj <- obj[-throwAwayGenes, ]
  obj
}

#' Filter genes that have less than a minimum threshold CPM for a given group/tissue
#'
#' @param obj ExpressionSet object.
#' @param groups Vector of labels for each sample or a column name of the phenoData slot.
#' for the ids to filter. Default is the column names.
#' @param threshold The minimum threshold for calling presence of a gene in a sample.
#' @param minSamples Minimum number of samples - defaults to half the minimum group size.
#' @param ... Options for \link[edgeR]{cpm}.
#' @seealso \link[edgeR]{cpm} function defined in the edgeR package.
#'
#' @return Filtered ExpressionSet object
#' @export
#'
#' @importFrom edgeR cpm
#' @importFrom Biobase exprs
#' @importFrom Biobase pData
#'
#' @examples
#' data(skin)
#' filterLowGenes(skin,'SMTSD')
#'
filterLowGenes <- function(obj, groups, threshold = 1, minSamples = NULL,
                           ...) {
  if (is.null(minSamples)) {
    if (length(groups) == 1) {
      minSamples <- min(table(pData(obj)[, groups]))/2
    } else {
      minSamples <- min(table(groups))/2
    }
  }
  counts <- cpm(exprs(obj), ...)
  keep <- rowSums(counts > threshold) >= minSamples
  obj <- obj[keep, ]
  obj
}

#' Filter genes not expressed in any sample
#'
#' The main use case for this function is the removal of missing genes.
#'
#' @param obj ExpressionSet object.
#' @param threshold Minimum sum of gene counts across samples -- defaults to zero.
#'
#' @return Filtered ExpressionSet object
#' @export
#'
#' @importFrom Biobase exprs
#' @importFrom Biobase fData
#'
#' @examples
#' data(skin)
#' filterMissingGenes(skin)
#'
filterMissingGenes <- function(obj, threshold = 0) {
  sumGenes <- rowSums(exprs(obj))
  throwAwayGenes <- which(sumGenes <= threshold)
  if (length(which(sumGenes <= 0)) > 0) {
    obj <- obj[-throwAwayGenes, ]
  }
  obj
}

#' Filter samples
#'
#' @param obj ExpressionSet object.
#' @param ids Names found within the groups labels corresponding to samples to be removed
#' @param groups Vector of labels for each sample or a column name of the phenoData slot
#' for the ids to filter. Default is the column names.
#' @param keepOnly Filter or keep only the samples with those labels.
#'
#' @return Filtered ExpressionSet object
#' @export
#'
#' @importFrom Biobase pData
#'
#' @examples
#' data(skin)
#' filterSamples(skin,ids = "Skin - Not Sun Exposed (Suprapubic)",groups="SMTSD")
#' filterSamples(skin,ids=c("GTEX-OHPL-0008-SM-4E3I9","GTEX-145MN-1526-SM-5SI9T"))
#'
filterSamples <- function(obj, ids, groups = colnames(obj), keepOnly = FALSE) {
  if (length(groups) == 1) {
    groups <- pData(obj)[, groups]
  }
  throwAway <- which(groups %in% ids)
  if (keepOnly) {
    obj <- obj[, throwAway]
  } else {
    obj <- obj[, -throwAway]
  }
  obj
}

#' Normalize in a tissue aware context
#'
#' This function provides a wrapper to various normalization methods developed.
#' Currently it only wraps qsmooth and quantile normalization returning a log-transformed
#' normalized matrix. qsmooth is a normalization approach that normalizes samples in
#' a condition aware manner.
#'
#' @param obj ExpressionSet object
#' @param groups Vector of labels for each sample or a column name of the phenoData slot
#' for the ids to filter. Default is the column names
#' @param normalizationMethod Choice of 'qsmooth' or 'quantile'
#' @param ... Options for \code{\link{qsmooth}} function or \code{\link[limma]{normalizeQuantiles}}
#'
#' @return ExpressionSet object with an assayData called normalizedMatrix
#' @export
#'
#' @source The function qsmooth comes from the qsmooth packages
#' currently available on github under user 'kokrah'.
#'
#' @importFrom limma normalizeQuantiles
#' @importFrom Biobase storageMode
#' @importFrom Biobase storageMode<-
#' @importFrom Biobase assayData
#' @importFrom Biobase assayData<-
#' @importFrom preprocessCore normalize.quantiles
#' @importClassesFrom Biobase eSet
#' @importClassesFrom Biobase ExpressionSet
#'
#' @examples
#' data(skin)
#' normalizeTissueAware(skin,"SMTSD")
#'
normalizeTissueAware <- function(obj, groups, normalizationMethod = c("qsmooth",
                                                                      "quantile"), ...) {
  normalizationMethod <- match.arg(normalizationMethod)
  if (length(groups) == 1) {
    groups <- factor(pData(obj)[, groups])
  }
  storageMode(obj) <- "environment"
  if (normalizationMethod == "qsmooth") {
    normalizedMatrix <- qsmooth(obj, groups = groups,
                                ...)
  } else if (normalizationMethod == "quantile") {
    if(length(unique(groups))>1){
      normalizedMatrix <- sapply(unique(groups), function(i) {
        cnts <- exprs(obj[, which(pData(obj)$our %in% i)])
        nmat <- normalize.quantiles(cnts)
        colnames(nmat) <- colnames(cnts)
        nmat
      })
      normalizedMatrix <- Reduce("cbind", normalizedMatrix)
      normalizedMatrix <- normalizedMatrix[, match(colnames(obj),
                                                   colnames(normalizedMatrix))]
    } else {
      normalizedMatrix <- normalize.quantiles(exprs(obj))
      colnames(normalizedMatrix) <- colnames(obj)
    }
  }
  assayData(obj)[["normalizedMatrix"]] <- normalizedMatrix
  storageMode(obj) <- "lockedEnvironment"
  obj
}

#' Plot classical MDS of dataset
#'
#' This function plots the MDS coordinates for the "n" features of interest. Potentially uncovering batch
#' effects or feature relationships.
#'
#' @param obj ExpressionSet object or objrix.
#' @param comp Which components to display.
#' @param normalized TRUE / FALSE, use the normalized matrix or raw counts.
#' @param distFun Distance function, default is dist.
#' @param distMethod The distance measure to be used. This must be one of "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski". Any unambiguous substring can be given.
#' @param n Number of features to make use of in calculating your distances.
#' @param samples Perform on samples or genes.
#' @param log TRUE/FALSE log2-transform raw counts.
#' @param plotFlag TRUE/FALSE whether to plot or not.
#' @param ... Additional plot arguments.
#' @return coordinates
#'
#' @importFrom matrixStats rowSds
#' @importFrom stats dist
#' @importFrom stats cmdscale
#' @importFrom graphics plot
#'
#' @export
#' @examples
#' data(skin)
#' res <- plotCMDS(skin,pch=21,bg=factor(pData(skin)$SMTSD))
#' \donttest{
#' # library(calibrate)
#' # textxy(X=res[,1],Y=res[,2],labs=rownames(res))
#' }
plotCMDS <- function(obj, comp = 1:2, normalized = FALSE, distFun = dist,
                     distMethod = "euclidian", n = NULL, samples = TRUE, log = TRUE,
                     plotFlag = TRUE, ...) {
  if (is.null(n))
    n <- min(nrow(obj), 1000)
  obj <- extractMatrix(obj, normalized, log)
  genesToKeep <- which(rowSums(obj) > 0)
  geneVars <- rowSds(obj[genesToKeep, ])
  geneIndices <- genesToKeep[order(geneVars, decreasing = TRUE)[seq_len(n)]]
  obj <- obj[geneIndices, ]
  
  if (samples == TRUE) {
    obj <- t(obj)
  }
  d <- distFun(obj, method = distMethod)
  ord <- cmdscale(d, k = max(comp))
  xl <- paste("MDS component:", comp[1])
  yl <- paste("MDS component:", comp[2])
  
  if (plotFlag == TRUE)
    plot(ord[, comp], ylab = yl, xlab = xl, ...)
  invisible(ord[, comp])
}

#' Density plots of columns in a matrix
#'
#' Plots the density of the columns of a matrix. Wrapper for \code{\link[quantro]{matdensity}}.
#'
#' @param obj ExpressionSet object
#' @param groups Vector of labels for each sample or a column name of the phenoData slot
#' for the ids to filter. Default is the column names.
#' @param normalized TRUE / FALSE, use the normalized matrix or log2-transformed raw counts
#' @param legendPos Legend title position. If null, does not create legend by default.
#' @param ... Extra parameters for \link[quantro]{matdensity}.
#'
#' @return A density plot for each column in the ExpressionSet object colored by groups
#' @export
#'
#' @importFrom quantro matdensity
#' @importFrom Biobase assayData
#' @importFrom Biobase storageMode
#' @importFrom graphics legend
#'
#' @examples
#' data(skin)
#' filtData <- filterLowGenes(skin,"SMTSD")
#' plotDensity(filtData,groups="SMTSD",legendPos="topleft")
#' # to remove the legend
#' plotDensity(filtData,groups="SMTSD")
#'
plotDensity <- function(obj, groups = NULL, normalized = FALSE,
                        legendPos = NULL, ...) {
  if (length(groups) == 1) {
    groups <- factor(pData(obj)[, groups])
  }
  mat <- extractMatrix(obj, normalized, log = TRUE)
  matdensity(mat, groupFactor = groups, ...)
  if (!is.null(legendPos)) {
    legend(legendPos, legend = levels(groups), fill = 1:length(levels(groups)),
           box.col = NA)
  }
}

#' Plot heatmap of most variable genes
#'
#' This function plots a heatmap of the gene expressions forthe "n" features of interest.
#'
#' @param obj ExpressionSet object or objrix.
#' @param n Number of features to make use of in plotting heatmap.
#' @param fun Function to sort genes by, default \code{\link[stats]{sd}}.
#' @param normalized TRUE / FALSE, use the normalized matrix or raw counts.
#' @param log TRUE/FALSE log2-transform raw counts.
#' @param ... Additional plot arguments for \code{\link[gplots]{heatmap.2}}.
#' @return coordinates
#'
#' @importFrom gplots heatmap.2
#' @importFrom RColorBrewer brewer.pal
#' @importFrom stats sd
#'
#' @export
#' @examples
#' data(skin)
#' tissues <- pData(skin)$SMTSD
#' plotHeatmap(skin,normalized=FALSE,log=TRUE,trace="none",n=10)
#' # Even prettier
#' \donttest{
#' # library(RColorBrewer)
#' data(skin)
#' tissues <- pData(skin)$SMTSD
#' heatmapColColors <- brewer.pal(12,"Set3")[as.integer(factor(tissues))]
#' heatmapCols <- colorRampPalette(brewer.pal(9, "RdBu"))(50)
#' plotHeatmap(skin,normalized=FALSE,log=TRUE,trace="none",n=10,
#'  col = heatmapCols,ColSideColors = heatmapColColors,cexRow = 0.6,cexCol = 0.6)
#'}
plotHeatmap <- function(obj, n = NULL, fun = stats::sd, normalized = TRUE,
                        log = TRUE, ...) {
  if (is.null(n))
    n <- min(nrow(obj), 100)
  mat <- extractMatrix(obj, normalized, log)
  genesToKeep <- which(rowSums(mat) > 0)
  geneStats <- apply(mat[genesToKeep, ], 1, fun)
  geneIndices <- genesToKeep[order(geneStats, decreasing = TRUE)[seq_len(n)]]
  mat <- mat[geneIndices, ]
  heatmap.2(mat, ...)
  invisible(mat)
}

#' Quantile shrinkage normalization
#'
#' This function was modified from github user kokrah.
#'
#' @param obj for counts use log2(raw counts + 1)), for MA use log2(raw intensities)
#' @param groups groups to which samples belong (character vector)
#' @param norm.factors scaling normalization factors
#' @param plot plot weights? (default=FALSE)
#' @param window window size for running median (a fraction of the number of rows of exprs)
#' @param log Whether or not the data should be log transformed before normalization, TRUE = YES.
#'
#' @importFrom stats ave
#' @importFrom graphics par
#' @importFrom graphics abline
#'
#' @return Normalized expression
#'
#' @source \href{https://raw.githubusercontent.com/kokrah/qsmooth/master/R/qsmooth.r}{Kwame Okrah's qsmooth R package}
#' @examples
#' data(skin)
#' head(netZooR:::qsmooth(skin,groups=pData(skin)$SMTSD))
#'
qsmooth <- function(obj, groups, norm.factors = NULL, plot = FALSE,
                    window = 0.05,log=TRUE) {
  stopifnot(class(obj)=="ExpressionSet")
  if(log==TRUE){
    exprs <- log2(exprs(obj)+1)
  } else {
    exprs <- exprs(obj)
  }
  # Stop if exprs contains any NA
  if (any(is.na(exprs)))
    stop("exprs contains NAs (K.Okrah)")
  # Scale normalization step
  if (is.null(norm.factors)) {
    dat <- exprs
  } else {
    dat <- t(t(exprs) - norm.factors)
  }
  # Compute quantile stats
  qs <- qstats(dat, groups, window = window)
  Qref <- qs$Qref
  Qhat <- qs$Qhat
  w <- qs$smoothWeights
  # Weighted quantiles
  normExprs <- w * Qref + (1 - w) * Qhat
  # Re-order normExprs by rank of exprs (columnwise)
  for (i in 1:ncol(normExprs)) {
    # Grab ref. i
    ref <- normExprs[, i]
    # Grab exprs column i
    x <- exprs[, i]
    # Grab ranks of x (using min rank for ties)
    rmin <- rank(x, ties.method = "min")
    # If x has rank ties then average the values of ref at those
    # ranks
    dups <- duplicated(rmin)
    if (any(dups)) {
      # Grab ranks of x (using random ranks for ties) (needed to
      # uniquely identify the indices of tied ranks)
      rrand <- rank(x, ties.method = "random")
      # Grab tied ranks
      tied.ranks <- unique(rmin[dups])
      for (k in tied.ranks) {
        sel <- rrand[rmin == k]  # Select the indices of tied ranks
        ref[sel] <- ave(ref[sel])
      }
    }
    # Re-order ref and replace in normExprs
    normExprs[, i] <- ref[rmin]
  }
  # Plot weights
  if (plot) {
    oldpar <- par(mar = c(4, 4, 1.5, 0.5))
    lq <- length(Qref)
    u <- (1:lq - 0.5)/lq
    if (length(u) > 10000) {
      # do not plot more than 10000 points
      sel <- sample(1:lq, 10000)
      plot(u[sel], w[sel], pch = ".", main = "qsmooth weights",
           xlab = " quantiles", ylab = "Weight", ylim = c(0,
                                                          1))
    } else {
      plot(u, w, pch = ".", main = "qsmooth weights", xlab = "quantiles",
           ylab = "Weight", ylim = c(0, 1))
    }
    abline(h = 0.5, v = 0.5, col = "red", lty = 2)
    par(oldpar)
  }
  rownames(normExprs) <- rownames(exprs)
  colnames(normExprs) <- colnames(exprs)
  normExprs
}

#' Compute quantile statistics
#'
#' This function was directly borrowed from github user kokrah.
#'
#' @param exprs for counts use log2(raw counts + 1)), for MA use log2(raw intensities)
#' @param groups groups to which samples belong (character vector)
#' @param window window size for running median as a fraction on the number of rows of exprs
#'
#' @importFrom stats runmed
#' @importFrom stats model.matrix
#'
#' @return list of statistics
#'
#' @source \href{https://raw.githubusercontent.com/kokrah/qsmooth/master/R/qstats.r}{Kwame Okrah's qsmooth R package}
#' Compute quantile statistics
#'
qstats <- function(exprs, groups, window) {
  # Compute sample quantiles
  Q <- apply(exprs, 2, sort)
  # Compute quantile reference
  Qref <- rowMeans(Q)
  # Compute SST
  SST <- rowSums((Q - Qref)^2)
  # Compute SSB
  f <- factor(as.character(groups))
  X <- model.matrix(~0 + f)
  QBETAS <- t(solve(t(X) %*% X) %*% t(X) %*% t(Q))
  Qhat <- QBETAS %*% t(X)
  SSB <- rowSums((Qhat - Qref)^2)
  # Compute weights
  roughWeights <- 1 - SSB/SST
  roughWeights[SST < 1e-06] <- 1
  # Compute smooth weights
  k <- floor(window * nrow(Q))
  if (k%%2 == 0)
    k <- k + 1
  smoothWeights <- runmed(roughWeights, k = k, endrule = "constant")
  list(Q = Q, Qref = Qref, Qhat = Qhat, QBETAS = QBETAS, SST = SST,
       SSB = SSB, SSE = SST - SSB, roughWeights = roughWeights,
       smoothWeights = smoothWeights)
}

#' Skin RNA-seq data from the GTEx consortium
#'
#' Skin RNA-seq data from the GTEx consortium. V6 release. Random selection of 20 skin samples.
#' 13 of the samples are fibroblast cells, 5 Skin sun exposed, 2 sun unexposed.
#'
#' @docType data
#'
#' @usage data(skin)
#'
#' @format An object of class \code{"ExpressionSet"}; see \code{\link[Biobase]{ExpressionSet}}.
#'
#' @keywords datasets
#'
#' @return ExpressionSet object
#'
#' @references GTEx Consortium, 2015. The Genotype-Tissue Expression (GTEx) pilot analysis: Multitissue gene regulation in humans. Science, 348(6235), pp.648-660.
#' (\href{http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4547484/}{PubMed})
#'
#' @source GTEx Portal
#' @name Skin_data
#' @examples
#' \donttest{data(skin);
#' checkMissAnnotation(skin,"GENDER");}
system('wget https://netzoo.s3.us-east-2.amazonaws.com/netZooR/unittest_datasets/yarn/skin.rdata')
system('mv skin.rdata data/')
"skin"



#' Bladder RNA-seq data from the GTEx consortium
#'
#' Bladder RNA-seq data from the GTEx consortium. V6 release.
#'
#' @docType data
#'
#' @usage data(bladder)
#'
#' @format An object of class \code{"ExpressionSet"}; see \code{\link[Biobase]{ExpressionSet}}.
#'
#' @keywords datasets
#'
#' @references GTEx Consortium, 2015. The Genotype-Tissue Expression (GTEx) pilot analysis: Multitissue gene regulation in humans. Science, 348(6235), pp.648-660.
#' (\href{http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4547484/}{PubMed})
#'
#' @source GTEx Portal
#'
#' @return ExpressionSet object
#' @name Bladder_data
#' @examples
#' \donttest{data(bladder);
#' checkMissAnnotation(bladder);}
system('wget https://netzoo.s3.us-east-2.amazonaws.com/netZooR/unittest_datasets/yarn/bladder.rdata')
system('mv bladder.rdata data/')
"bladder"
