#' Run SEAHORSE in R
#' Author(s): Enakshi Saha
#' Description:
#'               SEAHORSE computes gene-gene coexpression matrix,
#'               associations between a set of given phenotypes and each gene,
#'               performs gene set enrichment analysis (GSEA) for each phenotype,
#'               using the measures of association.
#'               GSEA is performed through R package "fgsea".
#'              
#' @param expression : gene expression matrix (normalized, and filtered) 
#'                     with rows as genes and columns as samples.
#'                     Row and column names must be present.
#'                     Row names must be HGNC symbols.
#'                     Column names must match the row names of the phenotype matrix.
#'                     We assume the data are adequately transformed for use with limma()
#'                     (i.e. log-scaled and/or transformed using VOOM or DeSeq2).
#'                    Users can also provide "NULL" (the default) if a sample-specific set
#'                    of networks is not available.
#' @param network : A weighted matrix where rows are samples and columns are network features.
#'                  Users can also provide "NULL" (the default) if a sample-specific
#'                  set of networks is not available.
#' @param phenotype : phenotype matrix
#'                    with rows as samples and columns as phenotype variables.
#' @param phenotype_dictionary : a vector of strings
#'                               containing type of each phenotype.
#'                               Types can be "dichotomous", "nominal", or "continuous" 
#' @param pathways : a list of pathways (e.g. KEGG, GO, Reactome etc. 
#'                   downloaded from http://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp)
#' @param compute_gene_cor : Whether or not to compute the gene-gene correlation matrix. Default is TRUE.
#' @param compute_phenotype_cor : Whether or not to compute the phenotype-phenotype association matrix. Default is TRUE.
#' @param compute_gene_phenotype_cor : Whether or not to compute the gene-phenotype association matrix. Default is TRUE.
#' @param compute_network_phenotype_cor : Whether or not to compute the network-phenotype association matrix. Default is FALSE.
#' Only computes if a network is provided.
#' @param assoc_method : The method used to infer associations between phenotypes and genes. Default
#' is "spearman". Other options are "pearson", "kendall", and "linear". The "pearson", "kendall", and "spearman" options
#' compute correlations between each phenotype and gene independently for numeric phenotypes and
#' compute an ANOVA for categorical phenotypes. linear" computes a linear regression
#' model of the form gene ~ phenotype1 + phenotype2 + ... + phenotypeN and returns the p-values
#' for each coefficient.
#' @param pval_adj_method Wrapper for p.adjust. Defaults to "none" (no adjustment). Other options are
#' "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", and "fdr".
#' @param usage_report_file Location where the usage report should be stored. Must be RDS format.
#' @param verbose Whether or not to print statements notifying user of run status.
#' Outputs:
#' @return results    : a list containing three objects
#'         results$coexpression: a gene x gene correlation matrix.
#'         results$phenotype_association : a list containing a vector for each phenotype
#'         results$GSEA: a list containing a matrix of GSEA results for each phenotype
#'         results$phenocor: an upper-triangular matrix of phenotype-phenotype associations
#'
#' @examples
#'
#' expression_data = data.frame(matrix(rexp(200, rate=.1), ncol=10, nrow = 20))
#' rownames(expression_data) = paste("gene", 1:20, sep = "")
#' colnames(expression_data) = paste("sample", 1:10, sep = "")
#' 
#' phenotype_data = data.frame(matrix(0, ncol=2, nrow = 10))
#' colnames(phenotype_data) = c("sex", "height")
#' rownames(phenotype_data) = colnames(expression_data)
#' phenotype_data$sex = c(rep("male", nrow(phenotype_data)/2), rep("female", nrow(phenotype_data)/2))
#' phenotype_data$height = 65 + sample.int(10, nrow(phenotype_data), replace = TRUE)
#' 
#' phenotype_dictionary = c("dichotomous", "continuous")
#' 
#' pathways = list()
#' pathways$pathway1 = sample(rownames(expression_data), 5)
#' pathways$pathway2 = sample(rownames(expression_data), 3)
#' pathways$pathway3 = sample(rownames(expression_data), 7)
#'
#' # Run seahorse
#' results <- seahorse(expression = expression_data, phenotype = phenotype_data, 
#'                     phenotype_dictionary = phenotype_dictionary, pathways = pathways)
#'  
#' @export
seahorse <- function(expression = NULL, network = NULL, phenotype, phenotype_dictionary, pathways, compute_gene_cor = TRUE,
                     compute_phenotype_cor = TRUE, compute_gene_phenotype_cor = TRUE,
                     compute_network_phenotype_cor = FALSE,
                     assoc_method = "spearman", pval_adj_method = "none",
                     usage_report_file = NULL, verbose = FALSE){
  
  # Check packages.
  if (!requireNamespace("fgsea", quietly = TRUE)) {
    stop("Package 'fgsea' is required but not installed.")
  }
  if (!requireNamespace("stats", quietly = TRUE)) {
    stop("Package 'stats' is required but not installed.")
  }
  if (!requireNamespace("matrixTests", quietly = TRUE)) {
    stop("Package 'matrixTests' is required but not installed.")
  }
  if (!requireNamespace("peakRAM", quietly = TRUE) && !is.null(usage_report_file)) {
    stop(paste("Package 'peakRAM' is required to monitor memory and time usage.",
               "Install it or set usage_report_file = NULL to turn off monitoring."))
  }
  if (!requireNamespace("limma", quietly = TRUE)) {
    stop("Package 'limma' is required but not installed.")
  }
  if(!requireNamespace("psych", quietly = TRUE)) {
    stop("Package 'psych' is required but not installed.")
  }
  set.seed(0)
  
  # Return an error if pval_adj_method is not an accepted method.
  if(!(pval_adj_method %in% stats::p.adjust.methods)){
    stop(paste(pval_adj_method, "is not a valid method for stats::p.adjust()."))
  }
  
  # Check that association method is valid.
  if(!assoc_method %in% c("pearson", "spearman", "kendall", "linear")){
    stop(paste(assoc_method, "is an invalid association method!"))
  }else if (assoc_method == "linear" && !requireNamespace("limma", quietly = TRUE)) {
    stop("Package 'limma' is required for linear regression. Install it with BiocManager::install('limma')")
  }
  
  # Check that file is valid.
  if(!assoc_method %in% c("pearson", "spearman", "kendall", "linear")){
    stop(paste(assoc_method, "is an invalid association method!"))
  }else if (assoc_method == "linear" && !requireNamespace("limma", quietly = TRUE)) {
    stop("Package 'limma' is required for linear regression. Install it with BiocManager::install('limma')")
  }
  results = list()
  
  # Ensure that the column names are valid for the phenotypic data.
  colnames(phenotype) <- make.names(colnames(phenotype))
  
  # Check that phenotype dictionary is correct.
  for(i in 1:length(phenotype_dictionary)){
    dictType <- phenotype_dictionary[i]
    pheno <- phenotype[,i]
    phenoName <- colnames(phenotype)[i]
    
    # Check that nominal have more than 2 levels.
    levels <- unique(pheno[which(!is.na(pheno))])
    if(dictType == "nominal" && length(levels) <= 2){
      stop(paste("Phenotype", phenoName, "set to nominal but has 2 levels or fewer."))
    }
    # Check that dichotomous don't have more than 2 levels.
    if(dictType == "dichotomous" && length(levels) > 2){
      stop(paste("Phenotype", phenoName, "set to dichotomous but has more than 2 levels."))
    }
    # Check that continuous can be converted to a numeric vector.
    if(dictType == "continuous"){
      phenoChar <- as.character(pheno)
      phenoNum <- suppressWarnings(as.numeric(phenoChar))
      if(length(setdiff(which(is.na(phenoNum)), which(is.na(phenoChar)))) > 0){
        stop(paste("Phenotype", phenoName, "set to continuous but cannot be converted to numeric."))
      }
    }
  }
  
  # Check that usage report file can be created.
  if(!is.null(usage_report_file)){
    tryCatch({
      saveRDS(list(1,2,3), usage_report_file)
      unlink(usage_report_file)
    }, error = function(cond){
      stop(paste(usage_report_file, "could not be created."))
    })
  }
  
  # Check that the network input is valid.
  if(compute_network_phenotype_cor == TRUE){
    if(is.null(network)){
      compute_network_phenotype_cor <- FALSE
      warning(paste("compute_network_phenotype_cor was set to TRUE but no network was provided.",
                    "Will not compute network-phenotype associations."))
    }else if(!is.matrix(network)){
      compute_network_phenotype_cor <- FALSE
      warning(paste("compute_network_phenotype_cor was set to TRUE but network is not a matrix.",
                    "Will not compute network-phenotype associations."))
    }else if(!identical(rownames(phenotype), rownames(network))){
      compute_network_phenotype_cor <- FALSE
      warning(paste("compute_network_phenotype_cor was set to TRUE but network and phenotype row names do not match.",
                    "Will not compute network-phenotype associations."))
    }
  }
  
  # Check that the expression input is valid.
  if(compute_gene_phenotype_cor == TRUE || compute_gene_cor == TRUE){
    if(is.null(expression)){
      compute_gene_phenotype_cor <- FALSE
      compute_gene_cor <- FALSE
      warning(paste("compute_gene_phenotype_cor or compute_gene_cor was set to TRUE but no expression data were provided.",
                    "Will not compute expression associations."))
    }
  }
  
  # Check whether the network has binary or numeric features.
  if(compute_network_phenotype_cor == TRUE){
    netFeatureType <- NA
    if(all(network %in% 0:1, na.rm = TRUE) == TRUE){
      netFeatureType <- "dichotomous"
    }else if(is.numeric(network)){
      netFeatureType <- "continuous"
    }else{
      compute_network_phenotype_cor <- FALSE
      warning(paste("compute_network_phenotype_cor was set to TRUE but network contains a value that is neither binary nor numeric.",
                    "Will not compute network-phenotype associations."))
    }
  }

  # Compute coexpression of genes
  usageGeneGene <- NA
  results$coexpression <- NA
  if(is.null(usage_report_file)){
    if(compute_gene_cor == TRUE){
      if(verbose == TRUE){
        print("Running gene co-expression")
      }
      results$coexpression = cor(t(expression), use="pairwise.complete.obs")
      if(verbose == TRUE){
        print("Gene co-expression complete")
      }
    }
  }else{
    if(compute_gene_cor == TRUE){
      usageGeneGene <- peakRAM::peakRAM({
        if(verbose == TRUE){
          print("Running gene co-expression")
        }
        results$coexpression = cor(t(expression), use="pairwise.complete.obs")
        if(verbose == TRUE){
          print("Gene co-expression complete")
        }
      })
    }
  }
  
  # Compute coexpression of phenotypes
  usagePhenPhen <- NA
  results$phenocor <- NA
  if(is.null(usage_report_file)){
    if(compute_phenotype_cor == TRUE){
      if(verbose == TRUE){
        print("Running phenotype associations")
      }
      results$phenocor = computePhenotypeCorrelations(phenotype = phenotype,
                                                      phenotype_dictionary = phenotype_dictionary,
                                                      method = assoc_method, pval_adj_method = pval_adj_method)
      if(verbose == TRUE){
        print("Phenotype associations complete")
      }
    }
  }else{
    if(compute_phenotype_cor == TRUE){
      usagePhenPhen <- peakRAM::peakRAM({
        if(verbose == TRUE){
          print("Running phenotype associations")
        }
        results$phenocor = computePhenotypeCorrelations(phenotype = phenotype,
                                                        phenotype_dictionary = phenotype_dictionary,
                                                        method = assoc_method, pval_adj_method = pval_adj_method)
        if(verbose == TRUE){
          print("Phenotype associations complete")
        }
      })
    } 
  }
  
  # Compute association of gene expression with phenotypes and run GSEA
  usagePhenGene <- NA
  results$phenotype_association = list()
  results$GSEA = list()
  if(is.null(usage_report_file)){
    if(compute_gene_phenotype_cor == TRUE){
      if(verbose == TRUE){
        print("Running gene-phenotype associations")
      }
      if(assoc_method %in% c("pearson", "spearman", "kendall")){
        corr <- computeCorrelations(expression = expression, phenotype = phenotype,
                                    pathways = pathways, phenotype_dictionary = phenotype_dictionary,
                                    method = assoc_method, pval_adj_method = pval_adj_method)
        results$phenotype_association <- corr$pheno
        results$GSEA <- corr$gsea
      }else{
        linReg <- computeLinearRegression(expression = expression, phenotype = phenotype,
                                          pathways = pathways, phenotype_dictionary = phenotype_dictionary,
                                          pval_adj_method = pval_adj_method)
        results$phenotype_association <- linReg$pheno
        results$GSEA <- linReg$gsea
      }
      if(verbose == TRUE){
        print("Gene-phenotype associations complete")
      }
    }
  }else{
    if(compute_gene_phenotype_cor == TRUE){
      usagePhenGene <- peakRAM::peakRAM({
        if(verbose == TRUE){
          print("Running gene-phenotype associations")
        }
        if(assoc_method %in% c("pearson", "spearman", "kendall")){
          corr <- computeCorrelations(expression = expression, phenotype = phenotype,
                                      pathways = pathways, phenotype_dictionary = phenotype_dictionary,
                                      method = assoc_method, pval_adj_method = pval_adj_method)
          results$phenotype_association <- corr$pheno
          results$GSEA <- corr$gsea
        }else{
          linReg <- computeLinearRegression(expression = expression, phenotype = phenotype,
                                            pathways = pathways, phenotype_dictionary = phenotype_dictionary,
                                            pval_adj_method = pval_adj_method)
          results$phenotype_association <- linReg$pheno
          results$GSEA <- linReg$gsea
        }
        if(verbose == TRUE){
          print("Gene-phenotype associations complete")
        }
      })
    }
  }
  
  # Compute association of network features with phenotypes.
  usagePhenNet <- NA
  results$phenotype_net_association = list()
  if(is.null(usage_report_file)){
    if(compute_network_phenotype_cor == TRUE){
      if(verbose == TRUE){
        print("Running network-phenotype associations")
      }
      if(assoc_method %in% c("pearson", "spearman", "kendall")){
        corr <- computeCorrelationsNetwork(network = network, phenotype = phenotype,
                                           phenotype_dictionary = phenotype_dictionary,
                                           method = assoc_method, pval_adj_method = pval_adj_method,
                                           featureType = netFeatureType)
        results$phenotype_net_association <- corr
      }else{
        linReg <- computeLinearRegressionNetwork(network = network, phenotype = phenotype,
                                          phenotype_dictionary = phenotype_dictionary,
                                          pval_adj_method = pval_adj_method,
                                          featureType = netFeatureType)
        results$phenotype_net_association <- linReg
      }
      if(verbose == TRUE){
        print("Network-phenotype associations complete")
      }
    }
  }else{
    if(compute_network_phenotype_cor == TRUE){
      usagePhenNet <- peakRAM::peakRAM({
        if(verbose == TRUE){
          print("Running network-phenotype associations")
        }
        if(assoc_method %in% c("pearson", "spearman", "kendall")){
          corr <- computeCorrelationsNetwork(network = network, phenotype = phenotype,
                                      phenotype_dictionary = phenotype_dictionary,
                                      method = assoc_method, pval_adj_method = pval_adj_method,
                                      featureType = netFeatureType)
          results$phenotype_net_association <- corr
        }else{
          linReg <- computeLinearRegressionNetwork(network = network, phenotype = phenotype,
                                            phenotype_dictionary = phenotype_dictionary,
                                            pval_adj_method = pval_adj_method,
                                            featureType = netFeatureType)
          results$phenotype_net_association <- linReg
        }
        if(verbose == TRUE){
          print("Network-phenotype associations complete")
        }
      })
    }
  }
  
  # Save the results to a usage file.
  if(!is.null(usage_report_file)){
    saveRDS(list(GeneToGeneCor = usageGeneGene, PhenToPhenCor = usagePhenPhen,
                 PhenToGeneCor = usagePhenGene, PhenToNetCor = usagePhenNet), usage_report_file)
  }
  
  # Return the results.
  return(results)
}

#' Computes linear regression for both numeric and categorical phenotypes. Returns p-values
#' computed from moderated t-statistics using ebayes() and uses t-statistics as input to the
#' GSEA analysis.
#' @param expression : gene expression matrix (normalized, and filtered) 
#'                     with rows as genes and columns as samples.
#'                     Row and column names must be present.
#'                     Row names must be HGNC symbols.
#'                     Column names must match the row names of the phenotype matrix.
#' @param phenotype : phenotype matrix
#'                    with rows as samples and columns as phenotype variables.
#' @param phenotype_dictionary : a vector of strings
#'                               containing type of each phenotype.
#'                               Types can be either "numeric" or "categorical" 
#' @param pathways : a list of pathways (e.g. KEGG, GO, Reactome etc. 
#'                   downloaded from http://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp)
#' @param pval_adj_method Wrapper for p.adjust. Defaults to "none" (no adjustment). Other options are
#' "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", and "fdr".
#' This parameter is only used for the "linear" association method. Use if your data are raw RNA-seq
#' counts. Default is TRUE.
computeLinearRegression <- function(expression, phenotype, phenotype_dictionary, pathways,
                                    pval_adj_method){
  
  # Initialize output.
  phenoAssoc <- list
  gsea <- list

  # Compute a design matrix for the phenotypic data.
  formula <- paste("~ ", paste(colnames(phenotype)[1:(ncol(phenotype) - 1)], 
                                   collapse = " + "), colnames(phenotype)[ncol(phenotype)],
                   sep = " + ")
  design <- model.matrix(object = stats::as.formula(formula), data = phenotype)
  
  # Check whether samples were removed and modify gene expression data accordingly.
  if(length(setdiff(colnames(expression), rownames(design)) > 0)){
    expressionSampsOld <- ncol(expression)
    expression <- expression[,rownames(design)]
    expressionSampsNew <- ncol(expression)
    message(paste("Out of", expressionSampsOld, "we retained", expressionSampsNew,
                  "samples."))
  }

  # Run the linear models.
  fit <- limma::lmFit(object = expression, design = design)

  # Compute empirical Bayes statistics.
  bayes <- limma::eBayes(fit = fit)

  # Extract p-values and t-statistics.
  t_stats <- bayes$t
  p_vals  <- bayes$p.value
  
  # Format p-values as a list for all covariates.
  p_list <- lapply(colnames(p_vals), function(c){
    p <- p_vals[,c]
    t <- t_stats[,c]
    padj <- stats::p.adjust(p_vals[,c], method = pval_adj_method)
    return(data.frame(stat = t, pval = p, padj = padj, testType = "linear", row.names = rownames(p_vals)))
  })
  names(p_list) <- colnames(p_vals)

  phenoAssoc = p_list

  # Run GSEA
  gsea_list <- lapply(colnames(t_stats), function(c){
    tscore <- t_stats[,c]
    t_rank = sort(tscore, decreasing = T)
    fgseaRes <- fgsea::fgsea(pathways = pathways, stats = t_rank, minSize=15, maxSize=500)
    return(fgseaRes)
  })
  names(gsea_list) <- colnames(t_stats)
  
  gsea = gsea_list
  return(list(pheno = phenoAssoc, gsea = gsea))
}

#' Computes linear regression for both numeric and categorical phenotypes with respect to network. 
#' Computes logistic regression when the outcome (numeric) is dichotomous.
#' @param network : A weighted matrix where rows are samples and columns are network features.
#' @param phenotype : phenotype matrix
#'                    with rows as samples and columns as phenotype variables.
#' @param phenotype_dictionary : a vector of strings
#'                               containing type of each phenotype.
#'                               Types can be either "numeric" or "categorical" 
#' @param pval_adj_method Wrapper for p.adjust. Defaults to "none" (no adjustment). Other options are
#' "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", and "fdr".
#' @param featureType "continuous" or "dichotomous" (for binary networks)
computeLinearRegressionNetwork <- function(network, phenotype, phenotype_dictionary,
                                    pval_adj_method, featureType){
  
  # Initialize output.
  phenoAssoc <- list

  # Compute a design matrix for the phenotypic data.
  formula <- paste("~ ", paste(colnames(phenotype)[1:(ncol(phenotype) - 1)], 
                               collapse = " + "), colnames(phenotype)[ncol(phenotype)],
                   sep = " + ")
  design <- model.matrix(object = stats::as.formula(formula), data = phenotype)

  # Check whether samples were removed and modify gene expression data accordingly.
  if(length(setdiff(rownames(network), rownames(design)) > 0)){
    networkSampsOld <- nrow(network)
    network <- networkSampsOld[rownames(design),]
    networkSampsNew <- nrow(network)
    message(paste("Out of", networkSampsOld, "we retained", networkSampsNew,
                  "samples."))
  }
  
  # Initialize pvals.
  pval <- NA
  thisstat <- NA
  thisTestType <- NA
  coef <- NA
  
  # Do linear or logistic regression depending on the network feature type.
  if(featureType == "continuous"){
    
    # Run the linear models. We can still use LIMMA, because the values aren't
    # adjusted until we do empirical Bayes adjustment.
    fit <- limma::lmFit(object = t(network), design = design)

    # Extract p-values and t-statistics.
    coef <- fit$coefficients
    se <- fit$stdev.unscaled * fit$sigma
    tstat <- coef / se
    pval <- 2 * pt(
      -abs(tstat),
      df = fit$df.residual
    )
    thisstat <- tstat
    thisTestType <- rep("linear", nrow(tstat))
  }else{
    
    # Run the linear models. We need to use GLM here because these are logistic regression models.
    resList <- lapply(colnames(network), function(feat){
      concatMat <- cbind(phenotype, network[,feat])
      colnames(concatMat)[ncol(concatMat)] <- "feat"
      fitFeat <- glm(formula = paste("feat", formula), data = concatMat, family = "binomial")
      zStatFeat <- coef(summary(fitFeat))[,3]
      pvalFeat <- coef(summary(fitFeat))[,4]
      fitFeatDF <- data.frame(z = zStatFeat, p = pvalFeat)
      rownames(fitFeatDF) <- rownames(coef(summary(fitFeat)))
      return(fitFeatDF)
    })
    pval <- do.call(rbind, lapply(1:length(resList), function(i){
      return(resList[[i]]$p)
    }))
    thisstat <- do.call(rbind, lapply(1:length(resList), function(i){
      return(resList[[i]]$z)
    }))
    rownames(pval) <- colnames(network)
    rownames(thisstat) <- colnames(network)
    colnames(pval) <- rownames(resList[[1]])
    colnames(thisstat) <- rownames(resList[[1]])
    thisTestType <- rep("logistic", nrow(thisstat))
  }

  # Format p-values as a list for all covariates.
  p_list <- lapply(colnames(pval), function(c){
    p <- pval[,c]
    stat <- thisstat[,c] 
    padj <- stats::p.adjust(pval[,c], method = pval_adj_method)
    return(data.frame(stat = stat, pval = p, padj = padj, testType = thisTestType, row.names = rownames(pval)))
  })
  names(p_list) <- colnames(pval)
  phenoAssoc = p_list

  return(phenoAssoc)
}

#' Computes gene-phenotype correlation for continuous phenotypes, ANOVA for nominal
#' phenotypes, and t-test for dichotomous phenotypes
#' @param expression : gene expression matrix (normalized, and filtered) 
#'                     with rows as genes and columns as samples.
#'                     Row and column names must be present.
#'                     Row names must be HGNC symbols.
#'                     Column names must match the row names of the phenotype matrix.
#' @param phenotype : phenotype matrix
#'                    with rows as samples and columns as phenotype variables.
#' @param pathways : Pathways used for GSEA calculation
#' @param phenotype_dictionary : a vector of strings
#'                               containing type of each phenotype.
#'                               Types can be either "numeric" or "categorical" 
#' @param method : One of "pearson", "spearman", or "kendall".
#' @param pval_adj_method Wrapper for p.adjust. Defaults to "none" (no adjustment). Other options are
#' "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", and "fdr".
computeCorrelations <- function(expression, phenotype, phenotype_dictionary, pathways, method,
                                pval_adj_method){
  
  phenoAssoc <- list()
  gsea <- list()
  output_seahorse <- data.frame(stat = NA, pval = NA, testType = NA)
  for (i in 1:ncol(phenotype)){
    pheno = phenotype[,i]
    pheno_name = colnames(phenotype)[i]

    if (phenotype_dictionary[i] == "continuous"){
      output_seahorse = gsea_continuous(expression, pheno, pathways, method = method)
    }else if (phenotype_dictionary[i] == "nominal") {
      output_seahorse = gsea_nominal(expression, pheno, pathways)
    }else{
    	output_seahorse = gsea_dichotomous(expression, pheno, pathways)
    }
    output_seahorse_padj <- stats::p.adjust(output_seahorse$p, method = pval_adj_method)
    phenoAssoc[[pheno_name]] = data.frame(stat = output_seahorse$stat,
                                          pval = output_seahorse$p,
                                          padj = output_seahorse_padj,
                                          testType = output_seahorse$testType,
                                          row.names = names(output_seahorse$stat))
    gsea[[pheno_name]] = output_seahorse$GSEA
  }
  return(list(pheno = phenoAssoc, gsea = gsea))
}

#' Computes network-phenotype correlation for continuous phenotypes with continuous network features,
#' t-test for continuous phenotypes with dichotomous network features, ANOVA for categorical
#' phenotypes with continuous network features, t-test for dichotomous phenotypes with continuous network features,
#' and FFH for categorical or dichotomous phenotypes with dichotomous network features.
#' @param network : A weighted matrix where rows are samples and columns are network features.
#' @param phenotype : phenotype matrix
#'                    with rows as samples and columns as phenotype variables.
#' @param phenotype_dictionary : a vector of strings
#'                               containing type of each phenotype.
#'                               Types can be either "numeric" or "categorical" 
#' @param pathways : Pathways used for GSEA calculation
#' @param method : One of "pearson", "spearman", or "kendall".
#' @param pval_adj_method : Wrapper for p.adjust. Defaults to "none" (no adjustment). Other options are
#' "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", and "fdr".
#' @param featureType : Whether network features are continuous or dichotomous.
computeCorrelationsNetwork <- function(network, phenotype, phenotype_dictionary, pathways, method,
                                pval_adj_method, featureType){
  
  phenoAssoc <- list()
  for (i in 1:ncol(phenotype)){
    pheno = phenotype[,i]
    pheno_name = colnames(phenotype)[i]
    
    if (phenotype_dictionary[i] == "continuous" && featureType == "continuous"){
      output_seahorse = networkContinuous(network, pheno, method = method)
    }else if (phenotype_dictionary[i] == "nominal" && featureType == "continuous") {
      output_seahorse = networkANOVA(network, pheno)
    }else if (phenotype_dictionary[i] == "dichotomous" && featureType == "continuous") {
      output_seahorse = networkTTest(network, pheno, dichotomousVar = "phenotype")
    }else if (phenotype_dictionary[i] == "continuous" && featureType == "dichotomous") {
      output_seahorse = networkTTest(network, pheno, dichotomousVar = "feature")
    }else if ((phenotype_dictionary[i] == "nominal" || phenotype_dictionary[i] == "dichotomous") && featureType == "dichotomous") {
      output_seahorse = networkChisq(network, pheno)
    }
    output_seahorse_padj <- stats::p.adjust(output_seahorse$p, method = pval_adj_method)
    phenoAssoc[[pheno_name]] = data.frame(stat = output_seahorse$stat,
                                          pval = output_seahorse$p,
                                          padj = output_seahorse_padj,
                                          testType = output_seahorse$testType,
                                          row.names = names(output_seahorse$stat))
  }
  return(phenoAssoc)
}

#' Function to find network correlations with continuous phenotypes.
#' @param network : A weighted matrix where rows are samples and columns are network features.
#' @param pheno : phenotype matrix
#'                    with rows as samples and columns as phenotype variables.
#' @param method : One of "pearson", "spearman", or "kendall".
#' @export
networkContinuous <- function(network, pheno, method){

  # Run correlations.
  # We use corr.test from the psych package rather than cor.test because cor.test
  # only allows vectorized input.
  phenotype_vector = as.numeric(pheno)
  cors <- psych::corr.test(network, phenotype_vector, use = "pairwise", method = method)
  corvals <- as.numeric(cors$r)
  pvals <- as.numeric(cors$p)
  names(corvals) <- colnames(network)
  names(pvals) <- colnames(network)
  
  return(data.frame(stat = corvals,
                    pvals = pvals,
                    V = rep(NA, ncol(network)),
                    testType = rep("Cor", ncol(network)),
                    row.names = colnames(network)))
}

#' Function to run ANOVA on nominal phenotypes with continuous network features.
#' @param network : A weighted matrix where rows are samples and columns are network features.
#' @param phenotype : phenotype matrix
#'                    with rows as samples and columns as phenotype variables.
#' @export
networkANOVA <- function(network, phenotype){
  
  phenotype_vector = factor(as.character(phenotype))
  anovares <- do.call(rbind, lapply(colnames(network), function(x){
    results <- data.frame(stat = NA, p = NA)
    tryCatch({
      res <- anova(lm(as.numeric(as.character(network[,x]))~phenotype_vector))
      results$p <- res$`Pr(>F)`[1]
      results$stat <- res$`F value`[1]
    }, error = function(cond){
      warning("In this phenotype pair, all continuous values are missing for all but one phenotype level - NA result will be returned")
    })
    return(results)
  }))
  
  return(data.frame(stat = anovares$stat,
                    p = anovares$p,
                    V = rep(NA, ncol(network)),
                    testType = rep("ANOVA", ncol(network)),
                    row.names = colnames(network)))
}

#' Function to compute Welch's t-test statistics on dichotomous phenotypes with
#' continuous network features or on continuous phenotypes with dichotomous network features.
#' @param network The network features to test.
#' @param phenotype A phenotype vector against which to run associations.
#' @param dichotomousVar Whether the phenotype is dichotomous or continuous.
#' @returns A data frame with two vectors: "cor", which lists the t-test p-values, 
#' and "V" which is set to NA.
networkTTest <- function(network, phenotype, dichotomousVar){
  
  # Case 1 - the phenotype is dichotomous and the network features are numeric.
  # We can vectorize this and use matrixTests.
  # Case 2 - the phenotype is numeric and the network features are dichotomous.
  # We cannot vectorize this and use t.test instead, which defaults to Welch's t-test.
  if(dichotomousVar == "phenotype"){
    
    # Split groups.
    phenotype_vector = factor(as.character(phenotype))
    levels <- unique(phenotype_vector)
    group1 <- network[which(phenotype_vector == levels[1]),]
    group2 <- network[which(phenotype_vector == levels[2]),]
    
    # Check that thresholds are met.
    whichLevel1 <- length(which(phenotype_vector == levels[1]))
    whichLevel2 <- length(which(phenotype_vector == levels[2]))
    tres <- data.frame(stat = NA, pval = NA)
    
    if(length(dim(group1)) >= 2 || length(dim(group2)) >= 2){
      meetsThreshold <- colSums(!is.na(group1)) >= 2 & colSums(!is.na(group2)) >= 2
      
      # Initialize result to NA.
      tres <- data.frame(stat = rep(NA, ncol(group1)), pval = rep(NA, ncol(group1)))
      rownames(tres) <- colnames(group1)
      
      # Only compute correlations if both levels are represented in the phenotype vector.
      if(whichLevel1 > 0 && whichLevel2 > 0 && length(which(meetsThreshold == TRUE)) > 0){
        # Compute t-test where thresholds are met.
        tresValid <- matrixTests::col_t_welch(
          group1[, meetsThreshold, drop = FALSE],
          group2[, meetsThreshold, drop = FALSE]
        )
        
        # Compute remaining t-tests.
        tres[meetsThreshold, "pval"] <- tresValid$pvalue
        tres[meetsThreshold, "stat"] <- tresValid$statistic
      }
    }else{
      meetsThreshold <- sum(!is.na(group1)) >= 2 & sum(!is.na(group2)) >= 2
      
      # Initialize result to NA.
      rownames(tres) <- names(group1)
      
      # Only compute correlations if both levels are represented in the phenotype vector.
      if(whichLevel1 > 0 && whichLevel2 > 0 && length(which(meetsThreshold == TRUE)) > 0){
        tryCatch({
          # Compute t-test where thresholds are met. Welch's t-test is computed by default.
          tresValid <- t.test(group1[meetsThreshold, drop = FALSE], 
                              group2[meetsThreshold, drop = FALSE])
          
          # Compute remaining t-tests.
          tres[meetsThreshold, "pval"] <- tresValid$p.value
          tres[meetsThreshold, "stat"] <- tresValid$statistic
        }, error = function(cond){
          warning("Could not compute t-test (it is possible that variance is too low). Returning NA.")
        })
      }
    }
    
  }else{
    tres <- do.call(rbind, lapply(1:ncol(network), function(i){
      network_vector = factor(as.character(network[,i]))
      levels <- unique(network_vector)
      group1 <- phenotype[which(network_vector == levels[1])]
      group2 <- phenotype[which(network_vector == levels[2])]
      result <- data.frame(stat = NA, p = NA)
      if(length(which(!is.na(group1))) > 2 && length(which(!is.na(group2))) > 2){
        tryCatch({
          testres <- t.test(group1, group2)
          result$stat <- testres$statistic
          result$p <- testres$p.value
        }, error = function(cond){
          warning("Could not compute t-test (it is possible that variance is too low). Returning NA.")
        })
      }
      rownames(result) <- colnames(network)[i]

      return(result)
    }))
  }
  if(length(which(is.na(tres$stat))) > 0){
    warning("Some network features did not have sufficient sample sizes to perform a t-test - NAs will be returned")
  }
  
  # Return the data frame.
  return(data.frame(stat = tres$stat,
                    p = tres$p,
                    V = rep(NA, ncol(network)),
                    testType = rep("T-Test", ncol(network)),
                    row.names = colnames(network)))
  
}

#' Function to run associations between categorical phenotypes and dichotomous network features.
#' If at least 5 samples exist in each group, run a Chi-square test and compute Cramer's V.
#' Otherwise, run a Fisher-Freeman-Halton Test.
#' @param network The network features to test.
#' @param phenotype The phenotype of interest
#' @returns A data frame with two vectors: "cor", which lists the Chi-square or Fisher-Freeman-Halton Test p-values, 
#' and "V" which lists the Cramer's V statistics (will be NA if Fisher-Freeman-Halton Test is computed)
networkChisq <- function(network, phenotype){
  
  # Generate the tables for all phenotype pairs.
  allPhenoPairTables <- lapply(1:ncol(network), function(i){
    grpCounts <- table(network[,i], phenotype)
  })
  names(allPhenoPairTables) <- colnames(network)
  
  # Do a chi-square test everywhere else.
  # Run the comparisons one at a time because Chi-square and FFH cannot be scaled.
  chisqRes <- lapply(allPhenoPairTables, function(chisqTable) {
    
    # Check that we have 2 or more non-zero marginals. Only run the test if we do.
    rowMarginals <- rowSums(chisqTable)
    colMarginals <- colSums(chisqTable) 
    nonzeroRowMarginalCount <- length(which(rowMarginals > 0))
    nonzeroColMarginalCount <- length(which(colMarginals > 0))
    pval <- NA
    stat <- NA
    V <- NA
    type <- "Chi-square"
    if(nonzeroRowMarginalCount >= 2 && nonzeroColMarginalCount >= 2){
      # Run the test. If a warning is thrown, switch to FFH.
      chisq <- withCallingHandlers(
        chisq.test(chisqTable),
        warning = function(w) {
          if (grepl("Chi-squared approximation may be incorrect", w$message)) {
            # Do not change only within the warning scope but also within the parent scope.
            type <<- "FFH"
            invokeRestart("muffleWarning")
          }
        }
      )
      
      # Get the p-value and Cramer's V.
      pval <- chisq$p.value
      stat <- chisq$statistic
      V <- sqrt(chisq$statistic / (sum(chisqTable) * min(nrow(chisqTable) - 1, 
                                                         ncol(chisqTable) - 1)))
    }else{
      warning("Less than 2 nonzero marginals in the contingency table - Chi-square will return NA")
    }
    return(list(pval = pval, stat = stat, v = V, testType = type))
  })
  chisqP <- unlist(lapply(chisqRes, function(res){
    return(res$pval)
  }))
  chisqStat <- unlist(lapply(chisqRes, function(res){
    return(res$stat)
  }))
  cramerV <- unlist(lapply(chisqRes, function(res){
    return(res$v)
  }))
  testType <- unlist(lapply(chisqRes, function(res){
    return(res$testType)
  }))
  corResNotLessFull <- data.frame(stat = chisqStat,
                                  pval = chisqP,
                                  V = cramerV,
                                  testType = testType,
                                  row.names = colnames(network))
  
  # Run FFS where warning was issued.
  corResNotLess <- corResNotLessFull[which(corResNotLessFull$testType == "Chi-square"),]
  switchedTables <- allPhenoPairTables[which(corResNotLessFull$testType == "FFH")]
  corResSwitched <- phenotype_ffs(switchedTables)
  
  # Bind together the results from the FFH and Chi-Square tests.
  corRes <- rbind(corResNotLess, corResSwitched)
  corRes <- corRes[colnames(network),]
  
  return(corRes)
}

#' Function to run GSEA for a continuous phenotype
#' @param expression : gene expression matrix (normalized, and filtered) 
#'                     with rows as genes and columns as samples.
#'                     Row and column names must be present.
#'                     Row names must be HGNC symbols.
#'                     Column names must match the row names of the phenotype matrix.
#' @param pheno : phenotype matrix
#'                    with rows as samples and columns as phenotype variables.
#' @param pathways : a list of pathways (e.g. KEGG, GO, Reactome etc. 
#'                   downloaded from http://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp)
#' @param method : One of "pearson", "spearman", or "kendall".
#' @export
gsea_continuous <- function(expression, pheno, pathways, method){

  output_seahorse <- list()
  output_seahorse$stat <- list()
  output_seahorse$p <- list()
  output_seahorse$GSEA <- list()
  output_seahorse$testType <- list()
  
  # Run correlations.
  phenotype_vector = as.numeric(pheno)
  cors <- psych::corr.test(t(expression), phenotype_vector, use = "pairwise", method = method)
  corvals <- as.numeric(cors$r)
  pvals <- as.numeric(cors$p)
  names(corvals) <- rownames(expression)
  names(pvals) <- rownames(expression)
  output_seahorse$stat <- corvals
  output_seahorse$p <- pvals
  output_seahorse$testType <- rep("Cor", length(pvals))
  
  # Run GSEA
  cor_rank = sort(output_seahorse$stat, decreasing = TRUE)
  fgseaRes <- fgsea::fgsea(pathways, cor_rank, minSize=15, maxSize=500)
  output_seahorse$GSEA = fgseaRes
  
  return(output_seahorse)
}

#' Function to run GSEA for a nominal phenotype
#' @param expression : gene expression matrix (normalized, and filtered) 
#'                     with rows as genes and columns as samples.
#'                     Row and column names must be present.
#'                     Row names must be HGNC symbols.
#'                     Column names must match the row names of the phenotype matrix.
#' @param pheno : phenotype matrix
#'                    with rows as samples and columns as phenotype variables.
#' @param pathways : a list of pathways (e.g. KEGG, GO, Reactome etc. 
#'                   downloaded from http://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp)
#'                   
#' @export
gsea_nominal <- function(expression, pheno, pathways){

  output_seahorse <- list()
  output_seahorse$stat <- list()
  output_seahorse$testType <- list()
  output_seahorse$p <- list()
  output_seahorse$GSEA <- list()
  
  # Run the linear models.
  phenotype_vector <- factor(as.character(pheno))
  # Remove NA values.
  hasVal <- which(!is.na(phenotype_vector))
  phenotype_vector <- phenotype_vector[hasVal]
  expression <- expression[,hasVal]
  design <- model.matrix(~ phenotype_vector)
  fit <- limma::lmFit(object = expression, design = design)
  
  # Compute empirical Bayes statistics.
  # If there are no residual degrees of freedom, catch it and do not return a value.
  tryCatch({
    bayes <- limma::eBayes(fit = fit)
    
    coef_to_test <- setdiff(colnames(design), "(Intercept)")
    anova_res <- limma::topTable(bayes, coef = coef_to_test, number = Inf, sort.by = "none")
    
    # Format p-values as a list for all covariates.
    p <- anova_res[,"P.Value"]
    output_seahorse$p <- p
    names(output_seahorse$p) <- rownames(anova_res)
    output_seahorse$stat <- anova_res[,"F"]
    names(output_seahorse$stat) <- rownames(anova_res)
    output_seahorse$testType <- rep("ANOVA", nrow(anova_res))
    
    # Run GSEA
    statSorted <- sort(-1 * log10(output_seahorse$p), decreasing = TRUE)
    fgseaRes <- fgsea::fgsea(pathways, statSorted, minSize=15, maxSize=500, scoreType = "pos")
    output_seahorse$GSEA = fgseaRes
  }, error = function(cond){
    warning("Could not compute empirical Bayes statistics for this phenotypic variable. Returning empty list.")
  })
  
  return(output_seahorse)
}

#' Function to run GSEA for a dichotomous phenotype
#' @param expression : gene expression matrix (normalized, and filtered) 
#'                     with rows as genes and columns as samples.
#'                     Row and column names must be present.
#'                     Row names must be HGNC symbols.
#'                     Column names must match the row names of the phenotype matrix.
#' @param pheno : phenotype matrix
#'                    with rows as samples and columns as phenotype variables.
#' @param pathways : a list of pathways (e.g. KEGG, GO, Reactome etc. 
#'                   downloaded from http://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp)
#' @export
gsea_dichotomous <- function(expression, pheno, pathways){
  
  output_seahorse = list()
  output_seahorse$stat = list()
  output_seahorse$p = list()
  output_seahorse$GSEA = list()
  output_seahorse$testType <- list()
  
  # Run the linear models.
  phenotype_vector = factor(as.character(pheno))

  # Remove NA values.
  hasVal <- which(!is.na(phenotype_vector))
  phenotype_vector <- phenotype_vector[hasVal]
  expression <- expression[,hasVal]
  
  # Check that we still have multiple values. If not, return NA for this covariate.
  output_seahorse$stat <- NA
  output_seahorse$p <- NA
  output_seahorse$testType <- NA
  output_seahorse$GSEA <- NA
  if(length(unique(phenotype_vector)) > 1){
    design <- model.matrix(~ phenotype_vector)
    fit <- limma::lmFit(object = expression, design = design)
    
    # Compute empirical Bayes statistics.
    tryCatch({
      bayes <- limma::eBayes(fit = fit)

      coef_to_test <- setdiff(colnames(design), "(Intercept)")
      t_res <- limma::topTable(bayes, coef = coef_to_test, number = Inf, sort.by = "none")
      
      # Format p-values as a list for all covariates.
      p <- t_res[,"P.Value"]
      t <- t_res[,"t"]
      output_seahorse$stat <- t
      output_seahorse$p <- p
      output_seahorse$testType <- rep("LIMMA Moderated t-test", nrow(t_res))
      names(output_seahorse$p) <- rownames(t_res)
      names(output_seahorse$stat) <- rownames(t_res)

      # Run GSEA
      fgseaRes <- NA
      tryCatch({
        statSorted <- sort(output_seahorse$stat, decreasing = TRUE)
        fgseaRes <- fgsea::fgsea(pathways, statSorted, minSize=15, maxSize=500)
      }, error = function(cond){
        print(cond)
      })
      output_seahorse$GSEA = fgseaRes
    }, error = function(cond){
      print(cond)
      warning("Could not compute empirical Bayes statistics for this phenotypic variable. Returning empty list.")
    })
  }else{
    warning("Phenotype has only one level. Returning NA for all gene associations.")
  }
  return(output_seahorse)
}


#' Computes phenotype-phenotype correlation for continuous phenotypes, ANOVA for nominal
#' phenotypes, and t-test for dichotomous phenotypes
#' @param phenotype : phenotype matrix
#'                    with rows as samples and columns as phenotype variables.
#' @param phenotype_dictionary : a vector of strings
#'                               containing type of each phenotype.
#'                               Types can be either "numeric" or "categorical" 
#' @param method : One of "pearson", "spearman", or "kendall".
#' @param pval_adj_method Wrapper for p.adjust. Defaults to "none" (no adjustment). Other options are
#' "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", and "fdr".
computePhenotypeCorrelations <- function(phenotype, phenotype_dictionary, method,
                                pval_adj_method){

  # For each phenotype, compute all of its associations.
  phenoAssoc <- lapply(1:(ncol(phenotype)-1), function(i){
    # Extract the phenotype.
    pheno = phenotype[,i]
    pheno_name = colnames(phenotype)[i]
    phenoType <- phenotype_dictionary[i]
    
    # For dichotomous phenotypes, calculate Chi-squared + Cramer's V OR
    # Fisher-Freeman-Halton Test for all categorical comparisons and t-test for
    # all continuous comparisons.
    #
    # For nominal phenotypes, calculate calculate Chi-squared + Cramer's V OR
    # Fisher-Freeman-Halton Test for all categorical comparisons and ANOVA for
    # all continuous comparisons.
    #
    # For continuous phenotypes, calculate t-test for all dichotomous comparisons,
    # ANOVA for all nominal comparisons, and correlation for all continuous comparisons.
    categorical <- c("dichotomous", "nominal")
    output_seahorse <- NULL
    
    # Only compute upper triangular matrix.
    rangeAfter <- (i+1):ncol(phenotype)
    phenosAfter <- as.data.frame(phenotype[,rangeAfter])
    colnames(phenosAfter) <- colnames(phenotype)[rangeAfter]
    typesAfter <- phenotype_dictionary[rangeAfter]
    
    # Initialize an empty data frame.
    emptyDF <- data.frame(stat = c(),
                          p = c(),
                          V = c(),
                          padj = c(),
                          testType = c(),
                          row.names = c())
    emptyPadj <- c()
    
    # Run associations.
    if(phenoType == "dichotomous"){
      
      # Do categorical phenotypes.
      whichCat <- which(typesAfter %in% categorical)
      phenoCat <- as.data.frame(phenosAfter[,whichCat])
      colnames(phenoCat) <- colnames(phenosAfter)[whichCat]
      rownames(phenoCat) <- rownames(phenosAfter)
      output_seahorse_cat <- emptyDF
      output_seahorse_padj_cat <- emptyPadj
      if(ncol(phenoCat) > 0){
        output_seahorse_cat <- phenotype_chisq(phenotype = pheno, 
                                               phenotypesToCompare = phenoCat, 
                                               phenotypeType = phenoType)
        output_seahorse_padj_cat <- rep(NA, nrow(output_seahorse_cat))
        whichFFH <- which(output_seahorse_cat$testType == "FFH")
        whichChisq <- which(output_seahorse_cat$testType == "Chi-square")
        if(length(whichFFH) > 0){
          output_seahorse_padj_cat[whichFFH] <- stats::p.adjust(output_seahorse_cat[whichFFH, "pval"], method = pval_adj_method)
        }
        if(length(whichChisq) > 0){
          output_seahorse_padj_cat[whichChisq] <- stats::p.adjust(output_seahorse_cat[whichChisq, "pval"], method = pval_adj_method)
        }
      }
      
      # Do continuous phenotypes.
      whichContinuous <- which(typesAfter == "continuous")
      phenoCon <- as.data.frame(phenosAfter[,whichContinuous])
      colnames(phenoCon) <- colnames(phenosAfter)[whichContinuous]
      rownames(phenoCon) <- rownames(phenosAfter)
      output_seahorse_con <- emptyDF
      output_seahorse_padj_con <- emptyPadj
      if(ncol(phenoCon) > 0){
        output_seahorse_con <- phenotype_ttest(phenotype = pheno, 
                                              phenotypesToCompare = phenoCon, 
                                              phenotypeType = phenoType)
        output_seahorse_padj_con <- stats::p.adjust(output_seahorse_con$pval, method = pval_adj_method)
      }
      
      # Concatenate categorical and continuous results.
      output_seahorse <- data.frame(stat = c(output_seahorse_cat$stat, output_seahorse_con$stat),
                                    pval = c(output_seahorse_cat$pval, output_seahorse_con$pval),
                               cramerV = c(output_seahorse_cat$V, output_seahorse_con$V),
                               padj = c(output_seahorse_padj_cat, output_seahorse_padj_con),
                               testType = c(output_seahorse_cat$testType, output_seahorse_con$testType),
                               row.names = c(rownames(output_seahorse_cat), rownames(output_seahorse_con)))

    }else if(phenoType == "nominal"){
      
      # Do categorical phenotypes.
      whichCat <- which(typesAfter %in% categorical)
      phenoCat <- as.data.frame(phenosAfter[,whichCat])
      colnames(phenoCat) <- colnames(phenosAfter)[whichCat]
      rownames(phenoCat) <- rownames(phenosAfter)
      output_seahorse_cat <- emptyDF
      output_seahorse_padj_cat <- emptyPadj
      if(ncol(phenoCat) > 0){
        output_seahorse_cat <- phenotype_chisq(phenotype = pheno, 
                                                    phenotypesToCompare = phenoCat, 
                                                    phenotypeType = phenoType)
        output_seahorse_padj_cat <- rep(NA, nrow(output_seahorse_cat))
        whichFFH <- which(output_seahorse_cat$testType == "FFH")
        whichChisq <- which(output_seahorse_cat$testType == "Chi-square")
        if(length(whichFFH) > 0){
          output_seahorse_padj_cat[whichFFH] <- stats::p.adjust(output_seahorse_cat[whichFFH, "pval"], method = pval_adj_method)
        }
        if(length(whichChisq) > 0){
          output_seahorse_padj_cat[whichChisq] <- stats::p.adjust(output_seahorse_cat[whichChisq, "pval"], method = pval_adj_method)
        }
      }
      
      # Do continuous phenotypes.
      whichContinuous <- which(typesAfter == "continuous")
      phenoCon <- as.data.frame(phenosAfter[,whichContinuous])
      colnames(phenoCon) <- colnames(phenosAfter)[whichContinuous]
      rownames(phenoCon) <- rownames(phenosAfter)
      output_seahorse_con <- emptyDF
      output_seahorse_padj_con <- emptyPadj
      if(ncol(phenoCon) > 0){
        output_seahorse_con <- phenotype_anova(phenotype = pheno, 
                                              phenotypesToCompare = phenoCon, 
                                              phenotypeType = phenoType)
        output_seahorse_padj_con <- stats::p.adjust(output_seahorse_con$p, method = pval_adj_method)
      }
      
      # Concatenate categorical and continuous results.
      output_seahorse <- data.frame(stat = c(output_seahorse_cat$stat, output_seahorse_con$stat),
                                    pval = c(output_seahorse_cat$p, output_seahorse_con$p),
                                   cramerV = c(output_seahorse_cat$V, output_seahorse_con$V),
                                   padj = c(output_seahorse_padj_cat, output_seahorse_padj_con),
                                   testType = c(output_seahorse_cat$testType, output_seahorse_con$testType),
                                   row.names = c(rownames(output_seahorse_cat), rownames(output_seahorse_con)))
      
    }else{

      # Do dichotomous phenotypes.
      whichDichotomous <- which(typesAfter == "dichotomous")
      phenoDich <- as.data.frame(phenosAfter[,whichDichotomous])
      colnames(phenoDich) <- colnames(phenosAfter)[whichDichotomous]
      rownames(phenoDich) <- rownames(phenosAfter)
      output_seahorse_dich <- emptyDF
      output_seahorse_padj_dich <- emptyPadj
      if(ncol(phenoDich) > 0){
        output_seahorse_dich <- phenotype_ttest(phenotype = pheno, 
                                                    phenotypesToCompare = phenoDich, 
                                                    phenotypeType = phenoType)
        output_seahorse_padj_dich <- stats::p.adjust(output_seahorse_dich$p, method = pval_adj_method)
      }
      
      # Do nominal phenotypes.
      whichNominal <- which(typesAfter == "nominal")
      phenoNom <- as.data.frame(phenosAfter[,whichNominal])
      colnames(phenoNom) <- colnames(phenosAfter)[whichNominal]
      rownames(phenoNom) <- rownames(phenosAfter)
      output_seahorse_nom <- emptyDF
      output_seahorse_padj_nom <- emptyPadj
      if(ncol(phenoNom) > 0){
        output_seahorse_nom <- phenotype_anova(phenotype = pheno, 
                                              phenotypesToCompare = phenoNom, 
                                              phenotypeType = phenoType)
        output_seahorse_padj_nom <- stats::p.adjust(output_seahorse_nom$p, method = pval_adj_method)
      }

      # Do continuous phenotypes.
      whichContinuous <- which(typesAfter == "continuous")
      phenoCon <- as.data.frame(phenosAfter[,whichContinuous])
      colnames(phenoCon) <- colnames(phenosAfter)[whichContinuous]
      rownames(phenoCon) <- rownames(phenosAfter)
      output_seahorse_con <- emptyDF
      output_seahorse_padj_con <- emptyPadj
      if(ncol(phenoCon) > 0){
        output_seahorse_con = phenotype_cor(phenotype = pheno, 
                                              phenotypesToCompare = phenoCon, 
                                            method = method)
        output_seahorse_padj_con <- rep(NA, nrow(output_seahorse_con))
      }
      
      # Concatenate categorical and continuous results.
      output_seahorse <- data.frame(stat = c(output_seahorse_dich$stat, output_seahorse_nom$stat, output_seahorse_con$stat),
                                    pval = c(output_seahorse_dich$p, output_seahorse_nom$p, output_seahorse_con$p),
                                   cramerV = c(output_seahorse_dich$V, output_seahorse_nom$V, output_seahorse_con$V),
                                   padj = c(output_seahorse_padj_dich, output_seahorse_padj_nom, output_seahorse_padj_con),
                                   testType = c(output_seahorse_dich$testType, output_seahorse_nom$testType, output_seahorse_con$testType),
                                   row.names = c(rownames(output_seahorse_dich), rownames(output_seahorse_nom), rownames(output_seahorse_con)))
      
    }
    # Add the NA values.
    varsNotAdded <- setdiff(colnames(phenotype), row.names(output_seahorse))
    additionalVars <- data.frame(stat = rep(NA, length(varsNotAdded)),
                                 pval = rep(NA, length(varsNotAdded)),
                                 cramerV = rep(NA, length(varsNotAdded)),
                                 padj = rep(NA, length(varsNotAdded)),
                                 testType = rep(NA, length(varsNotAdded)),
                                 row.names = varsNotAdded)

    output_seahorse <- rbind(output_seahorse, additionalVars)
    output_seahorse <- output_seahorse[colnames(phenotype),]
    return(output_seahorse)
  })
  names(phenoAssoc) <- colnames(phenotype)[1:(ncol(phenotype)-1)]

  # Compile each vector list into a matrix. Add an extra NA row at the end.
  statMat <- t(as.matrix(do.call(cbind, list(lapply(1:length(phenoAssoc), function(i){
    df <- data.frame(phenoAssoc[[i]]$stat)
    colnames(df) <- names(phenoAssoc)[i]
    rownames(df) <- names(phenotype)
    return(df)
  }), setNames(data.frame(X=rep(NA, ncol(phenotype))), colnames(phenotype)[ncol(phenotype)])))))
  pvalMat <- t(as.matrix(do.call(cbind, list(lapply(1:length(phenoAssoc), function(i){
    df <- data.frame(phenoAssoc[[i]]$pval)
    colnames(df) <- names(phenoAssoc)[i]
    rownames(df) <- names(phenotype)
    return(df)
  }), setNames(data.frame(X=rep(NA, ncol(phenotype))), colnames(phenotype)[ncol(phenotype)])))))
  vMat <- t(as.matrix(do.call(cbind, list(lapply(1:length(phenoAssoc), function(i){
    df <- data.frame(phenoAssoc[[i]]$cramerV)
    colnames(df) <- names(phenoAssoc)[i]
    rownames(df) <- names(phenotype)
    return(df)
  }), setNames(data.frame(X=rep(NA, ncol(phenotype))), colnames(phenotype)[ncol(phenotype)])))))
  adjMat <- t(as.matrix(do.call(cbind, list(lapply(1:length(phenoAssoc), function(i){
    df <- data.frame(phenoAssoc[[i]]$padj)
    colnames(df) <- names(phenoAssoc)[i]
    rownames(df) <- names(phenotype)
    return(df)
  }), setNames(data.frame(X=rep(NA, ncol(phenotype))), colnames(phenotype)[ncol(phenotype)])))))
  typeMat <- t(as.matrix(do.call(cbind, list(lapply(1:length(phenoAssoc), function(i){
    df <- data.frame(phenoAssoc[[i]]$testType)
    colnames(df) <- names(phenoAssoc)[i]
    rownames(df) <- names(phenotype)
    return(df)
  }), setNames(data.frame(X=rep(NA, ncol(phenotype))), colnames(phenotype)[ncol(phenotype)])))))

  # Return all matrices.
  phenoAssocList <- list(stat = statMat, pval = pvalMat, cramerV = vMat, padj = adjMat, testType = typeMat)
  return(phenoAssocList)
}

#' Function to run associations between categorical (dichotomous or nominal) phenotypes.
#' If at least 5 samples exist in each group, run a Chi-square test and compute Cramer's V.
#' Otherwise, run a Fisher-Freeman-Halton Test.
#' @param phenotype The phenotype vector to test.
#' @param phenotypesToCompare A data frame containing all phenotypes against which to run associations.
#' @param phenotypeType Whether the phenotype is dichotomous or nominal.
#' phenotypes against which to compare.
#' @returns A data frame with two vectors: "cor", which lists the Chi-square or Fisher-Freeman-Halton Test p-values, 
#' and "V" which lists the Cramer's V statistics (will be NA if Fisher-Freeman-Halton Test is computed)
phenotype_chisq <- function(phenotype, phenotypesToCompare, phenotypeType){
  
  # Generate the tables for all phenotype pairs.
  allPhenoPairTables <- lapply(phenotypesToCompare, function(phen){
    grpCounts <- table(phenotype, phen)
  })
  names(allPhenoPairTables) <- colnames(phenotypesToCompare)
  
  # Do a chi-square test everywhere else.
  # Run the comparisons one at a time because Chi-square and FFH cannot be scaled.
  chisqRes <- lapply(allPhenoPairTables, function(chisqTable) {
    
    # Check that we have 2 or more non-zero marginals. Only run the test if we do.
    rowMarginals <- rowSums(chisqTable)
    colMarginals <- colSums(chisqTable) 
    nonzeroRowMarginalCount <- length(which(rowMarginals > 0))
    nonzeroColMarginalCount <- length(which(colMarginals > 0))
    pval <- NA
    V <- NA
    stat <- NA
    type <- "Chi-square"
    if(nonzeroRowMarginalCount >= 2 && nonzeroColMarginalCount >= 2){
      # Run the test. If a warning is thrown, switch to FFH.
      chisq <- withCallingHandlers(
        chisq.test(chisqTable),
        warning = function(w) {
          if (grepl("Chi-squared approximation may be incorrect", w$message)) {
            # Do not change only within the warning scope but also within the parent scope.
            type <<- "FFH"
            invokeRestart("muffleWarning")
          }
        }
      )
      
      # Get the p-value and Cramer's V.
      stat <- chisq$statistic
      pval <- chisq$p.value
      V <- sqrt(chisq$statistic / (sum(chisqTable) * min(nrow(chisqTable) - 1, 
                                                         ncol(chisqTable) - 1)))
    }else{
      warning("Less than 2 nonzero marginals in the contingency table - Chi-square will return NA")
    }
    return(list(pval = pval, stat = stat, v = V, testType = type))
  })
  chisqP <- unlist(lapply(chisqRes, function(res){
    return(res$pval)
  }))
  chisqStat <- unlist(lapply(chisqRes, function(res){
    return(res$stat)
  }))
  cramerV <- unlist(lapply(chisqRes, function(res){
    return(res$v)
  }))
  testType <- unlist(lapply(chisqRes, function(res){
    return(res$testType)
  }))
  corResNotLessFull <- data.frame(stat = chisqStat,
                                  pval = chisqP,
                              V = cramerV,
                              testType = testType,
                              row.names = colnames(phenotypesToCompare))
  
  # Run FFS where warning was issued.
  corResNotLess <- corResNotLessFull[which(corResNotLessFull$testType == "Chi-square"),]
  switchedTables <- allPhenoPairTables[which(corResNotLessFull$testType == "FFH")]
  corResSwitched <- phenotype_ffs(switchedTables)

  # Bind together the results from the FFH and Chi-Square tests.
  corRes <- rbind(corResNotLess, corResSwitched)
  corRes <- corRes[colnames(phenotypesToCompare),]

  return(corRes)
}

#' Function to run associations between categorical (dichotomous or nominal) phenotypes using a Fisher-Freeman-Halton Test.
#' @param tables The list of tables to test, named by the phenotypes against which to compare.
#' @returns A data frame with two vectors: "cor", which lists Fisher-Freeman-Halton Test p-values, 
#' and "V" which lists the Cramer's V statistics (NA)
phenotype_ffs <- function(tables){
  
  # Run the comparisons one at a time because FFH cannot be scaled.
  pvals <- unlist(lapply(tables, function(fishTable) {
    fishP <- fisher.test(fishTable, simulate.p.value = TRUE, B = 10000)$p.value
    return(fishP)
  }))

  # Return the results.
  return(data.frame(stat = rep(NA, length(tables)),
                    pval = pvals,
                    V = rep(NA, length(tables)),
                    testType = rep("FFH", length(tables)),
                    row.names = names(tables)))
}

#' Function to compute ANOVA statistics between nominal and continuous phenotypes.
#' @param phenotype The phenotype vector to test.
#' @param phenotypesToCompare A data frame containing all phenotypes against which to run associations.
#' @param phenotypeType Whether the phenotype is nominal or continuous.
#' @returns A data frame with two vectors: "cor", which lists the ANOVA p-values, 
#' and "V" which is set to NA.
phenotype_anova <- function(phenotype, phenotypesToCompare, phenotypeType){
  
  res <- data.frame(p = rep(NA, ncol(phenotypesToCompare)),
                    stat = rep(NA, ncol(phenotypesToCompare)))

  # Case 1 - the phenotype is nominal and the phenotypes to compare are numeric.
  # Case 2 - the phenotype is numeric and the phenotypes to compare are nominal.
  if(phenotypeType == "nominal"){
    phenotype_vector = factor(as.character(phenotype))
    res <- do.call(rbind, lapply(colnames(phenotypesToCompare), function(x){
      results <- data.frame(p = NA, stat = NA)
      tryCatch({
        anovares <- anova(lm(as.numeric(as.character(phenotypesToCompare[,x]))~phenotype_vector))
        results <- data.frame(p = anovares$`Pr(>F)`[1],
                              stat = anovares$`F value`[1])
      }, error = function(cond){
        warning("In this phenotype pair, all continuous values are missing for all but one phenotype level - NA result will be returned")
      })
      return(results)
    }))
  }else{
    res <- do.call(rbind, lapply(1:ncol(phenotypesToCompare), function(i){
      results <- data.frame(p = NA, stat = NA)
      phenotype_vector <- factor(as.character(phenotypesToCompare[,i]))
      tryCatch({
        anovares <- anova(lm(as.numeric(as.character(phenotype))~phenotype_vector))
        results <- data.frame(p = anovares$`Pr(>F)`[1],
                              stat = anovares$`F value`[1])
      }, error = function(cond){
        print(cond)
        warning("In this phenotype pair, all continuous values are missing for all but one phenotype level - NA result will be returned")
      })
      return(results)
    }))
  }
  return(data.frame(stat = res$stat,
                    p = res$p,
                    V = rep(NA, ncol(phenotypesToCompare)),
                    testType = rep("ANOVA", ncol(phenotypesToCompare)),
                    row.names = colnames(phenotypesToCompare)))
}

#' Function to compute Welch's t-test statistics between dichotomous and continuous phenotypes.
#' @param phenotype The phenotype vector to test.
#' @param phenotypesToCompare A data frame containing all phenotypes against which to run associations.
#' @param phenotypeType Whether the phenotype is dichotomous or continuous.
#' @returns A data frame with two vectors: "cor", which lists the t-test p-values, 
#' and "V" which is set to NA.
phenotype_ttest <- function(phenotype, phenotypesToCompare, phenotypeType){
  
  # Case 1 - the phenotype is dichotomous and the phenotypes to compare are numeric.
  # We can vectorize this and use matrixTests.
  # Case 2 - the phenotype is numeric and the phenotypes to compare are dichotomous.
  # We cannot vectorize this and use t.test instead, which defaults to Welch's t-test.
  res <- data.frame(p = rep(NA, ncol(phenotypesToCompare)),
                    stat = rep(NA, ncol(phenotypesToCompare)))
  if(phenotypeType == "dichotomous"){
    
    # Split groups.
    phenotype_vector = factor(as.character(phenotype))
    levels <- unique(phenotype_vector)
    group1 <- phenotypesToCompare[which(phenotype_vector == levels[1]),]
    group2 <- phenotypesToCompare[which(phenotype_vector == levels[2]),]
    
    # Check that thresholds are met.
    whichLevel1 <- length(which(phenotype_vector == levels[1]))
    whichLevel2 <- length(which(phenotype_vector == levels[2]))
    
    if(length(dim(group1)) >= 2 || length(dim(group2)) >= 2){
      meetsThreshold <- colSums(!is.na(group1)) >= 2 & colSums(!is.na(group2)) >= 2
      
      # Initialize result to NA.
      res <- data.frame(p = rep(NA, ncol(group1)),
                         stat = rep(NA, ncol(group1)))
      rownames(res) <- colnames(group1)
      
      # Only compute correlations if both levels are represented in the phenotype vector.
      if(whichLevel1 > 0 && whichLevel2 > 0 && length(which(meetsThreshold == TRUE)) > 0){
        # Compute t-test where thresholds are met.
        tresValid <- matrixTests::col_t_welch(
          group1[, meetsThreshold, drop = FALSE],
          group2[, meetsThreshold, drop = FALSE]
        )
        
        # Compute remaining t-tests.
        res[meetsThreshold, "p"] <- tresValid$pvalue
        res[meetsThreshold, "stat"] <- tresValid$statistic
      }
    }else{
      meetsThreshold <- sum(!is.na(group1)) >= 2 & sum(!is.na(group2)) >= 2
      
      # Initialize result to NA.
      rownames(res) <- names(group1)
      
      # Only compute correlations if both levels are represented in the phenotype vector.
      if(whichLevel1 > 0 && whichLevel2 > 0 && length(which(meetsThreshold == TRUE)) > 0){
        tryCatch({
          # Compute t-test where thresholds are met.
          tresValid <- t.test(group1[meetsThreshold, drop = FALSE], 
                              group2[meetsThreshold, drop = FALSE])

          # Compute remaining t-tests.
          res[meetsThreshold, "p"] <- tresValid$p.value
          res[meetsThreshold, "stat"] <- tresValid$statistic
        }, error = function(cond){
          warning("Could not compute t-test (it is possible that variance is too low). Returning NA.")
        })
      }
    }
  }else{
    res <- do.call(rbind, lapply(1:ncol(phenotypesToCompare), function(i){
      phenotype_vector = factor(as.character(phenotypesToCompare[,i]))
      levels <- unique(phenotype_vector)
      group1 <- phenotype[which(phenotype_vector == levels[1])]
      group2 <- phenotype[which(phenotype_vector == levels[2])]
      stat <- NA
      p <- NA
      if(length(which(!is.na(group1))) > 2 && length(which(!is.na(group2))) > 2){
        tryCatch({
          tres <- t.test(group1, group2)
          p = tres$p.value
          stat = tres$statistic
        }, error = function(cond){
          warning("Could not compute t-test (it is possible that variance is too low). Returning NA.")
        })
      }
      toreturn <- data.frame(p = p, stat = stat)
      rownames(toreturn) <- colnames(phenotypesToCompare)[i]
      return(toreturn)
    }))
  }
  if(length(which(is.na(res))) > 0){
    warning("Some phenotypes did not have sufficient sample sizes to perform a t-test - NAs will be returned")
  }
  
  # Return the data frame.
  return(data.frame(stat = res$stat,
                    pval = res$p,
                    V = rep(NA, ncol(phenotypesToCompare)),
                    testType = rep("T-Test", ncol(phenotypesToCompare)),
                    row.names = colnames(phenotypesToCompare)))
  
}

#' Function to compute correlation between continuous phenotypes.
#' @param phenotype The phenotype vector to test.
#' @param phenotypesToCompare A data frame containing all phenotypes against which to run associations.
#' phenotypes against which to compare.
#' @param method : One of "pearson", "spearman", or "kendall".
#' @returns A data frame with two vectors: "cor", which lists the t-test p-values, 
#' and "V" which is set to NA.
phenotype_cor <- function(phenotype, phenotypesToCompare, method){
  corRes <- psych::corr.test(phenotype, phenotypesToCompare, use = "pairwise", method = method)
  corvals <- as.numeric(corRes$r)
  pvals <- as.numeric(corRes$p)
  names(corvals) <- colnames(phenotype)
  names(pvals) <- colnames(phenotype)
  return(data.frame(stat = corvals,
                    p = pvals,
                    testType = rep("Cor", ncol(phenotypesToCompare)),
                    V = rep(NA, ncol(phenotypesToCompare)),
                    row.names = colnames(phenotypesToCompare)))
}

#' Generates the following files which are needed for use with the SEAHORSE UI:
#' - all_gsea_results.tsv.gz: A list of pathway enrichment scores with respect to all
#' phenotypes and all files, formatted as pathway	pval	padj	log2err	ES	NES	size	ranks	leadingEdge	varname	tissue
#' - data_dictionary.tsv.gz: A list of phenotype names, descriptions, and variable types formatted as
#' VARNAME VARDESC VARMETA TYPE
#' - geneexpression2geneexpression.tsv.gz: Gene-gene association statistics, formatted as Gene A	Gene B	Tissue	Correlation
#' - geneexpression_data.tsv.gz: Gene expression data, formatted as ENSG    SAMPID  GENE_EXPRESSION
#' - human_ensembl2symbol_map.tsv.gz: A mapping from Ensembl IDs to HGNC symbols and Entrez IDs,
#' formatted as ALIAS	ENSEMBL	SYMBOL	ENTREZID
#' - metadata.tsv.gz: The phenotype data, formatted as SAMPID	tissue	VARNAME	VALUE
#' - metadata2expression.tsv.gz: Gene-phenotype association statistics, formatted as 
#' VARNAME GENE    tissue  TEST    TESTSTAT        TESTPVALUE
#' - metadata2metadata.tsv.gz: Phenotype-phenotype association statistics, formatted as 
#' VARNAME1        VARNAME2        tissue  TEST    TESTSTAT        TESTPVALUE
#' @param input_directory : A directory containing all inputs to SEAHORSE, stored as RDS files,
#' where the input values contain slots for expression, phenotype, and phenotype dictionary. 
#' @param result_directory : A directory containing all SEAHORSE results from the expression data,
#' phenotype data, phenotype dictionary, and pathways provided, stored as RDS files. Files must follow
#' the same naming convention as in the result directory.
#' @param output_directory : The directory where you wish to store the files generated by this function.
#' @param data_dictionary : The data dictionary.
#' @param pathways : The pathways used to run SEAHORSE
#' @param geneGeneCutoff : The correlation cutoff for gene-gene associations. Only correlations greater
#' than this cutoff will be saved.
#' @return NULL
#' @export
seahorseFormatForUI <- function(input_directory, result_directory, output_directory,
                                data_dictionary, pathways, geneGeneCutoff = -2){
  
  # Check packages.
  if (!requireNamespace("AnnotationDbi", quietly = TRUE)) {
    stop("Package 'AnnotationDbi' is required but not installed.")
  }
  if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    stop("Package 'org.Hs.eg.db' is required but not installed.")
  }
  
  # Create output directory if it doesn't already exist.
  if(!dir.exists(output_directory)){
    dir.create(output_directory)
  }
  
  # Check format of data dictionary.
  if(is(data_dictionary, "data.frame")){
    if(!identical(colnames(data_dictionary), c("VARNAME", "VARDESC", "VARMETA", "TYPE"))){
      stop("Incorrect format for data dictionary")
    }
  }else{
    stop("Incorrect format for data dictionary")
  }
  
  # Loop through all files.
  resultfiles <- list.files(result_directory)
  inputfiles <- list.files(input_directory)
  tissue_names_result <- formatTissueNames(resultfiles)
  tissue_names_input <- formatTissueNames(inputfiles)
  
  # Check that names match.
  if(length(tissue_names_result) == length(tissue_names_input)){
    if(!identical(tissue_names_result, tissue_names_input)){
      stop("Result and input file names do not match")
    }
  }else{
    stop("Number of result and input files do not match")
  }
  
  # If the output directory does not exist, create it.
  if(!dir.exists(output_directory)){
    dir.create(output_directory)
  }
  
  # Write results.
  message("Preparing to write expression...")
  writeExpression(inFiles = paste(input_directory, inputfiles, sep = "/"),
                  file = paste(output_directory, "geneexpression_data", sep = "/"))
  message("Expression file done.")
  message("Preparing to write gene mapping...")
  writeMapping(inFiles = paste(input_directory, inputfiles, sep = "/"),
               file = paste(output_directory, "human_ensembl2symbol_map", sep = "/"))
  message("Mapping file done.")
  message("Preparing to write data dictionary...")
  writeToGz(data = data_dictionary, file = paste(output_directory, "data_dictionary", sep = "/"))
  message("Data dictionary done.")
  message("Preparing to write GSEA results...")
  writePathwayEnrichment(inFiles = paste(result_directory, resultfiles, sep = "/"),
                         file = paste(output_directory, "all_gsea_results", sep = "/"),
                         pathways = pathways)
  message("GSEA result file done.")
  message("Preparing to write phenotype data...")
  writePhenotype(inFiles = paste(input_directory, inputfiles, sep = "/"),
                         file = paste(output_directory, "metadata", sep = "/"))
  message("Phenotype data done.")
  message("Preparing to write phenotype-gene associations...")
  writePhenotypeGeneAssociations(inFiles = paste(input_directory, inputfiles, sep = "/"),
                  resultFiles = paste(result_directory, resultfiles, sep = "/"),
                 file = paste(output_directory, "metadata2expression", sep = "/"))
  message("Phenotype-gene associations done.")
  message("Preparing to write phenotype associations...")
  writePhenotypePhenotypeAssociations(inFiles = paste(input_directory, inputfiles, sep = "/"),
                                      resultFiles = paste(result_directory, resultfiles, sep = "/"),
                                 file = paste(output_directory, "metadata2metadata", sep = "/"))
  message("Phenotype associations done")
  message("Preparing to write gene associations...")
  writeGeneGene(inFiles = paste(result_directory, resultfiles, sep = "/"),
                                      file = paste(output_directory, "geneexpression2geneexpression", sep = "/"),
                cutoff = geneGeneCutoff)
  message("Gene associations done.")
  
}

#' Generates the pathway enrichment scores across all tissues and writes them.
#' Format is pathway	pval	padj	log2err	ES	NES	size	ranks	leadingEdge	varname	tissue
#' @param inFiles The SEAHORSE files
#' @param file The file where results should be stored. .tsv.gz will be appended.
#' @param pathways The pathways used to run SEAHORSE
#' @param cutoff The adjusted p-value cutoff for pathway enrichment. Default is 2 (all retained).
#' @return NULL
writePathwayEnrichment <- function(inFiles, file, pathways, cutoff = 2){
  
  # Get input file names.
  tissueNames <- formatTissueNames(inFiles)
  
  # Open the file.
  con <- gzfile(paste0(file, ".tsv.gz"), "wt")
  on.exit(close(con), add = TRUE)
  
  # Write the pathway enrichment results.
  pathwayEnrichment <- tryCatch({for(i in 1:length(inFiles)){
    
    # Get tissue name and read file.
    results = readRDS(inFiles[i])

    enrichmentVar <- data.frame(pathway = character(), pval = numeric(), 
                                padj = numeric(), log2err = numeric(),
                                ES = numeric(), NES = numeric(),
                                size = numeric(), ranks = character(),
                                leadingEdge = character(), varname = character(),
                                tissue = character())
    
    if(length(results$GSEA) > 0){
      # Loop through variable names.
      enrichmentVar <- do.call(rbind, lapply(names(results$GSEA), function(var){
        
        # Get full ranking of genes for this var based on whether correlation or p-value was used.
        # P-values are log-transformed and absolute values are used for correlations.
        # The largest values should then receive the lowest (best) ranking.
        transformedStats <- results$phenotype_association[[var]]$stat
        if(!all(is.null(transformedStats))){
          transformedStats[which(!is.na(results$phenotype_association[[var]]$padj))] <- -log10(transformedStats[which(!is.na(results$phenotype_association[[var]]$padj))])
          transformedStats[which(is.na(results$phenotype_association[[var]]$padj))] <- abs(transformedStats[which(is.na(results$phenotype_association[[var]]$padj))])
          pvalRanking <- order(-transformedStats)
          names(pvalRanking) <- rownames(results$phenotype_association[[var]])
        }
        
        # Build basic pathway data.
        returnDat <- data.frame(pathway = character(), pval = numeric(), 
                                padj = numeric(), log2err = numeric(),
                                ES = numeric(), NES = numeric(),
                                size = numeric(), ranks = character(),
                                leadingEdge = character(), varname = character(),
                                tissue = character())
        if(!all(is.na(results$GSEA[[var]]))){
          returnDat <- results$GSEA[[var]][,c("pathway", "pval", "padj", "log2err", "ES", "NES", "size")]
          # Get ranking of pathway genes.
          varResults <- as.data.frame(results$GSEA[[var]])
          rankVec <- unlist(lapply(1:nrow(varResults), function(pw){
            # Get pathway genes.
            pathwayGenes <- pathways[[varResults[pw, "pathway"]]]
            
            # Get adjusted p-value rankings.
            rankForLeadingEdge <- pvalRanking[pathwayGenes]
            return(paste0("{", paste(rankForLeadingEdge, collapse = ","), "}"))
          }))
          returnDat$ranks <- rankVec
          # Format leading edges as string.
          leadingEdgeVec <- unlist(lapply(results$GSEA[[var]][,"leadingEdge"][[1]], function(edge){
            return(paste0("{", paste(edge, collapse = ","), "}"))
          }))
          returnDat$leadingEdge <- leadingEdgeVec
          
          # Add variable name.
          returnDat$varname <- rep(var, nrow(returnDat))
          
          # Add tissue.
          splitPath <- strsplit(inFiles[i], "/")[[1]]
          localFile <- splitPath[length(splitPath)]
          withoutFileExt <- strsplit(localFile, ".RDS")[[1]][1]
          returnDat$tissue <- rep(withoutFileExt, nrow(returnDat))
          
          # Subset by cutoff.
          returnDat <- returnDat[which(returnDat$padj < cutoff),]
        }
        return(returnDat)
      }))
    }
    
    # Write the chunk.
    write.table(enrichmentVar, file = con, sep = "\t", row.names = FALSE, quote = FALSE,
                col.names = i == 1)
    message(paste("Wrote data for", i, "out of", length(inFiles), "tissues"))
    
  }}, error = function(cond){
    print(cond)
    stop(paste("Could not format pathway enrichment"))})
}

#' Formats the tissue names by removing the .RDS extension from the input files.
#' @param names The list of tissue names.
#' @return A list of names.
formatTissueNames <- function(names){
  return(unlist(lapply(names, function(name){
    nameSplit <- strsplit(name, "/")[[1]]
    localName <- nameSplit[length(nameSplit)]
    result <- strsplit(localName, ".RDS")[[1]][1]
    return(result)
  })))
}

#' Generates the gene correlation results in the format expected by the UI.
#' Format is Gene A	Gene B	Tissue	Correlation
#' @param inFiles The SEAHORSE files
#' @param file The file where results should be stored. .tsv.gz will be appended.
#' @param removeLowerTri Whether or not to remove the lower triangular matrix
#' @param cutoff The correlation cutoff for filtering. Default is -2 (all retained).
#' (default is TRUE)
#' @return NULL 
writeGeneGene <- function(inFiles, file, removeLowerTri = TRUE, cutoff = -2){
  
  # Open the file.
  con <- gzfile(paste0(file, ".tsv.gz"), "wt")

  # Get input file names.
  tissueNames <- formatTissueNames(inFiles)

  # Write the gene correlation results.
  geneCor <- for(i in 1:length(tissueNames)){
    
    # Get tissue name and read file.
    results = readRDS(inFiles[i])

    # Remove duplicate pairs.
    geneCor <- results$coexpression
    flatTable <- data.frame(A = character(), B = character(), 
                            Tissue = character(), Correlation = numeric())
    colnames(flatTable) <- c("Gene A", "Gene B", "Tissue", "Correlation")
    if(length(which(!is.na(geneCor) == TRUE)) > 0){
      if(removeLowerTri == TRUE){
        geneCor[lower.tri(geneCor, diag = TRUE)] <- NA
      }
      flatTable <- data.frame(
        A = rep(colnames(geneCor), times = nrow(geneCor)),
        B = rep(rownames(geneCor), each = ncol(geneCor)),
        Correlation = as.vector(geneCor)
      )

      # Remove NA values. To improve memory, remove unneeded variables and run
      # garbage collection.
      rm(results)
      rm(geneCor)
      gc()
      keep <- complete.cases(flatTable)
      newFlatTable <- flatTable[keep, , drop = FALSE]
      rm(flatTable, keep)
      gc()
      flatTable <- newFlatTable
      rm(newFlatTable)
      gc()
      
      # Add names of variables.
      colnames(flatTable) <- c("Gene A", "Gene B", "Correlation")

      # Add tissue.
      flatTable$Tissue <- rep(tissueNames[i], nrow(flatTable))

      # Rearrange columns.
      flatTable <- flatTable[,c("Gene A", "Gene B", "Tissue", "Correlation")]

      # Subset by cutoff.
      flatTable <- flatTable[which(abs(flatTable$Correlation) > cutoff),]
    }
    
    # Write the chunk.
    write.table(flatTable, file = con, sep = "\t", row.names = FALSE, quote = FALSE,
                col.names = i == 1)
    message(paste("Wrote data for", i, "out of", length(inFiles), "tissues"))
  }
  
  # Close the file.
  close(con)
}

#' Generates the expression data in the format expected by the UI. Writes a
#' data frame formatted as ENSG    SAMPID  GENE_EXPRESSION as .tsv.gz.
#' @param inFiles The files to input to SEAHORSE, formatted as RDS.
#' @param file The name of the file to write. .tsv.gz. will be appended.
#' @return NULL
writeExpression <- function(inFiles, file){
  
  # Open the file.
  con <- gzfile(paste0(file, ".tsv.gz"), "wt")
  on.exit(close(con))
  
  # Get input file names.
  tissueNames <- formatTissueNames(inFiles)

  # Write the pathway enrichment results.
  expression <- tryCatch({for(i in 1:length(inFiles)){
    
    # Get tissue name and read file.
    results = readRDS(inFiles[i])

    # Flatten it.
    exprMat <- as.matrix(results$expression)
    flatTable <- data.frame(
      ENSG = rep(rownames(exprMat), times = ncol(exprMat)),
      SAMPID = rep(colnames(exprMat), each = nrow(exprMat)),
      GENE_EXPRESSION = as.vector(exprMat),
      tissue = tissueNames[i]
    )
    
    write.table(flatTable, file = con, sep = "\t", row.names = FALSE, quote = FALSE,
                col.names = i == 1)
    message(paste("Wrote data for", i, "out of", length(inFiles), "tissues"))
  }
  }, error = function(cond){
    print(cond)
    stop("Could not format expression")})
}

#' Writes a data frame as .csv.gz.
#' @param data The data frame
#' @param file The name of the file to write. .tsv.gz. will be appended.
#' @return NULL
writeToGz <- function(data, file){
  gz <- gzfile(paste0(file, ".tsv.gz"), "w")
  on.exit(close(gz))
  write.table(data, file = gz, sep = "\t", row.names = FALSE, quote = FALSE)
}

#' Generates the gene mapping data in the format expected by the UI. Writes as
#' .tsv.gz formatted as ALIAS	ENSEMBL	SYMBOL	ENTREZID.
#' @param inFiles The files to input to SEAHORSE, formatted as RDS.
#' @param file The name of the file to write. .tsv.gz. will be appended.
#' @return NULL
writeMapping <- function(inFiles, file){
  
  # Get all unique genes.
  ensembl <- tryCatch({unique(unlist(lapply(1:length(inFiles), function(i){
    
    # Return genes from each file.
    results = readRDS(inFiles[i])
    return(rownames(results$expression))
  })))}, error = function(cond){
    stop("Incorrect format for input file")})

  # Map genes.
  mappedSymbols <- AnnotationDbi::mapIds(
    x = org.Hs.eg.db::org.Hs.eg.db,
    keys = ensembl,
    column = "SYMBOL",
    keytype = "ENSEMBL",
    multiVals = "first" # Returns one value for each ENSEMBL ID.
  )
  mappedEntrez <- AnnotationDbi::mapIds(
    x = org.Hs.eg.db::org.Hs.eg.db,
    keys = ensembl,
    column = "ENTREZID",
    keytype = "ENSEMBL",
    multiVals = "first" # Returns one value for each ENSEMBL ID.
  )
  mappedGenes <- data.frame(ENSEMBL = ensembl, SYMBOL = mappedSymbols, ENTREZID = mappedEntrez)
  mappedGenes$ALIAS <- mappedGenes$SYMBOL
  mappedGenes <- mappedGenes[,c("ALIAS", "ENSEMBL", "SYMBOL", "ENTREZID")]

  # Write the file.
  writeToGz(data = mappedGenes, file = file)
}

#' Generates the phenotype data in the format expected by the UI.
#' Format is SAMPID	tissue	VARNAME	VALUE
#' @param inFiles The files to input to SEAHORSE, formatted as RDS.
#' @param file The file where results should be stored. .tsv.gz will be appended.
#' @return NULL
writePhenotype <- function(inFiles, file){
  
  # Open the file.
  con <- gzfile(paste0(file, ".tsv.gz"), "wt")
  on.exit(close(con))

  # Get input file names.
  tissueNames <- formatTissueNames(inFiles)

  # Write the pathway enrichment results.
  pheno <- tryCatch({for(i in 1:length(tissueNames)){
    
    # Get tissue name and read file.
    results = readRDS(inFiles[i])

    # Convert to a table.
    phenoMat <- as.matrix(results$phenotype)
    flatTable <- data.frame(
      VARNAME = rep(colnames(phenoMat), each = nrow(phenoMat)),
      SAMPID = rep(rownames(phenoMat), times = ncol(phenoMat)),
      VALUE = as.vector(phenoMat)
    )

    # Add tissue
    flatTable$tissue <- rep(tissueNames[i], nrow(flatTable))

    # Rearrange columns.
    flatTable <- flatTable[,c("SAMPID", "tissue", "VARNAME", "VALUE")]
    
    write.table(flatTable, file = con, sep = "\t", row.names = FALSE, quote = FALSE,
                col.names = i == 1)
    message(paste("Wrote data for", i, "out of", length(inFiles), "tissues"))

  }}, error = function(cond){
    stop("Incorrect format for input file")})
}

#' Generates the phenotype data in the format expected by the UI.
#' Format is VARNAME GENE    tissue  TEST    TESTSTAT        TESTPVALUE
#' @param inFiles The SEAHORSE inputs
#' @param resultFiles The SEAHORSE results
#' @param file The file where results should be stored. .tsv.gz will be appended.
#' @param corCutoff The cutoff for correlation. Default is -2 (all values will be retained.)
#' @param padjCutoff The cutoff for adjusted p-value. Default is 2 (all values will be retained.)
#' @return NULL
writePhenotypeGeneAssociations <- function(inFiles, resultFiles, file, corCutoff = -2, padjCutoff = 2){
  
  # Open the file.
  con <- gzfile(paste0(file, ".tsv.gz"), "wt")
  on.exit(close(con))

  # Get input file names.
  tissueNames <- formatTissueNames(inFiles)
  
  # Write the pathway enrichment results.
  phenotypeGene <- for(i in 1:length(tissueNames)){
    
    # Get tissue name and read file.
    results <- readRDS(resultFiles[i])
    inputs <- readRDS(inFiles[i])
    
    phenotypeToGeneDf <- data.frame(VARNAME = character(), GENE = character(), 
                                    tissue = character(), TEST = character(),
                                    TESTSTAT = numeric(), TESTPVALUE = numeric())
    if(length(results$phenotype_association) > 0){
        phenotypeToGeneList <- lapply(names(results$phenotype_association), function(var){
          outDf <- data.frame(VARNAME = character(), GENE = character(), 
                              tissue = character(), TEST = character(),
                              TESTSTAT = numeric(), TESTPVALUE = numeric())
          if(var %in% colnames(inputs$phenotype)){
            outDf <- results$phenotype_association[[var]][,c("stat", "padj")]
            if(length(outDf) > 1){
              colnames(outDf) <- c("TESTSTAT", "TESTPVALUE")
              outDf$GENE <- rownames(outDf)
              outDf$VARNAME <- rep(var, nrow(outDf))
              outDf$tissue <- tissueNames[i]
              outDf$TEST <- "Correlation"
              dictVal <- inputs$dict[which(colnames(inputs$phenotype) == var)]
              if(dictVal == "dichotomous"){
                outDf$TEST <- "LIMMA Moderated t-test"
              }else if(dictVal == "nominal"){
                outDf$TEST <- "ANOVA"
              }
              outDf <- outDf[,c("VARNAME", "GENE", "tissue", "TEST", "TESTSTAT", "TESTPVALUE")]
            }
          }
          
          return(outDf)
        })
        phenotypeToGeneDf <- do.call(rbind, phenotypeToGeneList)
        
        # Subset by cutoffs.
        whichCorCutoff <- intersect(which(abs(phenotypeToGeneDf$TESTSTAT) > corCutoff),
                                    which(phenotypeToGeneDf$TEST == "Correlation"))
        whichPadjCutoff <- intersect(which(phenotypeToGeneDf$TESTPVALUE < padjCutoff),
                                    which(phenotypeToGeneDf$TEST != "Correlation"))
        whichPadjCorCutoff <- intersect(which(phenotypeToGeneDf$TESTPVALUE < padjCutoff),
                                     which(phenotypeToGeneDf$TEST == "Correlation"))
        phenotypeToGeneDf <- phenotypeToGeneDf[sort(c(intersect(whichCorCutoff, whichPadjCorCutoff), 
                                                      whichPadjCutoff)),]
    }
    write.table(phenotypeToGeneDf, file = con, sep = "\t", row.names = FALSE, quote = FALSE,
                col.names = i == 1)
    message(paste("Wrote data for", i, "out of", length(inFiles), "tissues"))
  }
}

#' Approximates the p-value given correlation and the number of samples
#' @param rho The correlation vector for a given phenotype
#' @param n The vector of sample counts used to calculate rho
#' @return A vector of p-values
spearmanPFromRho <- function(rho, n) {
  t_stat <- rho * sqrt((n - 2) / (1 - rho^2))
  p_val <- 2 * pt(abs(t_stat), df = n - 2, lower.tail = FALSE)
  return(p_val)
}

#' Generates the phenotype results in the format expected by the UI.
#' Format is VARNAME1        VARNAME2        tissue  TEST    TESTSTAT        TESTPVALUE
#' @param inFiles The SEAHORSE input files
#' @param resultFiles The SEAHORSE results
#' @param file The file where results should be stored. .tsv.gz will be appended.
#' @param corCutoff The cutoff for correlation. Default is -2 (all values will be retained.)
#' @param padjCutoff The cutoff for adjusted p-value. Default is 2 (all values will be retained.)
#' Default is FALSE.
#' @return NULL
writePhenotypePhenotypeAssociations <- function(inFiles, resultFiles, file, corCutoff = -2,
                                                padjCutoff = 2){
  
  # Open the file.
  con <- gzfile(paste0(file, ".tsv.gz"), "wt")

  # Get input file names.
  tissueNames <- formatTissueNames(inFiles)
  
  # Write the pathway enrichment results.
  phenotype <- for(i in 1:length(tissueNames)){
    
    # Get tissue name and read file.
    results = readRDS(resultFiles[i])
    inputs = readRDS(inFiles[i])
    
    phenotypeToPhenotypeDf <- data.frame(VARNAME1 = character(), VARNAME2 = character(), tissue = character(),
                                         TEST = character(), TESTSTAT = numeric(), TESTPVALUE = numeric())
    
    if(all(!is.na(results$phenocor)) == TRUE){
      phenotypeToPhenotypeList <- lapply(1:(nrow(results$phenocor$stat)-1), function(j){
        VARNAME2 <- colnames(results$phenocor$stat)[(j+1):ncol(results$phenocor$stat)]
        TEST <- results$phenocor$testType[j,(j+1):ncol(results$phenocor$stat)]
        TESTSTAT <- results$phenocor$stat[j,(j+1):ncol(results$phenocor$stat)]
        TESTPVALUE <- results$phenocor$padj[j,(j+1):ncol(results$phenocor$padj)]
        var <- rownames(results$phenocor$stat)[j]
        VARNAME1 <- rep(var, length((j+1):ncol(results$phenocor$padj)))
        tissue <- tissueNames[i]
        
        return(data.frame(VARNAME1 = VARNAME1, VARNAME2 = VARNAME2, tissue = tissue,
                          TEST = TEST, TESTSTAT = TESTSTAT, TESTPVALUE = TESTPVALUE))
      })
      phenotypeToPhenotypeDf <- do.call(rbind, phenotypeToPhenotypeList)
      
      # Subset by cutoffs.
      whichCorCutoff <- intersect(which(abs(phenotypeToPhenotypeDf$TESTSTAT) > corCutoff),
                                  which(phenotypeToPhenotypeDf$TEST == "Correlation"))
      whichPadjCutoff <- intersect(which(phenotypeToPhenotypeDf$TESTPVALUE < padjCutoff),
                                   which(phenotypeToPhenotypeDf$TEST != "Correlation"))
      whichPadjCorCutoff <- intersect(which(phenotypeToPhenotypeDf$TESTPVALUE < padjCutoff),
                                      which(phenotypeToPhenotypeDf$TEST == "Correlation"))
      phenotypeToPhenotypeDf <- phenotypeToPhenotypeDf[sort(c(intersect(whichCorCutoff, whichPadjCorCutoff), 
                                                    whichPadjCutoff)),]
    }
    write.table(phenotypeToPhenotypeDf, file = con, sep = "\t", row.names = FALSE, quote = FALSE,
                col.names = i == 1)
    message(paste("Wrote data for", i, "out of", length(inFiles), "tissues"))
  }
  
  # Close the file.
  close(con)
}

#' Generates a rug plot for a single pathway-phenotype association.
#' @param result The SEAHORSE result
#' @param pathwayName The pathway of interest
#' @param phenotypeName The phenotype of interest
#' @param pathways The pathways used as input to SEAHORSE
#' @return NULL
rugPlot <- function(result, pathwayName, pathways, phenotypeName) {

  # Check for formatting errors.
  if(is.null(names(result))){
    stop("Result format is incorrect")
  }else if(!("GSEA" %in% names(result)) ||
     !("phenotype_association" %in% names(result))){
    stop("Result format is incorrect")
  }
  if(is.null(names(pathways))){
    stop("Pathway format is incorrect")
  }
  
  # Extract result.
  fgseaResAll <- result$GSEA
  if(!(phenotypeName %in% names(fgseaResAll))){
    stop(paste("Phenotype", phenotypeName, "does not have a pathway enrichment result"))
  }
  fgseaRes <- fgseaResAll[[phenotypeName]]
  
  # Extract genes.
  if(!(pathwayName %in% names(pathways))){
    stop(paste("Pathway", pathwayName, "does not exist in the list of pathways"))
  }
  pathwayGenes <- pathways[[pathwayName]]
  
  # Extract gene association data.
  result$phenotype_association[[phenotypeName]]$padj <- as.numeric(result$phenotype_association[[phenotypeName]]$padj)
  transformedStats <- result$phenotype_association[[phenotypeName]]$stat
  names(transformedStats) <- rownames(result$phenotype_association[[phenotypeName]])
  
  # Only proceed if gene association data are available.
  if(!all(is.null(transformedStats))){
    
    # Transform to absolute valued correlation or negative log p-value, then sort to obtain
    # ranks.
    transformedStats[which(!is.na(result$phenotype_association[[phenotypeName]]$padj))] <- 
      -log10(pmax(transformedStats[which(!is.na(result$phenotype_association[[phenotypeName]]$padj))],
                  .Machine$double.xmin))
    transformedStats[which(is.na(result$phenotype_association[[phenotypeName]]$padj))] <- 
      abs(transformedStats[which(is.na(result$phenotype_association[[phenotypeName]]$padj))])
    sortedStats <- transformedStats[order(-transformedStats)]
    N <- length(sortedStats)
    
    # Assign rank based on sorting.
    hitPositions <- which(names(sortedStats) %in% pathwayGenes)
    hitPositions <- sort(unique(hitPositions))
    
    par(mar = c(4, 4, 4, 1))
    plot(
      seq_len(N),
      rep(0, N),
      type = "n",
      axes = FALSE,
      xlab = "Gene rank position",
      ylab = "",
      ylim = c(0, 1),
      main = paste("Association Between", pathwayName, "and", phenotypeName)
    )
    
    graphics::segments(
      x0 = hitPositions,
      y0 = 0,
      x1 = hitPositions,
      y1 = 1
    )
    
    axis(1)
    box()
  }else{
    stop(paste("No gene association data for", phenotypeName))
  }
}

#' Generates a summary histogram for a gene.
#' @param expression The gene expression data, formatted as for SEAHORSE
#' @param geneName The name of the gene
#' @param breaks The number of breaks in the histogram
#' @return NULL
summaryHistogramGene <- function(expression, geneName, breaks = 20){
  
  # Formatting check
  if(!is.data.frame(expression) && !is.matrix(expression)){
    stop("Expression data format must be a data frame or matrix")
  }
  if(!(geneName %in% rownames(expression))){
    stop(paste("Gene", geneName, "not found in expression data"))
  }
  if(breaks < 1){
    stop("Number of breaks must be at least 1")
  }
  
  # Plot histogram.
  hist(as.numeric(expression[geneName,]), breaks = breaks, xlab = paste(geneName, "Expression"),
       main = paste("Histogram of", geneName, "Expression"))
}

#' Generates a summary histogram for a phenotype.
#' @param phenotype The phenotype data, formatted as for SEAHORSE
#' @param phenotypeName The name of the phenotype
#' @param phenotype_dictionary The phenotype dictionary, formatted as for SEAHORSE
#' @param breaks The number of breaks in the histogram
#' @return NULL
summaryHistogramPhenotype <- function(phenotype, phenotypeName, phenotype_dictionary, breaks = 20){
  
  # Formatting check
  if(!is.data.frame(phenotype) && !is.matrix(phenotype)){
    stop("Phenotype data format must be a data frame or matrix")
  }
  if(!(phenotypeName %in% colnames(phenotype))){
    stop(paste("Phenotype", phenotypeName, "not found in phenotype data"))
  }
  if(breaks < 1){
    stop("Number of breaks must be at least 1")
  }
  
  # Plot histogram.
  if(phenotype_dictionary[which(colnames(phenotype) == phenotypeName)] == "continuous"){
    hist(as.numeric(phenotype[,phenotypeName]), breaks = breaks, xlab = phenotypeName,
         main = paste("Histogram of", phenotypeName))
  }else{
    graphics::barplot(table(phenotype[,phenotypeName]), xlab = phenotypeName, ylab = "frequency",
            main = paste("Bar plot of", phenotypeName))
  }
  
}

#' Plots the association between two variables
#' @param expression The gene expression data, formatted as for SEAHORSE
#' @param phenotype The phenotype data, formatted as for SEAHORSE
#' @param phenotype_dictionary The phenotype dictionary, formatted as for SEAHORSE
#' @param variableX The variable of interest to plot on the x-axis (a phenotype or gene)
#' @param variableY The variable of interest to plot on the y-axis (a phenotype or gene)
#' @param variableTypeX Either "gene" or "phenotype" (x-axis variable)
#' @param variableTypeY Either "gene" or "phenotype" (y-axis variable)
#' @return NULL
plotAssociation <- function(expression, phenotype, phenotype_dictionary,
                            variableX, variableY, variableTypeX, variableTypeY){
  
  # Check packages.
  if (!requireNamespace("graphics", quietly = TRUE)) {
    stop("Package 'graphics' is required but not installed.")
  }
  
  # Check formatting.
  if(!is.data.frame(expression) && !is.matrix(expression)){
    stop("Expression data are in the wrong format")
  }
  if(!is.data.frame(phenotype) && !is.matrix(phenotype)){
    stop("Phenotype data are in the wrong format")
  }
  if(!(variableTypeX %in% c("gene", "phenotype"))){
    stop(paste("Type", variableTypeX, "is invalid"))
  }
  if(!(variableTypeY %in% c("gene", "phenotype"))){
    stop(paste("Type", variableTypeY, "is invalid"))
  }
  if(variableTypeX == "gene"){
    if(!variableX %in% rownames(expression)){
      stop(paste("Gene", variableX, "is not in the expression data"))
    }
  }
  if(variableTypeY == "gene"){
    if(!variableY %in% rownames(expression)){
      stop(paste("Gene", variableY, "is not in the expression data"))
    }
  }
  if(variableTypeX == "phenotype"){
    if(!variableX %in% colnames(phenotype)){
      stop(paste("Phenotype", variableX, "is not in the phenotype data"))
    }
  }
  if(variableTypeY == "phenotype"){
    if(!variableY %in% colnames(phenotype)){
      stop(paste("Phenotype", variableY, "is not in the phenotype data"))
    }
  }
  
  # Extract the variables.
  xData <- NULL
  yData <- NULL
  if(variableTypeX == "gene"){
    xData <- expression[variableX,]
  }else{
    xData <- phenotype[,variableX]
  }
  if(variableTypeY == "gene"){
    yData <- expression[variableY,]
  }else{
    yData <- phenotype[,variableY]
  }
  xData <- unname(unlist(xData))
  yData <- unname(unlist(yData))
  
  # Determine which function to call
  isContinuousX <- FALSE
  isContinuousY <- FALSE
  if(variableTypeX == "phenotype"){
    if(phenotype_dictionary[which(colnames(phenotype) == variableX)] == "continuous"){
      isContinuousX <- TRUE
    }
  }
  if(variableTypeY == "phenotype"){
    if(phenotype_dictionary[which(colnames(phenotype) == variableY)] == "continuous"){
      isContinuousY <- TRUE
    }
  }

  if(variableTypeX == "gene" || isContinuousX){
    
    # Plot scatter or box.
    if(variableTypeY == "gene" || isContinuousY){
      plotScatter(x = xData, y = yData, xName = variableX, yName = variableY)
    }else{
      plotBoxplot(x = xData, y = yData, xName = variableX, yName = variableY)
    }
  }else{
    
    # Plot box or heatmap
    if(variableTypeY == "gene" || isContinuousY){
      plotBoxplot(x = xData, y = yData, xName = variableX, yName = variableY)
    }else{
      plotHeat(x = xData, y = yData, xName = variableX, yName = variableY)
    }
  }
}

#' A helper function for plotAssociation(). Called when both variables are dichotomous phenotypes.
#' @param x The data for the x axis
#' @param y The data for the y axis
#' @param xName The variable name to use for the x axis
#' @param yName The variable name to use for the y axis
#' @return NULL
plotHeat <- function(x, y, xName, yName){
  
  # Make contingency table.
  tbl <- table(x, y)

  # Plot.
  heatmap(t(tbl), Rowv = NA, Colv = NA, scale = "none", xlab = xName, ylab = yName)
}

#' A helper function for plotAssociation(). Called when one is a dichotomous phenotype and
#' the other is a continuous variable.
#' @param x The data for the x axis
#' @param y The data for the y axis
#' @param xName The variable name to use for the x axis
#' @param yName The variable name to use for the y axis 
#' @return NULL
plotBoxplot <- function(x, y, xName, yName){
  if(is.numeric(x)){
    graphics::boxplot(x ~ y, xlab = xName, ylab = yName, horizontal = TRUE)
  }else{
    graphics::boxplot(y ~ x, xlab = xName, ylab = yName)
  }
}

#' A helper function for plotAssociation(). Called when both variables are continuous.
#' @param x The data for the x axis
#' @param y The data for the y axis
#' @param xName The variable name to use for the x axis
#' @param yName The variable name to use for the y axis
#' @return NULL
plotScatter <- function(x, y, xName, yName){
  plot(x = x, y = y, xlab = xName, ylab = yName)
}

#' Write the table containing the top results across tissues.
#' @param seahorseResultDir The directory containing the SEAHORSE results across tissues
#' @param seahorseInputDir The directory containing input files for SEAHORSE, with the file
#' names corresponding to the file names in seahorseResultDir and each file formatted as
#' a named list of the following: expression, phenotype, dict.
#' @param pathways The pathways in the format used by SEAHORSE
#' @param resultType The name of the result component ("coexpression", "phenotype_association",
#' "phenocor", or "GSEA")
#' @param consolidatedResultFile The path where you wish to save your consolidated file.
#' @param corCutoffSignificance The correlation cutoff significance value. Default is 0.9.
#' @param padjCutoffSignificance The significance cutoff for adjusted p-value.
#' Default is 0.05.
#' @return NULL
writeSigTable <- function(seahorseResultDir, seahorseInputDir, resultType, 
                          consolidatedResultFile, pathways = NULL,
                          corCutoffSignificance = 0.9, padjCutoffSignificance = 0.05){
  
  # Check significance cutoffs.
  if(!is.numeric(corCutoffSignificance)){
    stop("Correlation cutoff must be a numeric value.")
  }
  if(!is.numeric(padjCutoffSignificance)){
    stop("Adjusted p-value significance cutoff must be a numeric value.")
  }
  
  # List all files.
  allFiles <- list.files(seahorseResultDir)
  allFilesInput <- list.files(seahorseInputDir)

  # Save the data.
  if(resultType == "coexpression"){
    
    # Write coexpression data.
    writeGeneGene(inFiles = paste(seahorseResultDir, allFiles, sep = "/"),
                  file = consolidatedResultFile, cutoff = corCutoffSignificance)
    
  }else if(resultType == "phenotype_association"){
    
    # Write the phenotype-gene associations.
    writePhenotypeGeneAssociations(resultFiles = paste(seahorseResultDir, allFiles, sep = "/"), 
                                   inFiles = paste(seahorseInputDir, allFilesInput, sep = "/"), 
                                   corCutoff = corCutoffSignificance,
                                   padjCutoff = padjCutoffSignificance,
                                   file = consolidatedResultFile)
    
  }else if(resultType == "GSEA"){
    
    # Write GSEA associations.
    writePathwayEnrichment(inFiles = paste(seahorseResultDir, allFiles, sep = "/"), 
                           file = consolidatedResultFile, 
                           pathways = pathways, cutoff = padjCutoffSignificance)
    
  }else if(resultType == "phenocor"){
    
    # Write phenotype associations.
    writePhenotypePhenotypeAssociations(inFiles = paste(seahorseInputDir, allFilesInput, sep = "/"),
                                        resultFiles = paste(seahorseResultDir, allFiles, sep = "/"), 
                                        file = consolidatedResultFile,
                                        padjCutoff = padjCutoffSignificance,
                                        corCutoff = corCutoffSignificance)
  }
}

#' Write table of SEAHORSE results filtered by a variable 
#' @param result The full SEAHORSE result
#' @param phenotype The phenotype data in the format used by SEAHORSE
#' @param dictionary The phenotype dictionary in the format used by SEAHORSE
#' @param pathways The pathways in the format used by SEAHORSE
#' @param variable The variable of interest by which to filter (a phenotype or gene)
#' @param variableType Either "gene" or "phenotype"
#' @param resultType The name of the result component ("coexpression", "phenotype_association",
#' "phenocor", or "GSEA")
#' @param tmpFile The temporary file you wish to create for constructing the table.
#' @param resultFile The location where you wish to save your result file.
#' @return NULL
writeTable <- function(result, phenotype = NULL, variable, variableType, resultType, tmpFile,
                       resultFile, dictionary = NULL, pathways = NULL){
  
  # Check formatting.
  if(variableType == "gene" && !(variable %in% colnames(result$coexpression))){
    stop(paste(variable, "not in expression data"))
  }else if(variableType == "phenotype" && !variable %in% names(result$phenotype_association)){
    stop(paste(variable, "not in phenotype data"))
  }
  
  # Obtain the filtered data.
  filteredDat <- result
  if(variableType == "gene" && resultType == "coexpression"){
    
    # Write coexpression data filtered by gene.
    corMat <- result$coexpression[variable,]
    filteredDat$coexpression <- as.matrix(as.data.frame(x = corMat))
    colnames(filteredDat$coexpression) <- variable
    saveRDS(filteredDat, tmpFile)
    writeGeneGene(tmpFile, resultFile, removeLowerTri = FALSE)
    
  }else if(variableType == "phenotype" && resultType == "coexpression"){
    
    # Error if user tries to filter coexpression data by phenotype.
    stop("Cannot filter coexpression data by phenotype")
    
  }else if(variableType == "gene" && resultType == "phenotype_association"){
    
    # Write gene-phenotype associations filtered by gene.
    assoc <- lapply(result$phenotype_association, function(pheno){
      return(pheno[variable,])
    })
    names(assoc) <- names(result$phenotype_association)
    filteredDat$phenotype_association <- assoc
    saveRDS(filteredDat, tmpFile)
    
    # Write input data.
    inputDat <- list(phenotype = phenotype, dict = dictionary)
    tmpFileInput <- paste0(tmpFile, "_phenotype.RDS")
    saveRDS(inputDat, tmpFileInput)
    
    # Write the phenotype-gene associations.
    writePhenotypeGeneAssociations(resultFiles = tmpFile, inFiles = tmpFileInput, resultFile)
    
  }else if(variableType == "phenotype" && resultType == "phenotype_association"){
    
    # Write gene-phenotype associations filtered by phenotype.
    filteredDat$phenotype_association <- list(x = result$phenotype_association[[variable]])
    names(filteredDat$phenotype_association) <- variable
    saveRDS(filteredDat, tmpFile)
    
    # Write input data.
    inputDat <- list(phenotype = phenotype, dict = dictionary)
    tmpFileInput <- paste0(tmpFile, "_phenotype.RDS")
    saveRDS(inputDat, tmpFileInput)
    
    # Write the phenotype-gene associations.
    writePhenotypeGeneAssociations(resultFiles = tmpFile, inFiles = tmpFileInput, resultFile)
    
  }else if(variableType == "gene" && resultType == "GSEA"){
    
    # Error if user tries to filter GSEA results by gene
    stop("Cannot filter GSEA results by gene")
    
  }else if(variableType == "phenotype" && resultType == "GSEA"){
    
    # Write GSEA associations filtered by phenotype.
    filteredDat$GSEA <- list(x = result$GSEA[[variable]])
    names(filteredDat$GSEA) <- variable
    saveRDS(filteredDat, tmpFile)
    writePathwayEnrichment(tmpFile, resultFile, pathways)
    
  }else if(variableType == "gene" && resultType == "phenocor"){
    
    # Error if user tries to filter phenotype correlation results by gene
    stop("Cannot filter phenocor results by gene")
    
  }else if(variableType == "phenotype" && resultType == "phenocor"){
    
    # Subset filtered data.
    otherVars <- setdiff(rownames(result$phenocor$stat), variable)
    filteredDat$phenocor$stat <- rbind(as.matrix(result$phenocor$stat[variable,otherVars]),
                                       as.matrix(result$phenocor$stat[otherVars, variable]))
    filteredDat$phenocor$cramerV <- rbind(as.matrix(result$phenocor$cramerV[variable,otherVars]),
                                          as.matrix(result$phenocor$cramerV[otherVars, variable]))
    filteredDat$phenocor$pval <- rbind(as.matrix(result$phenocor$pval[variable,otherVars]),
                                       as.matrix(result$phenocor$pval[otherVars, variable]))
    filteredDat$phenocor$padj <- rbind(as.matrix(result$phenocor$padj[variable,otherVars]),
                                       as.matrix(result$phenocor$padj[otherVars, variable]))
    filteredDat$phenocor$testType <-rbind(as.matrix(result$phenocor$testType[variable,otherVars]),
                                          as.matrix(result$phenocor$testType[otherVars, variable]))
    row <- rownames(filteredDat$phenocor$stat)
    whichNotNA <- which(!is.na(filteredDat$phenocor$stat))
    filteredDat$phenocor$stat <- as.matrix(filteredDat$phenocor$stat[whichNotNA])
    filteredDat$phenocor$cramerV <- as.matrix(filteredDat$phenocor$cramerV[whichNotNA])
    filteredDat$phenocor$pval <- as.matrix(filteredDat$phenocor$pval[whichNotNA])
    filteredDat$phenocor$padj <- as.matrix(filteredDat$phenocor$padj[whichNotNA])
    filteredDat$phenocor$testType <- as.matrix(filteredDat$phenocor$testType[whichNotNA])
    rowFinal <- row[whichNotNA]
    rownames(filteredDat$phenocor$stat) <-
      rownames(filteredDat$phenocor$cramerV) <-
      rownames(filteredDat$phenocor$pval) <-
      rownames(filteredDat$phenocor$padj) <-
      rownames(filteredDat$phenocor$testType) <- rowFinal
    colnames(filteredDat$phenocor$stat) <-
      colnames(filteredDat$phenocor$cramerV) <-
      colnames(filteredDat$phenocor$pval) <-
      colnames(filteredDat$phenocor$padj) <-
      colnames(filteredDat$phenocor$testType) <- variable
    
    # Build data frame.
    VARNAME2 <- rownames(filteredDat$phenocor$stat)
    TEST <- unname(filteredDat$phenocor$testType)
    TESTSTAT <- unname(filteredDat$phenocor$stat)
    TESTPVALUE <- unname(filteredDat$phenocor$padj)
    VARNAME1 <- variable
    tissueSplit <- strsplit(tmpFile, "/")[[1]]
    tissueDot <- strsplit(tissueSplit[length(tissueSplit)], ".", fixed = TRUE)[[1]]
    tissue <- tissueDot[1]
    dfReturn <- data.frame(VARNAME1 = VARNAME1, VARNAME2 = VARNAME2, tissue = tissue,
                           TEST = TEST, TESTSTAT = TESTSTAT, TESTPVALUE = TESTPVALUE)
    
    # Write File
    con <- gzfile(paste0(resultFile, ".tsv.gz"), "wt")
    on.exit(close(con))
    write.table(dfReturn, file = con, sep = "\t", row.names = FALSE, quote = FALSE)
  }
  
  # Delete temporary file.
  if(file.exists(paste0(tmpFile, ".RDS"))){
    unlink(paste0(tmpFile, ".RDS"))
  }
  if(file.exists(paste0(tmpFile, "_phenotype.RDS"))){
    unlink(paste0(tmpFile, "_phenotype.RDS"))
  }
}
