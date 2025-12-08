#' Run SEAHORSE in R
#' Author(s): Enakshi Saha
#' Description:
#'               SEAHORSE computes gene-gene coexpression matrix,
#'               associations between a set of given phenotypes and each gene,
#'               performs gene set enrichment analysis (GSEA) for each phenotype,
#'               using the measures of association.
#'               GSEA is performed through R package "fgsea".
#'               The measures of association 
#'               for a numerical phenotype is Pearson correlation and
#'               for a categorical phenotype is the p-value of an ANOVA test
#'                
#'
#' Inputs:
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
#'
#' Outputs:
#' @return results    : a list containing three objects
#'         results$coexpression: a gene x gene Pearson correlation matrix.
#'         results$phenotype_association : a list containing a vector for each phenotype
#'         results$GSEA: a list containing a matrix of GSEA results for each phenotype
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
#' phenotype_data$height = 65 + sample.int(10, nrow(phenotype_data), replace = T)
#' 
#' phenotype_dictionary = c("categorical", "numeric")
#' 
#' pathways = list()
#' pathways$pathway1 = sample(rownames(expression_data), 5)
#' pathways$pathway2 = sample(rownames(expression_data), 3)
#' pathways$pathway1 = sample(rownames(expression_data), 7)
#'
#' # Run seahorse
#' results <- seahorse(expression_data, phenotype_data, phenotype_dictionary, pathways)
#'  
#' 
#'  


# Function to run GSEA for a numeric phenotype
#' @export
gsea_numeric <- function(expression, pheno, pathways, results){
  output_seahorse = list()
  output_seahorse$cor = list()
  output_seahorse$GSEA = list()
  
  phenotype_vector = as.numeric(pheno)
  cor = unlist(apply(expression, MARGIN=1, function(x){cor(as.numeric(x), phenotype_vector, use="pairwise.complete.obs")}))
  output_seahorse$cor = cor
  
  # Run GSEA
  cor_rank = sort(cor, decreasing = T)
  fgseaRes <- fgsea(pathways, cor_rank, minSize=15, maxSize=500)
  output_seahorse$GSEA = fgseaRes
  
  return(output_seahorse)
}

# Function to run GSEA for a categorical phenotype
#' @export
gsea_categorical <- function(expression, pheno, pathways, results){
  output_seahorse = list()
  output_seahorse$cor = list()
  output_seahorse$GSEA = list()
  
  phenotype_vector = factor(as.character(pheno))
  cor = unlist(apply(expression, MARGIN=1, function(x){anova(lm(as.numeric(x)~phenotype_vector))$`Pr(>F)`[1]}))
  output_seahorse$cor = cor
  
  # Run GSEA
  cor_rank = sort(cor, decreasing = T)
  fgseaRes <- fgsea(pathways, cor_rank, minSize=15, maxSize=500, scoreType = "pos")
  output_seahorse$GSEA = fgseaRes
  
  return(output_seahorse)
}

# Main SEAHORSE function
#' @export
seahorse <- function(expression, phenotype, phenotype_dictionary, pathways){
  set.seed(0)
  library(fgsea)
  
  results = list()
  
  # Compute coexpression of genes
  results$coexpression = cor(t(expression), use="pairwise.complete.obs")
  
  # Compute association of gene expression with phenotypes and run GSEA
  results$phenotype_association = list()
  results$GSEA = list()
  
  for (i in 1:ncol(phenotype)){
    pheno = phenotype[,i]
    pheno_name = colnames(phenotype)[i]
    
    if (phenotype_dictionary[i] == "numeric"){
      output_seahorse = gsea_numeric(expression, pheno, pathways, results)
    }else {output_seahorse = gsea_categorical(expression, pheno, pathways, results)}
    results$phenotype_association[[pheno_name]] = output_seahorse$cor
    results$GSEA[[pheno_name]] = output_seahorse$GSEA
  }
  return(results)
}