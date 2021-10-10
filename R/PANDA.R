#' Run Python implementation PANDA in R
#' 
#' \strong{PANDA}(Passing Attributes between Networks for Data Assimilation) is a message-passing model 
#' to reconstruct gene regulatory network, which integrates multiple sources of biological data-including protein-protein interaction data,
#' gene expression data, and transcription factor binding motifs data to reconstruct genome-wide, condition-specific regulatory networks.
#' \href{http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0064832}{[(Glass et al. 2013)])}
#' This function is designed to run the a derived PANDA implementation in Python Library "netZooPy" \href{https://github.com/netZoo/netZooPy}{netZooPy}.
#'
#' @param expr_file Character string indicating the file path of expression values file, with each gene(in rows) across samples(in columns).
#' @param motif_file An optional character string indicating the file path of a prior transcription factor binding motifs dataset.
#'          When this argument is not provided, analysis will continue with Pearson correlation matrix.
#' @param ppi_file An optional character string indicating the file path of protein-protein interaction edge dataset.
#'          Also, this can be generated with a list of proteins of interest by \code{\link{source.PPI}}.
#'          
#' @param computing 'cpu' uses Central Processing Unit (CPU) to run PANDA; 'gpu' use the Graphical Processing Unit (GPU) to run PANDA. The default value is "cpu".
#' 
#' @param precision 'double' computes the regulatory network in double precision (15 decimal digits); 'single' computes the regulatory network in single precision (7 decimal digits) which is fastaer, requires half the memory but less accurate. The default value is 'double'. 
#' @param save_memory 'TRUE' removes temporary results from memory. The result network is weighted adjacency matrix of size (nTFs, nGenes); 'FALSE' keeps the temporary files in memory. The result network has 4 columns in the form gene - TF - weight in motif prior - PANDA edge. PANDA indegree/outdegree of panda network, only if save_memory = FALSE. The default value is 'FALSE'.
#' @param save_tmp 'TRUE' saves middle data like expression matrix and normalized networks; 'FALSE' deletes the middle data. The default value is 'TURE'.
#' @param keep_expression_matrix 'TRUE' keeps the input expression matrix as an attribute in the result Panda object.'FALSE' deletes the expression matrix attribute in the Panda object. The default value is 'FALSE'.
#' @param modeProcess 'legacy' refers to the processing mode in netZooPy<=0.5, 'union': takes the union of all TFs and genes across priors and fills the missing genes in the priors with zeros; 'intersection': intersects the input genes and TFs across priors and removes the missing TFs/genes. Default values is 'union'.
#' @param remove_missing Only when modeProcess='legacy': remove_missing='TRUE' removes all unmatched TF and genes; remove_missing='FALSE' keeps all tf and genes. The default value is 'FALSE'.
#' 
#' @return When save_memory=FALSE(default), this function will return a list of three items: 
#'          Use \code{$panda} to access the standard output of PANDA as data frame, which consists of four columns: 
#'          "TF", "Gene", "Motif" using 0 or 1 to indicate if this edge belongs to prior motif dataset, and "Score".
#' 
#'          Use \code{$indegree} to access the indegree of PANDA network as data frame, which consists of two columns: "Gene", "Score".
#' 
#'          Use \code{$outdegree} to access the outdegree of PANDA network as data frame, which consists of two columns: "TF", "Score".
#'          
#'          When save_memory=TRUE, this function will return a weigheted adjacency matirx of size (nTFs, nGenes), use \code{$WAMpanda} to access.
#' 
#' @examples 

#' # take the treated TB dataset as example here.
#' # refer to the datasets files path in inst/extdat
#' 
#' treated_expression_file_path <- system.file("extdata", "expr4_matched.txt", 
#' package = "netZooR", mustWork = TRUE)
#' treated_expression_file_path <- system.file("extdata", "expr4_matched.txt",
#'  package = "netZooR", mustWork = TRUE)
#' motif_file_path <- system.file("extdata", "chip_matched.txt", package = "netZooR", mustWork = TRUE)
#' ppi_file_path <- system.file("extdata", "ppi_matched.txt", package = "netZooR", mustWork = TRUE)
#' 
#' 
#' # Run PANDA for treated and control network
#' treated_all_panda_result <- panda.py(expr_file = treated_expression_file_path, 
#' motif_file = motif_file_path, ppi_file = ppi_file_path, 
#' modeProcess="legacy", remove_missing = TRUE )
#' 
#' # access PANDA regulatory network
#' treated_net <- treated_all_panda_result$panda
#' 
#' # access PANDA regulatory indegree network.
#' indegree_net <- treated_all_panda_result$indegree
#' 
#' # access PANDA regulatory outdegree networks
#' outdegree_net <- treated_all_panda_result$outdegree
#' 
#' @import reticulate
#' @export
#'


panda.py <- function(expr_file, motif_file=NULL, ppi_file=NULL, computing="cpu", precision="double",save_memory=FALSE, save_tmp=TRUE, keep_expression_matrix=FALSE, modeProcess="union", remove_missing=FALSE){
  
  if(missing(expr_file)){
    stop("Please provide the path of gene expression data file to 'expr_file' variable") }
  else{ expr.str <- paste("\'", expr_file, "\'", sep = '') }
  
  if(is.null(motif_file)){
    motif.str <- 'None'
    message("The prior motif network is not provided, so analysis continues with Pearson correlation matrix and gene coexpression mastrix is returned as a result network") 
  } else{ motif.str <- paste("\'", motif_file,"\'", sep = '') }
  
  if(is.null(ppi_file)){
    ppi.str <- 'None'
    message("No protein-protein interaction network provided.")
  } else{ ppi.str <- paste("\'", ppi_file, "\'", sep = '') }
  
  # computing variable
  if(computing=="gpu"){
    computing.str <- "computing='gpu'"
  } else{ computing.str <- "computing='cpu'"}
  
  # precision variable
  if(precision == "single" ){
    precision.str <- "precision='single'"
  } else{ precision.str <- "precision='double'"}
  
  # save_memory variable
  if(save_memory==TRUE){
    savememory.str <- "save_memory= True"
  } else{ savememory.str <- "save_memory=False" }
  
  # save_tmp variable
  if(save_tmp==FALSE){
    savetmp.str <- "save_tmp= False"
  } else{ savetmp.str <- "save_tmp=True" }
  
  # keep_expression_matrix variable
  if(keep_expression_matrix==TRUE){
    keepexpression.str <- "keep_expression_matrix=True"
  } else{ keepexpression.str <- "keep_expression_matrix=False" }
  
  # when pre-processing mode is legacy
  if(modeProcess == "legacy"){
    
    if(remove_missing == TRUE){
      message("Use the legacy mode to pre-process the input dataset and keep only the matched TFs or Genes")
      mode.str <- "modeProcess ='legacy', remove_missing = True"  
    } else{ message("Use the legacy mode (netZooPy version <= 0.5) to pre-process the input dataset and keep all the unmatched TFs or Genes")
      mode.str <- "modeProcess ='legacy', remove_missing = False"}
  }
  
  # when pre-processing mode is union
  else if(modeProcess == "union"){
    mode.str <- "modeProcess ='union'"
  }
  
  else if(modeProcess == "intersection"){
    mode.str <- "modeProcess ='intersection'"
  }
  
  # source the pypanda from github raw website.
  pandapath <- system.file("extdata", "panda.py", package = "netZooR", mustWork = TRUE)
  reticulate::source_python(pandapath,convert = TRUE)
  
  # invoke Python script to create a Panda object
  obj.str <-  paste("panda_obj=Panda(", expr.str, ",", motif.str,",", ppi.str, ",", computing.str, ",", precision.str, ",", savememory.str, ",", savetmp.str, "," , keepexpression.str, ",",  mode.str, ")", sep ='')
  
  # run Python code
  py_run_string(obj.str)
  # run PAMDA
  if(save_memory == FALSE){
    
    py_run_string("panda_network=panda_obj.export_panda_results",local = FALSE, convert = TRUE)
    # convert python object to R vector
    panda_net <- py$panda_network
    
    # re-assign data type
    panda_net$tf <- as.character(panda_net$tf)
    panda_net$gene <- as.character(panda_net$gene)
    panda_net$motif <- as.numeric(panda_net$motif)
    panda_net$force <- as.numeric(panda_net$force)
    # adjust column order
    panda_net <- panda_net[,c("tf","gene","motif","force")]
    # rename the PANDA output colnames
    colnames(panda_net) <- c("TF","Gene","Motif","Score")
    
    
    # in-degree of panda network
    py_run_string(paste("indegree=panda_obj.return_panda_indegree()"))
    indegree_net <- py$indegree
    indegree_net <- as.data.frame(cbind(Target = rownames(indegree_net), Target_Score = indegree_net$force), stringsAsFactors =FALSE)
    indegree_net$`Target_Score` <- as.numeric(indegree_net$`Target_Score`)
    
    # out-degree of panda netwook
    py_run_string(paste("outdegree=panda_obj.return_panda_outdegree()"))
    outdegree_net <- py$outdegree
    outdegree_net <- as.data.frame(cbind(Regulator = rownames(outdegree_net), Regulator_Score = outdegree_net$force), stringsAsFactors =FALSE)
    outdegree_net$`Regulator_Score` <- as.numeric(outdegree_net$`Regulator_Score`)
    
    if( length(intersect(panda_net$Gene, panda_net$TF))>0){
      panda_net$TF <- paste('reg_', panda_net$TF, sep='')
      panda_net$Gene <- paste('tar_', panda_net$Gene, sep='')
      message("Rename the content of first two columns with prefix 'reg_' and 'tar_' as there are some duplicate node names between the first two columns" )
    }
    
    output <- list("panda" = panda_net, "indegree" = indegree_net, "outdegree" = outdegree_net)
    
    
  } else{ py_run_string("panda_network=panda_obj.panda_network",local = FALSE, convert = TRUE) 
    panda_net <- py$panda_network
    # weighted adjacency matrix of PANDA network 
    output <- list("WAMpanda" = panda_net)
  }
  
  message ("...Finish PANDA...")
  return(output)
}  


