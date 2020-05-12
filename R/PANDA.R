#' Run PANDA in R
#' 
#' \strong{PANDA}(Passing Attributes between Networks for Data Assimilation) is a message-passing model 
#' to reconstruct gene regulatory network, which integrates multiple sources of biological data-including protein-protein interaction data,
#' gene expression data, and transcription factor binding motifs data to reconstruct genome-wide, condition-specific regulatory networks.
#' \href{http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0064832}{[(Glass et al. 2013)])}
#' This function is designed to run the a derived PANDA implement in Python "netZooPy" \href{https://github.com/netZoo/netZooPy}{netZooPy},
#'  which is also a Python library of netZooR.
#'
#' @param expr Character string indicating the file path of expression values file, with each gene(in rows) across samples(in columns).
#' @param motif An optional character string indicating the file path of a prior transcription factor binding motifs dataset.
#'          When this argument is not provided, analysis will continue with Pearson correlation matrix.
#' @param ppi An optional character string indicating the file path of protein-protein interaction edge dataset.
#'          Also, this can be generated with a list of proteins of interest by \code{\link{source.PPI}}.
#' @param mode_process An optional character to define the pre-processing of input dataset, options include "legacy" to the original behavior of PANDA in Python implement,
#'          "interaction" to constrain the TF and gene as interaction across all input datasets, and "union" to take TF union and gene union across all input datasets.
#'          The default value is "union".
#' @param rm_missing When mode_process = "legacy", rm_missing is an optional boolean indicating whether to remove genes and TFs not present in all input files. If TRUE, remove all unmatched TF and genes.
#'         if FALSE, keep all tf and genes. The default value is FALSE.
#'
#' @return A list of three items：
#'          Use \code{$panda} to access the standard output of PANDA as data frame, which consists of four columns: 
#'          "TF", "Gene", "Motif" using 0 or 1 to indicate if this edge belongs to prior motif dataset, and "Score".
#' 
#'          Use \code{$indegree} to access the indegree of PANDA network as data frame, which consists of two columns: "Gene", "Score".
#' 
#'          Use \code{$outdegree} to access the outdegree of PANDA network as data frame, which consists of two columns: "TF", "Score".
#' 
#' @examples 

#' # take the treated TB dataset as example here.
#' # refer to the datasets files path in inst/extdat
#' 
#' treated_expression_file_path <- system.file("extdata", "expr4_matched.txt", package = "netZooR", mustWork = TRUE)
#' motif_file_path <- system.file("extdata", "chip_matched.txt", package = "netZooR", mustWork = TRUE)
#' ppi_file_path <- system.file("extdata", "ppi_matched.txt", package = "netZooR", mustWork = TRUE)
#' 
#' 
#' # Run PANDA for treated and control network
#' treated_all_panda_result <- panda.py(expr = treated_expression_file_path, motif = motif_file_path, ppi = ppi_file_path, mode_process="legacy", rm_missing = TRUE )
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


panda.py <- function( expr, motif=NULL, ppi=NULL, mode_process="union", rm_missing = FALSE){
  
  if(missing(expr)){
    stop("Please provide the gene expression value with option e, e.g. e=\"expression.txt\"") }
  else{ expr.str <- paste("\'", expr, "\'", sep = '') }
  
  if(is.null(motif)){
    motif.str <-  paste('None')
    message("prior motif network is not provided, analysis continues with Pearson correlation matrix.") }
  else{ motif.str <- paste("\'", motif,"\'", sep = '') }
  
  if(is.null(ppi)){
    ppi.str <- paste('None')
    message("No PPI provided.") }
  else{ ppi.str <- paste("\'", ppi, "\'", sep = '') }
  
  # when pre-processing mode is legacy
  if(mode_process == "legacy"){
    
    if(rm_missing == FALSE | missing(rm_missing)){
      message("Use the legacy mode to pre-processing dataset and keep the unmatched tfs or genes")
      arg.str <- "modeProcess ='legacy',remove_missing = False"}
    else { 
      message("Use the legacy mode to pre-processing dataset and keep the matched tfs or genes")
      arg.str <- "modeProcess ='legacy',remove_missing = True"  }
  }
  if(mode_process == "union"){
    arg.str <- "modeProcess ='union'"
  }
  if(mode_process == "intersection"){
    arg.str <- "modeProcess ='intersection'"
  }
  
  # source the pypanda from github raw website.
  reticulate::source_python("https://raw.githubusercontent.com/netZoo/netZooPy/netZoo/panda.py",convert = TRUE)
  
  # invoke py code to create a Panda object
  str <-  paste("panda_obj=Panda(", expr.str, ",", motif.str,",", ppi.str, ",", arg.str, ")", sep ='')
  
  # call py
  py_run_string(str)
  py_run_string("panda_network=pd.DataFrame(panda_obj.export_panda_results,columns=['tf','gene','motif','force'])",local = FALSE, convert = TRUE)
  
  # in-degree of panda network
  py_run_string(paste("indegree=panda_obj.return_panda_indegree()"))
  
  # out-degree of panda netwook
  py_run_string(paste("outdegree=panda_obj.return_panda_outdegree()"))
  
  # return a list with three items-- panda all output data frame, indegree (gene nodes) data frame,
  # and outdegree (tf nodes) data frame.
  # use $panda, $indegree and $outdegree to access each item.
  # 
  # assign the output into three data frames
  panda_net <- py$panda_network
  
  
  # convert the factor to character & character to numeric
  panda_net$tf <- as.character(panda_net$tf)
  panda_net$gene <- as.character(panda_net$gene)
  panda_net$motif <- as.numeric(panda_net$motif)
  panda_net$force <- as.numeric(panda_net$force)
  
  # indegree network
  indegree_net <- py$indegree
  indegree_net <- as.data.frame(cbind(Target = rownames(indegree_net), Target_Score = indegree_net$force), stringsAsFactors =F)
  indegree_net$`Target_Score` <- as.numeric(indegree_net$`Target_Score`)
  
  # outdegree network
  outdegree_net <- py$outdegree
  outdegree_net <- as.data.frame(cbind(Regulator = rownames(outdegree_net), Regulator_Score = outdegree_net$force), stringsAsFactors =F)
  outdegree_net$`Regulator_Score` <- as.numeric(outdegree_net$`Regulator_Score`)
  # rename the PANDA output colnames
  colnames(panda_net) <- c("TF","Gene","Motif","Score")
  
  # check if there is duplicate name of nodes in first two columns
  # if true, prefix the content in regulator column with "reg_" and content in target column with"tar_"
  
  if( length(intersect(panda_net$Gene, panda_net$TF))>0){
    panda_net$TF <- paste('reg_', panda_net$TF, sep='')
    panda_net$Gene <- paste('tar_', panda_net$Gene, sep='')
    message("Rename the content of first two columns with prefix 'reg_' and 'tar_' as there are some duplicate node names between first two columns" )
  }
  output <- list("panda" = panda_net)
  # assign all three network into a list.
  output <- list("panda" = panda_net, "indegree" = indegree_net, "outdegree" = outdegree_net)
  message ("...Finish PANDA...")
  return(output)
}
