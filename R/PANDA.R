#' Run PANDA in R
#' 
#' \strong{PANDA}(Passing Attributes between Networks for Data Assimilation) is a message-passing model 
#' to reconstruct gene regulatory network. It integrates multiple sources of biological data-including protein-protein interaction,
#' gene expression data, and transccription factor binding motifs data-to reconstruct genome-wide, condition-specific regulatory networks.
#' \href{http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0064832}{[(Glass et al. 2013)])}
#' This function is designed to run the a derived PANDA algorithm in Python from \href{https://github.com/netZoo/netZooPy}{netZooPy},
#'  which is also a Python library of netZooR.
#'
#' @param e Character string indicatining the file path of expression values file, with each gene (in row) by samples (in columns)
#' @param m An optional character string indicatining the file path of pair file of a prior transcription factor binding motifs dataset.
#'          When this argument is not provided, analysis will continue with Pearson correlation matrix.
#' @param ppi An optional character string indicatining the file path of protein-protein interaction edge dataset.
#'          Also this can be generated given a list of proteins of interest by \code{\link{source.PPI}}.
#' @param rm_missing an optional boolean indicatining whether to remove genes and tfs not present in all input files. If TRUE, remove all unmatched tf and genes.
#'         if FALSE, keep all tf and genes. The default value is FALSE.
#'
#' @return A list of three items：
#'          Use \code{$panda} to access the standard output of PANDA network in data.frame, which consists of four columns: 
#'          "tf", "gene", "motif" 0 or 1 to indicate if this edge belongs to prior motif dataset, and "force".
#' 
#'          Use \code{$indegree} to access the indegree of PANDA network in data.frame, consisting of two columns: "gene", "force".
#' 
#'          Use \code{$outdegree} to access the outdegree of PANDA network in data.fra,e, consisting of two columns: "tf", "force".
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
#' treated_all_panda_result <- panda(e = treated_expression_file_path, m = motif_file_path, ppi = ppi_file_path, rm_missing = TRUE )
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


panda <- function( e = expression, m = motif, ppi = ppi, rm_missing = FALSE){
  
  if(missing(e)){
    stop("Please provide the gene expression value with option e, e.g. e=\"expression.txt\"") }
  else{ str1 <- paste("\'", e, "\'", sep = '') }
  
  if(missing(m)){
    str2 <-  paste('None')
    message("Pair file of motif edges is not provided, analysis continues with Pearson correlation matrix.") }
  else{ str2 <- paste("\'", m,"\'", sep = '') }
  
  if(missing(ppi)){
    str3 <- paste('None')
    message("No PPI provided.") }
  else{ str3 <- paste("\'", ppi, "\'", sep = '') }
  
  if(rm_missing == FALSE | missing(e)){
    str4 <- paste('False')
    message("Not removing missing values") }
  else { str4 <- paste('True') }
  
  # source the pypanda from github raw website.
  reticulate::source_python("https://raw.githubusercontent.com/netZoo/netZooPy/netZoo/panda.py",convert = TRUE)
  
  # invoke py code to create a Panda object
  str <-  paste("panda_obj=Panda(", str1, ",", str2,",", str3, ",", "remove_missing=", str4, ")", sep ='')
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
  
  # assign the output into three data frames
  panda_net <- py$panda_network
  # convert the character to numeric
  panda_net$motif <- as.numeric(panda_net$motif)
  panda_net$force <- as.numeric(panda_net$force)
  indegree_net <- py$indegree
  outdegree_net <- py$outdegree
  # check if there is duplicate name of nodes in first two columns
  # if true, prefix the content in regulator column with "reg_" and content in target column with"tar_"
  
  if( length(intersect(panda_net[, 1], panda_net[, 2]))>0){
    panda_net[,1] <-paste('reg_', panda_net[,1], sep='')
    panda_net[,2] <-paste('tar_', panda_net[,2], sep='')
    colnames(indegree_net)<- paste("tar_",colnames(indegree_net), sep='')
    colnames(outdegree_net)<- paste("reg_",colnames(outdegree_net), sep='')
    message("Rename the context of first two columns with prefix 'reg_' and 'tar_', as there are some duplicate node name between first two columns" )
  }
  
  # assign all three network into a list.
  output <- list("panda" = panda_net, "indegree" = indegree_net, "outdegree" = outdegree_net)
  message ("...Finish PANDA run...")
  return(output)
}

