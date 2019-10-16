#' Run pypanda in R
#' 
#' \strong{PANDA}(Passing Attributes between Networks for Data Assimilation) is a message-passing model 
#' to gene regulatory network reconstruction. It integrates multiple sources of biological data, 
#' including protein-protein interaction, gene expression, and sequence motif information,
#' in order to reconstruct genome-wide, condition-specific regulatory networks.
#' \href{http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0064832}{[(Glass et al. 2013)])}
#' This function is able to run \href{https://github.com/aless80/pypanda}{pypanda} -- Python implementation of PANDA in R enviroment.
#'
#' @param e Character String indicatining the file path of expression values file, as each gene (row) by samples (columns) \emph{required}
#' @param m Character String indicatining the file path of pair file of motif edges,
#'          when not provided, analysis continues with Pearson correlation matrix. \emph{optional}
#' @param ppi Character String indicatining the pair file path of Protein-Protein interaction dataset. \emph{optional}
#' @param rm_missing Boolean indicatining whether to remove missing values. If TRUE, removes missing values.
#'         if FALSE, keep missing values. THe default value is FALSE. \emph{optional}
#'
#' @return List of three items：
#'  Use \code{$panda} to access the entire data frame of PANDA output network consisting of four columns: 
#' "tf", "gene", "motif", and "force".
#' 
#' Use \code{$indegree} to access the data frame of indegree of PANDA network, consisting of two columns: "gene", "force".
#' 
#' Use \code{$outdegree} to access the data frame of outdegree of PANDA network, consisting of two columns: "tf", "force".
#' 
#' \strong{Note}: If there is any duplicate name of nodes in first two columns,
#' prefix the content in regulator column with "reg_" and content in target column with"tar_"
#' 
#' Please re-name the colnames if needed.
#' 
#' @examples 
#' # refer to four input datasets files in inst/extdat
#' treated_expression_file_path <- system.file("extdata", "expr4_matched.txt", package = "netZooR", mustWork = TRUE)
#' control_expression_file_path <- system.file("extdata", "expr10_matched.txt", package = "netZooR", mustWork = TRUE)
#' motif_file_path <- system.file("extdata", "chip_matched.txt", package = "netZooR", mustWork = TRUE)
#' ppi_file_path <- system.file("extdata", "ppi_matched.txt", package = "netZooR", mustWork = TRUE)
#' 
#' 
#' # Run PANDA for treated and control network
#' treated_all_panda_result <- panda(e = treated_expression_file_path, m = motif_file_path, ppi = ppi_file_path, rm_missing = TRUE )
#' control_all_panda_result <- panda(e = control_expression_file_path, m = motif_file_path, ppi = ppi_file_path, rm_missing = TRUE )
#' 
#' # access PANDA regulatory network
#' treated_net <- treated_all_panda_result$panda
#' control_net <- control_all_panda_result$panda
#' 
#' # access PANDA regulatory indegree network result
#' indegree_net <- treated_all_panda_result$indegree
#' 
#' # access PANDA regulatory outdegree network result
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
    message("Rename the context of first two columns with prefix 'reg_' and 'tar_'" )
  }
  
  # assign all three network into a list.
  output <- list("panda" = panda_net, "indegree" = indegree_net, "outdegree" = outdegree_net)
  message ("...Finish PANDA run...")
  return(output)
}
