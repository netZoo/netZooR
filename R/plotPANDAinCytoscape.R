#' Plot PANDA network in Cytoscape
#'
#'This function is able to modify PANDA network and plot in Cytoscape.
#'
#' @param panda.net Character string indicating the input PANDA network in data frame structure type.
#' @param network.name Character string indicating the name of Cytoscape network. 
#' @examples
#' # refer to the input datasets files of control TB dataset in inst/extdat as example
#' control_expression_file_path <- system.file("extdata", "expr10.txt", package = "netZoo", mustWork = TRUE)
#' motif_file_path <- system.file("extdata", "chip.txt", package = "netZoo", mustWork = TRUE)
#' ppi_file_path <- system.file("extdata", "ppi.txt", package = "netZoo", mustWork = TRUE)
#' 
#' # Run PANDA algorithm
#' control_all_panda_result <- runPanda(e = control_expression_file_path, m = motif_file_path, ppi = ppi_file_path, rm_missing = TRUE )
#' 
#' # access PANDA regulatory network
#' control_net <- control_all_panda_result$panda
#' 
#' # select top 1000 edges in PANDA network by edge weight.
#' panda.net <- head(control_net[order(control_net$force,decreasing = T),], 1000)
#' 
#' # run this function to create a network in Cytoscape.
#' plotPANDAinCytoscape(panda.net, network.name="PANDA")
#' 
#' @return NULL
#' @import RCy3
#' @export 
plotPANDAinCytoscape <- function(panda.net, network.name="PANDA"){
  # launch Cytoscape 3.6.1 or greater
  cytoscapePing ()
  cytoscapeVersionInfo ()
  # change colnames of input PANDA network
  colnames(panda.net) <- c("source","target","interaction","weight")
  # convert the weight to the range between 0 to 10000.
  panda.net$weight_transformed <- rank(panda.net$weight)/length(panda.net$weight)*10000
  if(nrow(panda.net)>6666){
    message("The maximum network objects (nodes+edges) that Cytoscape could support varies in different memories.
            0-20000 objects are suggested by 512M memory size for Cytoscape to view. 
            For better view function in Cytoscape please reduce the size of PANDA network imported")
  }
  
  #create nodes arg for creating a cytoscape plot
  panda_nodes <- data.frame(id=c(panda.net$source,panda.net$target), group=rep(c("TF","Gene"), each=nrow(panda.net)),stringsAsFactors=FALSE)
  # creat cytoscape from DataFrames
  createNetworkFromDataFrames(nodes=panda_nodes,edges=panda.net, title=network.name, collection="DataFrame Example")
  
}

