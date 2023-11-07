#' Plot PANDA network in Cytoscape
#'
#' This function is able to modify PANDA network and plot in Cytoscape. Please make sure that Cytoscape 
#' is installed and open it before calling this function.
#'
#' @param panda.net Character string indicating the input PANDA network in data frame structure type.
#' @param network_name Character string indicating the name of Cytoscape network. 
#' @examples
#' \donttest{
#' # refer to the input datasets files of control TB dataset in inst/extdat as example
#' control_expression_file_path <- system.file("extdata", "expr10_matched.txt", 
#' package = "netZooR", mustWork = TRUE)
#' motif_file_path <- system.file("extdata", "chip_matched.txt", package = "netZooR", mustWork = TRUE)
#' ppi_file_path <- system.file("extdata", "ppi_matched.txt", package = "netZooR", mustWork = TRUE)
#' 
#' # Run PANDA algorithm
#' control_all_panda_result <- panda.py(expr = control_expression_file_path, motif = motif_file_path, 
#' ppi = ppi_file_path, mode_process="legacy", rm_missing = TRUE )
#' 
#' # access PANDA regulatory network
#' control_net <- control_all_panda_result$panda
#' 
#' # select top 1000 edges in PANDA network by edge weight.
#' panda.net <- head(control_net[order(control_net$force,decreasing = TRUE),], 1000)
#' 
#' # run this function to create a network in Cytoscape.
#' visPandaInCytoscape(panda.net, network.name="PANDA")
#' }
#' @return PANDA network in Cytoscape
#' @export 
visPandaInCytoscape <- function(panda.net, network_name="PANDA"){
  # launch Cytoscape 3.6.1 or greater
  cytoscapePing ()
  cytoscapeVersionInfo ()
  # change colnames of input PANDA network
  names(panda.net)[names(panda.net) == 'TF'] <- "source"
  names(panda.net)[names(panda.net) == 'Gene'] <- "target"
  names(panda.net)[names(panda.net) == 'Motif'] <- "interaction"
  names(panda.net)[names(panda.net) == 'Score'] <- "weight"

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
  createNetworkFromDataFrames(nodes=panda_nodes,edges=panda.net, title=network_name, collection="DataFrame Example")
  
}

