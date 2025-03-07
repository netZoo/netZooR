#' Plot PANDA network in Cytoscape
#'
#' This function is able to modify PANDA network and plot in Cytoscape. Please make sure that Cytoscape 
#' is installed and open it before calling this function.
#'
#' @param panda.net Character string indicating the input PANDA network in data frame structure type.
#' @param network_name Character string indicating the name of Cytoscape network. 
#' @return PANDA network in Cytoscape
#' @import RCy3
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

