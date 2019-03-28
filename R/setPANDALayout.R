#' Plot PANDA network in Cytoscape 
#'
#'This function is able to plot specified amount of egdes in weight descreasing order of PANDA network in Cytoscape.
#'
#' @param panda.net Character string indicating the input PANDA network in data frame structure type.
#' @param top Numeric vector indicating the specified amount of edges to plot. Defaults to top 100 edges.
#' @param network.name Character string indicating the name of Cytoscape network. Defaults to "PANDA".
#'
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
#' # run this function to create a network in Cytoscape.
#' runCytoscapePlot(control_net, top = 200, network.name="TB_control")
#' 
#' @return NULL
#' @import RCy3
#' @export 
runCytoscapePlot <- function(panda.net, top = 100, network.name="PANDA",style="PANDAStyle"){
  # launch Cytoscape 3.6.1 or greater
  cytoscapePing ()
  cytoscapeVersionInfo ()
  # change colnames of input PANDA network
  colnames(panda.net) <- c("source","target","interaction","weight")
  # sort PANDA netwoke and select top "top" edges.(default top 100 edges by weight)
  panda_sorted <- head(panda.net[order(panda.net$weight,decreasing = T),], top)
  #create nodes arg for creating a cytoscape plot
  panda_nodes <- data.frame(id=c(panda_sorted$source,panda_sorted$target), group=rep(c("TF","Gene"), each=top),stringsAsFactors=FALSE)
  # creat cytoscape from DataFrames
  createNetworkFromDataFrames(nodes=panda_nodes,edges=panda_sorted, title=network.name, collection="DataFrame Example")
  setVisualStyle(style)
  #deleteVisualStyle("PANDAStyle")
}

runCytoscapePlot(control_net,top = 1000)

deleteVisualStyle("PANDAStyle")




createPANDAStyle <- function(panda.net){
  colnames(panda.net) <- c("source","target","interaction","weight")
  #sort PANDA netwoke and select top "top" edges.(default top 100 edges by weight)
  #panda_sorted <- head(panda.net[order(panda.net$weight,decreasing = T),], top)
  panda_sorted <- panda.net[order(panda.net$weight,decreasing = T),]
  #create nodes arg for creating a cytoscape plot
  #panda_nodes <- data.frame(id=c(panda_sorted$source,panda_sorted$target), group=rep(c("TF","Gene"), each=top),stringsAsFactors=FALSE)
  panda_nodes <- data.frame(id=c(panda_sorted$source,panda_sorted$target), group=rep(c("TF","Gene"), each=nrow(panda.net)),stringsAsFactors=FALSE)
  max <- max(panda.net$weight)
  min <- min(panda.net$weight)
  
  defaults <- list(NODE_SIZE=50,
                 EDGE_TRANSPARENCY=255,
                 NODE_LABEL_POSITION="center")
  createVisualStyle("PANDAStyle", defaults)
  # set node properties
  setNodeShapeMapping('group',c("TF","Gene"),c("ELLIPSE","RECTANGLE"),style.name = "PANDAStyle")
  setNodeBorderColorMapping('group', c("TF","Gene"), mapping.type = "d",c('#FD7622', '#499df3'),style.name = "PANDAStyle")
  setNodeColorMapping('group', c("TF","Gene"), mapping.type = "d",c('#FD7622', '#499df3'),style.name = "PANDAStyle")
  setNodeLabelMapping("id", style.name = "PANDAStyle")
  # set edge properties
  lockNodeDimensions(F)
  setEdgeColorMapping('weight', c(min,max), c('#dcdcdc', '#808080', '#000000'),style.name = "PANDAStyle")
  setEdgeLineStyleMapping('interaction',c("1","0"), c("SOLID","LONG_DASH"), style.name = "PANDAStyle")
  setEdgeLineWidthMapping('weight',c(min,max), c(3,7,11),style.name = "PANDAStyle")
  setEdgeSourceArrowShapeDefault('ARROW',style.name = "PANDAStyle")
  setEdgeSourceArrowColorMapping('weight',c(min,mean,max), c('#dcdcdc', '#808080', '#000000'),style.name = "PANDAStyle")
}


search()
createPANDAStyle(control_net)



