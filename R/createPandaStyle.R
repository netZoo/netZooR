#' Create a Cytoscape visual style for PANDA network
#'
#'This function is able to create a Cytoscape visual style for any PANDA network output.
#' @param style_name Character string indicating the style name. Defaults to "PandaStyle"
#'
#' @return A visual style in Cytoscape Control Panel under "Style" button.
#' @import RCy3
#' @examples
#' # Here we will load a customized visual style for our network, in which TF 
#' # nodes are orange circles, target gene nodes are blue squares, and edges 
#' # shade and width are the edge weight (likelyhood of regulatory interaction 
#' # between the TF and gene). You can further customize the network style 
#' # directly from Cytoscape.
#' \donttest{
#' sampleNet <- data.frame("TF"=c("TF1", "TF2", "TF3"),
#'   "Gene"=c("gene1", "gene2", "gene3"),"Motif"=NA,
#'   "Score"=c(1,2,3),stringsAsFactors = FALSE)
#' visPandaInCytoscape(sampleNet, network_name="sample")
#' createPandaStyle(style_name="PandaStyle")
#' }
#' @export

createPandaStyle <- function(style_name="PandaStyle"){
  
  # node properties
  nodeShape <- mapVisualProperty('node shape','group','d',c("TF","Gene"),c("ELLIPSE","RECTANGLE"), network = "current")
  nodeColor <- mapVisualProperty('node fill color', 'group','d', c("TF","Gene"), c('#FD7622', '#499df3'))
  nodeBorderColor <- mapVisualProperty("Node Border Paint", 'group', 'd', c("TF","Gene"), c('#FD7622', '#499df3'))
  nodeLabel <- mapVisualProperty("node label","id",'p')

  # edge properties
  edgeLineType <- mapVisualProperty('edge line type','interaction','d',c("1","0"), c("SOLID","LONG_DASH"))
  edgeTargetArrowShape <- mapVisualProperty("edge target arrow shape",'interaction','d',c(0,1),c("DELTA","DELTA"))
  edgeTargetArrowUnSelectedPaint <- mapVisualProperty("Edge Target Arrow Unselected Paint","weight_transformed",'c',c(0,10000),c('#d3d3d3', '#000000'))
  # edge width and edge color shade both represent the edge weight.
  edgeWidth <- mapVisualProperty("edge width", "weight_transformed",'c', c(0,10000), c(0.5,7))
  edgeStrokeUnselectedPaint <- mapVisualProperty('Edge Stroke Unselected Paint',"weight_transformed",'c',c(0,10000), c('#d3d3d3', '#000000'))
  
  # default setting
  defaults <- list(NODE_SIZE=50,
                 EDGE_TRANSPARENCY=255,
                 NODE_LABEL_POSITION="center")
  # create visual style working for PANDA network
  createVisualStyle(style_name,defaults,list(nodeShape, nodeColor, nodeBorderColor,nodeLabel,
                                           edgeLineType, edgeTargetArrowShape, edgeTargetArrowUnSelectedPaint,
                                           edgeWidth, edgeStrokeUnselectedPaint))

}

