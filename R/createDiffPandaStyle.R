#' Create a Cytoscape visual style for two PANDA networks
#'
#'This function is able to create a Cytoscape visual style for two PANDA networks output.
#' @param style_name Character string indicating the style name. Defaults to "Diff.PandaStyle"
#' @param condition_name string vector indicating the same condition name used in \code{\link{panda.diff.edges}} and \code{\link{vis.diff.panda.in.cytoscape}},
#'        and all edges belong to this condition or phenotype will be represented as purple, otherwise, green.
#'
#' @return a visual style in Cytoscape Control Panel under "Style" button.
#' @import RCy3
#' @export

createDiffPandaStyle <- function(style_name="Diff.PandaStyle", condition_name="cond.1"){
  
  # node properties
  nodeShape <- mapVisualProperty('node shape','group','d',c("TF","Gene"),c("ELLIPSE","RECTANGLE"))
  nodeColor <- mapVisualProperty('node fill color', 'group','d', c("TF","Gene"), c('#FD7622', '#499df3'))
  nodeBorderColor <- mapVisualProperty("Node Border Paint", 'group', 'd', c("TF","Gene"), c('#FD7622', '#499df3'))
  nodeLabel <- mapVisualProperty("node label","id",'p')
  
  # edge properties
  edgeLineType <- mapVisualProperty('edge line type','motif','d',c(1,0), c("SOLID","LONG_DASH"))
  edgeTargetArrowShape <- mapVisualProperty("edge target arrow shape",'motif','d',c(1,0),c("DELTA","DELTA"))
  edgeTargetArrowUnSelectedPaint <- mapVisualProperty("Edge Target Arrow Unselected Paint",condition_name,'d',c("T","F"), c('#92278f', '#60bf71'))
  # edge width and edge color shade both represent the edge weight.
  edgeWidth <- mapVisualProperty("edge width", "weight_transformed",'c', c(0,10000), c(0.5,7))
  edgeStrokeUnselectedPaint <- mapVisualProperty('Edge Stroke Unselected Paint',condition_name,'d',c("T","F"), c('#92278f', '#60bf71'))
  
  # default setting
  defaults <- list(NODE_SIZE=50,
                   EDGE_TRANSPARENCY=255,
                   NODE_LABEL_POSITION="center")
  # create visual style working for PANDA network
  createVisualStyle(style_name,defaults,list(nodeShape, nodeColor, nodeBorderColor,nodeLabel,
                                             edgeLineType, edgeTargetArrowShape, edgeTargetArrowUnSelectedPaint,
                                             edgeWidth, edgeStrokeUnselectedPaint))
  
}


