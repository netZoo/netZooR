#' Identify differential edges in two PANDA networks
#'
#'To determine the probability that an edge is "different" between the networks, 
#'we first subtracted the z-score weight values estimated by PANDA for the two networks and then determined the value of the inverse cumulative distribution for this difference. 
#'The product of these two probabilities represents the probability than an edge is both "supported" and "different."
#'We select edges for which this combined probability is greater than a threshold probability (default value is 0.8).
#'
#' @param panda.net1 vector indicating the PANDA networks of one condition or phenotype.
#' @param panda.net2 vector indicating the PANDA networks of another compared condition orphenotype.
#' @param threshold numerical vector indicating a threshold probability to select select edges.
#' @param condition_name string vector indicating the condition name of net1
#'
#' @return a data.frame with five columns: tf, gene, motif, Score and defined condition name(the row with "T" in this column means this egde belongs to first condition or phenotype,
#' "F" means this edge belongs to the second condition or phenotype)
#' @examples 
#' # refer to four input datasets files in inst/extdat
#' treated_expression_file_path <- system.file("extdata", "expr4_matched.txt", 
#' package = "netZooR", mustWork = TRUE)
#' control_expression_file_path <- system.file("extdata", "expr10_matched.txt", 
#' package = "netZooR", mustWork = TRUE)
#' motif_file_path <- system.file("extdata", "chip_matched.txt", package = "netZooR", mustWork = TRUE)
#' ppi_file_path <- system.file("extdata", "ppi_matched.txt", package = "netZooR", mustWork = TRUE)
#' 
#' 
#' # Run PANDA for treated and control network
#' #treated_all_panda_result <- pandaPy(expr_file = treated_expression_file_path, 
#' #motif_file= motif_file_path, ppi_file = ppi_file_path, modeProcess="legacy", remove_missing = TRUE )
#' #control_all_panda_result <- pandaPy(expr_file = control_expression_file_path, 
#' #motif_file= motif_file_path, ppi_file= ppi_file_path, modeProcess="legacy", remove_missing = TRUE )
#' 
#' # access PANDA regulatory network
#' #treated_net <- treated_all_panda_result$panda
#' #control_net <- control_all_panda_result$panda
#' 
#' #merged.panda <- pandaDiffEdges(treated_net, control_net, condition_name="treated")
#' 
#' @export
#'

pandaDiffEdges <- function(panda.net1, panda.net2, threshold=0.8, condition_name="cond.1"){
  
  # reshape two PANDA networks
  merged.net <- merge(panda.net1, panda.net2, by=c("TF","Gene"))
  # edge-weight differences
  merged.net$'x-y' <- merged.net$Score.x-merged.net$Score.y
  merged.net$'y-x' <- merged.net$Score.y-merged.net$Score.x
  
  # CDF
  fnx <- ecdf(merged.net$'Score.x')
  fny <- ecdf(merged.net$'Score.y')
  fnxy <- ecdf(merged.net$'x-y')
  fnyx <- ecdf(merged.net$'y-x')
  
  # filter edge weights
  sub.net1 <- merged.net[fnx(merged.net$Score.x)*(fnxy(merged.net$'x-y'))>threshold,]
  sub.net1[,c(condition_name)] <- "T"
  sub.net1 <- sub.net1[,c("TF","Gene","Motif.x","Score.x",condition_name)]

  colnames(sub.net1) <- c("TF","Gene","Motif","Score",condition_name)
  
  sub.net2 <- merged.net[fny(merged.net$Score.y)*(fnyx(merged.net$'y-x'))>threshold,]
  sub.net2[,c(condition_name)] <- "F"
  sub.net2 <- sub.net2[,c("TF","Gene","Motif.y","Score.y", condition_name)]
  colnames(sub.net2) <- c("TF","Gene","Motif","Score",condition_name)
  
  # merge two panda networks
  merge.sub.net <- rbind(sub.net1,sub.net2)
  return(merge.sub.net)
}






