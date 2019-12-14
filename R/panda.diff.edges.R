
#' Identify differential edges in two PANDA networks
#'
#'To determine the probability that an edge is “different” between the networks, 
#'we first subtracted the z-score weight values estimated by PANDA for the two networks and then determined the value of the inverse cumulative distribution for this difference. 
#'The product of these two probabilities represents the probability than an edge is both “supported” and “different.” 
#'We select edges for which this combined probability is greater than a threshold probability (default value is 0.8).
#'
#' @param net1 vector indicating the PANDA networks of one condition or phenotype.
#' @param net2 vector indicating the PANDA networks of another compared condition orphenotype.
#' @param threshold numerical vector indicating a threshold probability to select select edges.
#' @param condition_name string vector indicating the condition name of net1
#'
#' @return a data.frame with five columns: tf, gene, motif, force and defined condition name("T"in this column means egde belongs to first condition or phenotype,
#' "F" means edge belongs to the second condition or phenotype)
#' @export
#'
panda.diff.edges <- function(net1, net2, threshold=0.8, condition_name="cond.1"){
  
  # reshape two PANDA networks
  merged.net <- merge(net1, net2, by=c("tf","gene"))
  # edge-weight differences
  merged.net$'x-y' <- merged.net$force.x-merged.net$force.y
  merged.net$'y-x' <- merged.net$force.y-merged.net$force.x
  
  # CDF
  fnx <- ecdf(merged.net$force.x)
  fny <- ecdf(merged.net$force.y)
  fnxy <- ecdf(merged.net$'x-y')
  fnyx <- ecdf(merged.net$'y-x')
  
  # filter edge weights
  sub.net1 <- merged.net[fnx(merged.net$force.x)*(fnxy(merged.net$'x-y'))>threshold,]
  sub.net1$cond1 <- "T"
  sub.net1 <- sub.net1[,c(1,2,3,4,9)]
  colnames(sub.net1) <- c("tf","gene","motif","force",condition_name)
  
  sub.net2 <- merged.net[fny(merged.net$force.y)*(fnyx(merged.net$'y-x'))>threshold,]
  sub.net2$cond1 <- "F"
  sub.net2 <- sub.net2[,c(1,2,3,6,9)]
  colnames(sub.net2) <- c("tf","gene","motif","force",condition_name)
  merge.sub.net <- rbind(sub.net1,sub.net2)
  return(merge.sub.net)
}







