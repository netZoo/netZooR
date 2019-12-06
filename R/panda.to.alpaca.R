#' panda.to.alpaca
#'
#' \strong{ALPACA}(ALtered Partitions Across Community Architectures) is a method for comparing two genome-scale networks derived from different phenotypic states to identify condition-specific modules.
#' \href{https://www.nature.com/articles/s41540-018-0052-5}{[(Padi and Quackenbush 2018)])}
#' This function compares two networks generate by \code{\link{panda.fast}} in this package and finds the sets of nodes that best characterize the change in modular structure.
#'
#' @param panda_net1 Data Frame indicating one entire network generate by \code{\link{panda.fast}}
#' @param panda_net2 Data Frame indicating another entire network generate by \code{\link{panda.fast}}
#' @param file.stem The folder location and title under which all results will be stored.
#' @param verbose Indicates whether the full differential modularity matrix should also be written to a file. Defaults to FALSE.
#'
#' @return A string message showing the location of output file if file.stem is given, 
#'        or a List where first element is the membership vector and second element is the contribution score of each node to its module's total differential modularity
#' @import Matrix
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
#' treated_panda_net <- panda.fast(e = treated_expression_file_path, m = motif_file_path, ppi = ppi_file_path, rm_missing = TRUE )$panda
#' control_panda_net <- panda.fast(e = control_expression_file_path, m = motif_file_path, ppi = ppi_file_path, rm_missing = TRUE )$panda
#'  
#' # Run ALPACA
#' alpaca<- panda.to.alpaca(treated_panda_net, control_panda_net, "./TB", verbose=TRUE)
#' 
#' @export
panda.to.alpaca <- function(panda_net1, panda_net2, file.stem = "./alpaca", verbose = F){
  # remove "motif" column
  panda_net1 <- panda_net1[, -3]
  colnames(panda_net1) <- c("tf","gene","force")
  panda_net2 <- panda_net2[, -3]
  colnames(panda_net2) <- c("tf","gene","force")
  # ****merge two PANDA network by "tf","motif" column to generate a four-columns data.frame.
  # ****rows in panda_net1 that has no matching row in panda_net2 will have NAs in those columns that are usually filled with values from panda_net2.  
  # net <- merge(panda_net1, panda_net2, by=c("tf","gene","motif"), all.x = T, all.y=T)
  
  # merge two PANDA network by "tf","motif" column to generate a four-columns data.frame.
  net <- merge(panda_net1, panda_net2, by=c("tf","gene"))
  # run ALPACA.
  alp <- alpaca(net, file.stem, verbose)
  # full differential modularity matrix has been written to a file, print out the location.
  message(paste("the ALPACA output located in", file.stem, sep=""))
  return(alp)
  }
  
