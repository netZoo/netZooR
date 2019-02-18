#' Run ALPACA
#'
#' \strong{ALPACA}(ALtered Partitions Across Community Architectures) is a method for comparing two genome-scale networks derived from different phenotypic states to identify condition-specific modules.
#' \href{https://www.nature.com/articles/s41540-018-0052-5}{[(Padi and Quackenbush 2018)])}
#' This function compares two networks generate by \code{\link{runPanda}} in this package and finds the sets of nodes that best characterize the change in modular structure.
#'
#' @param panda_net1 Data Frame indicating one entire network generate by \code{\link{runPanda}}
#' @param panda_net2 Data Frame indicating another entire network generate by \code{\link{runPanda}}
#' @param file.stem The folder location and title under which all results will be stored.
#' @param verbose Indicates whether the full differential modularity matrix should also be written to a file. Defaults to FALSE.
#'
#' @return A string message showing the location of output file if file.stem is given, 
#'        or a List where first element is the membership vector and second element is the contribution score of each node to its module's total differential modularity
#' @import condor
#' @import Matrix
#' @import GOstats
#' @import org.Hs.eg.db
#' @import GO.db
#'
#' @examples
#' # refer to four input datasets files in inst/extdat
#' treated_expression_file_path <- system.file("extdata", "expr4.txt", package = "netZoo", mustWork = TRUE)
#' control_expression_file_path <- system.file("extdata", "expr10.txt", package = "netZoo", mustWork = TRUE)
#' motif_file_path <- system.file("extdata", "chip.txt", package = "netZoo", mustWork = TRUE)
#' ppi_file_path <- system.file("extdata", "ppi.txt", package = "netZoo", mustWork = TRUE)
#' 
#' 
#' # Run PANDA for treated and control network
#' treated_panda_net <- runPanda(e = treated_expression_file_path, m = motif_file_path, ppi = ppi_file_path, rm_missing = TRUE )$panda
#' control_panda_net <- runPanda(e = control_expression_file_path, m = motif_file_path, ppi = ppi_file_path, rm_missing = TRUE )$panda
#'  
#' # Run ALPACA
#' alpaca<- runAlpaca(treated_panda_net, control_panda_net, "~/Desktop/TB", verbose=T)
#' 
#' @export
runAlpaca <- function(panda_net1, panda_net2, file.stem, verbose=F){
  # remove "motif" column
  panda_net1 <- panda_net1[, -3]
  panda_net2 <- panda_net2[, -3]
  # ****merge two PANDA network by "tf","motif" column to generate a four-columns data.frame.
  # ****rows in panda_net1 that has no matching row in panda_net2 will have NAs in those columns that are usually filled with values from panda_net2.  
  # net <- merge(panda_net1, panda_net2, by=c("tf","gene","motif"), all.x = T, all.y=T)
  
  # merge two PANDA network by "tf","motif" column to generate a four-columns data.frame.
  net <- merge(panda_net1, panda_net2, by=c("tf","gene"))
  # run ALPACA.
  alp <- alpaca(net, file.stem, verbose)
  # if full differential modularity matrix has been written to a file, print out the location.
  if(exists(file.stem)) {
    return(message(passte("the ALPACA output located in", file.stem, sep="")))
  }
  
  # if not return the list of ALPACA output.
  else{
    return(alp)
  }
  
}