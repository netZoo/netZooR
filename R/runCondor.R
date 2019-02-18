#' Turn PANDA network into a CONDOR object
#'
#' \strong{CONDOR} (COmplex Network Description Of Regulators) implements methods for clustering biapartite networks
#' and estimatiing the contribution of each node to its community's modularity, 
#' \href{http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005033}{[(Platig et al. 2016)])}
#' This function uses the result of PANDA algorithm as the input dataset to run CONDOR algorithm. More about \href{https://github.com/jplatig/condor}{condor} package and usage.
#'  
#' @param df Data Frame indicating the result of PANDA regulatory network, created by \code{\link{runPanda}}
#' @param threshold Numeric vector of the customered threshold to select edges. Default value is the the midpoint between 
#' the median edge-weight of prior ( 3rd column "Motif" is 1.0) edges 
#' and the median edge-weight of non-prior edges (3rd column "Motif" is 0.0) in PANDA network, see \code{\link{calculateThreshold}}.
#'
#' @return a CONDOR object, see \code{\link[condor]{create.condor.object}}.
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
#' treated_all_panda_result <- runPanda(e = treated_expression_file_path, m = motif_file_path, ppi = ppi_file_path, rm_missing = TRUE )
#' control_all_panda_result <- runPanda(e = control_expression_file_path, m = motif_file_path, ppi = ppi_file_path, rm_missing = TRUE )
#' 
#' # access PANDA regulatory network
#' treated_net <- treated_all_panda_result$panda
#' control_net <- control_all_panda_result$panda
#' 
#' # Run CONDOR
#' treated_condor_object <- runCondor(treated_net, threshold = 0)
#' control_condor_object <- runCondor(control_net, threshold = 0)
#' 
#' # plot communities
#' # package igraph and package viridisLite are already loaded with this package.
#' 
#' treated_color_num <- max(treated_condor_object$red.memb$com)
#' treated_color <- viridis(treated_color_num, alpha = 1, begin = 0, end = 1, direction = 1, option = "D")
#' condor.plot.communities(treated_condor_object, color_list=treated_color, point.size=0.04, xlab="Target", ylab="Regulator")
#' 
#' control_color_num <- max(control_condor_object$red.memb$com)
#' control_color <- viridis(control_color_num, alpha = 1, begin = 0, end = 1, direction = 1, option = "D")
#' condor.plot.communities(control_condor_object, color_list=control_color , point.size=0.04, xlab="Target", ylab="Regulator")
#' 
#' 
#' @export

#'

runCondor <- function(df, threshold = m){
  
  # *** SELECT EDGE ***
  
  # if the threshold (cutoff) of edge-weight is undefined.
  if (missing(threshold)){
    
    threshold <- calculateThreshold(df)
    message("Using the midway of [median weight of non-prior edges] and [median weight of prior edges], 
            all weights mentioned above are transformationed with formula w'=ln(e^w+1) first")
    
    # transform the edge weight with formula w'=ln(e^w+1) to generate a new column of original data frame
    newdf <- cbind(df,log(exp(df[,4])+1))
    
    # rename the data frame and use cutoff to select edge-weights.
    colnames(newdf)[5] <- c("modifiedForce")
    newdf <- newdf[newdf$modifiedForce >= threshold,c(-3,-5)]
    }
  
  # if the threshold (cutoff) of edge-weight is defined. 
  # when the customed threshold is out of range, print out error message.
   if (threshold > max(df[,4]) || threshold < min(df[,4]) ) {
    stop(paste("Please provide the edge-weight threshold between ", min(df[,4])," and ", max(df[,4])))
  }
  else {
    newdf <- df[df[,4] >= threshold,-3]
  }
  
  
  # *** RUN CONDOR ***
  n_reg <- length(unique(newdf[,1]))
  n_tar <- length(unique(newdf[,2]))
  if(n_reg < n_tar) {
    
    condor.object <- create.condor.object(newdf[,c(2,1)])
  } else { condor.object <- create.condor.object(newdf[,c(1,2)])}
  
  condor.object <- condor.cluster(condor.object, project=F)
  colnames(condor.object$edges)[c(1,2)] <- c ("red","blue")

  return(condor.object)
}



