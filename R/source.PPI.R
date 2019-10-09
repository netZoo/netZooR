#' Source the Protein-Protein interaction in STRING database
#' 
#' This function is able to use a list of Transcription Factors(TF) of interest to source the Protein-Protein interactions (PPI)in STRING database 
#'
#' @param TF Data frame with one column indicating the TF of interest
#' @param species.index Numeric vector indicating NCBI taxonomy identifiers 
#' @param ... any dditional arguments passing to
#'
#' @examples
#' # the example motif file
#' motif_file_path <- system.file("extdata", "chip_matched.txt", package = "netZooR", mustWork = TRUE)
#' motif <- read.table(motif_file_path, sep="\t")
#' # create a TF data frame with one column
#' TF <-data.frame(motif[,1])
#' # create PPI data frame
#' PPI <- source.PPI(TF,species.index=83332, score_threshold=0)
#' # write out locally then can be used in \code{\link{runPanda}}.
#' 
#' @return A PPI data frame
#' @import STRINGdb
#' @export

source.PPI <- function(TF, species.index, ...){
  # creat a new STRINGdb class.
  string_db=STRINGdb$new(version="10",species=species.index,...)
  # change the colname to "TF"
  colnames(TF) <- c("TF")
  # map the TF to STRINGdb dataset
  TF_mapped <-  string_db$map(TF,"TF",removeUnmappedRows=F)
  # collect the interactions between the TF of interest
  PPI <- string_db$get_interactions(TF_mapped$STRING_id)[,c(1,2)]
  # remove the species index in the string identifiers.
  PPI$from <- gsub(paste(species.index,".",sep=""), "", PPI$from)
  PPI$to <- gsub(paste(species.index,".",sep=""), "", PPI$to)
  # return a PPI network based on TF list of interest
  return(PPI)
}

