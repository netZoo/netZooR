#' Source the Protein-Protein interaction in STRING database
#' 
#' This function is able to use a list of Transcription Factors(TF) of interest to source the Protein-Protein interactions (PPI)in STRING database 
#'
#' @param TF a data frame with one column indicating the TF of interest
#' @param STRING.version a numeric vector indicating the STRING version. Default valuve is 10
#' @param species.index a numeric vector indicating NCBI taxonomy identifiers 
#' @param ... any dditional arguments passed to
#'
#' @examples
#' # the example motif file
#' motif_file_path <- system.file("extdata", "chip_matched.txt", package = "netZooR", mustWork = TRUE)
#' motif <- read.table(motif_file_path, sep="\t")
#' # create a TF data frame with one column
#' TF <-data.frame(motif[,1])
#' # create PPI data frame by searching in STRING version 10 
#' # and specifying specie to "Mycobacterium tuberculosis H37Rv".
#' # STRING verison 11 is only accessible to R 4.0.
#' if(R.Version()$major=="3"){PPI <- source.PPI(TF, STRING.version="10", 
#' species.index=83332, score_threshold=0)}
#' if(R.Version()$major=="4"){PPI <- source.PPI(TF, STRING.version="11", 
#' species.index=83332, score_threshold=0)}
#' # write out locally then can be used in \code{\link{panda.py}}.
#' 
#' @return A PPI data.frame which contains three columns: "from" and "to" indicating the direction of protein-protein interaction, and "score" indicating the interaction score between two proteins.
#' @import STRINGdb
#' @export

source.PPI <- function(TF, STRING.version="10", species.index, ...){
  # creat a new STRINGdb class.
  string_db=STRINGdb$new(version=STRING.version, species=as.numeric(species.index),...)
  # change the colname to "TF"
  colnames(TF) <- "TF"
  # map the TF to STRINGdb dataset
  TF_mapped <-  string_db$map(TF,"TF",removeUnmappedRows=F)
  # collect the interactions between the TF of interest
  ppi_tmp <- string_db$get_interactions(TF_mapped$STRING_id)[,c(1,2)]
  # store the PPI by using original identifier.
  PPI <- data.frame(from=TF_mapped[match(ppi_tmp$from,TF_mapped$STRING_id),1], to=TF_mapped[match(ppi_tmp$to,TF_mapped$STRING_id),1])
  # assign "score"column  to "1"
  PPI$score <- "1"
  return(PPI)
}

