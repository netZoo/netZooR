#' Plot PANDA network
#'
#'This function is able to plot specified amount of egdes in PANDA network, after generate PANDA network \code{\link{runPanda}}.
#'
#' @param top Numeric vector indicating the amount of edges selected to plot by decreasting order of egde weights. Defaults to 100.
#' @param file Character string indicating the name of output .png file. Defaults to 'panda_top_100.png'.
#'
#' @return a message showing the path of output plot file. 
#'
#' @examples 
#' # refer to the input datasets files of control in inst/extdat as example
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
#' # plot the top 100 edges of the PANDA network
#' plotPanda(top =100, file="top100_panda.png")
#' 
#' @import reticulate
#' @export
plotPanda <- function(top = 100, file = 'panda_top_100.png'){
  # source analyze_panda.py from GitHub
  reticulate::source_python("https://raw.githubusercontent.com/twangxxx/pypanda/master/pypanda/analyze_panda.py",convert = TRUE)
  # run py code to create an instance named "plot" of AnalyzePanda class.
  py_run_string(paste("plot = AnalyzePanda(p)"))
  # invoke the method "top_network_plot" to plot the network. 
  py_run_string(paste("plot.top_network_plot(top=", top,",","file='",file,"\')",sep=""))
  # print out a message to indicate the path of output .png file.
  message(paste("The plot located in ", getwd(), "/", file, sep = ""))
}


# plotPanda(top =10, file="test3.png")
