% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LIONESS.R
\name{lioness}
\alias{lioness}
\title{Compute LIONESS (Linear Interpolation to Obtain Network Estimates for Single Samples)}
\usage{
lioness(
  expr,
  motif = NULL,
  ppi = NULL,
  network.inference.method = "panda",
  ncores = 1,
  ...
)
}
\arguments{
\item{expr}{A mandatory expression dataset, as a genes (rows) by samples (columns) data.frame}

\item{motif}{A motif dataset, a data.frame, matrix or exprSet containing 3 columns.
Each row describes an motif associated with a transcription factor (column 1) a
gene (column 2) and a score (column 3) for the motif.}

\item{ppi}{A Protein-Protein interaction dataset, a data.frame containing 3 columns.
Each row describes a protein-protein interaction between transcription factor 1(column 1),
transcription factor 2 (column 2) and a score (column 3) for the interaction.}

\item{network.inference.method}{String specifying choice of network inference method. Default is "panda".
Options include "pearson".}

\item{ncores}{int specifying the number of cores to be used. Default is 1. 
(Note: constructing panda networks can be memory-intensive, and the number of cores should take into consideration available memory.)}

\item{...}{additional arguments for panda analysis}
}
\value{
A list of length N, containing objects of class "panda" 
corresponding to each of the N samples in the expression data set.\cr
"regNet" is the regulatory network\cr
"coregNet" is the coregulatory network\cr
"coopNet" is the cooperative network
}
\description{
Compute LIONESS (Linear Interpolation to Obtain Network Estimates for Single Samples)
}
\examples{
data(pandaToyData)
lionessRes <- lioness(expr = pandaToyData$expression[,1:3], motif = pandaToyData$motif, 
  ppi = pandaToyData$ppi,hamming=1,progress=FALSE)
}
\references{
Kuijjer, M.L., Tung, M., Yuan, G., Quackenbush, J. and Glass, K., 2015. 
Estimating sample-specific regulatory networks. arXiv preprint arXiv:1505.06440.
Kuijjer, M.L., Hsieh, PH., Quackenbush, J. et al. lionessR: single sample network inference in R. BMC Cancer 19, 1003 (2019). https://doi.org/10.1186/s12885-019-6235-7
}
\keyword{keywords}
