% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BLOBFISH.R
\name{FindSignificantEdgesForHop}
\alias{FindSignificantEdgesForHop}
\title{Find the subnetwork of significant edges n / 2 hops away from each gene.}
\usage{
FindSignificantEdgesForHop(
  geneSet,
  combinedNetwork,
  hopConstraint,
  pValues,
  verbose = FALSE,
  topX = NULL
)
}
\arguments{
\item{geneSet}{A character vector of genes comprising the targets of interest.}

\item{combinedNetwork}{A concatenation of n PANDA-like networks with the following format:
tf,gene,score_net1, score_net2, ... , score_netn}

\item{hopConstraint}{The maximum number of hops to be considered for a gene.}

\item{pValues}{The p-values for all edges.}

\item{verbose}{Whether or not to print detailed information about the run.}

\item{topX}{Select the X lowest significant p-values for each gene. NULL by default.}
}
\value{
A bipartite subnetwork in the same format as the original networks.
}
\description{
Find the subnetwork of significant edges n / 2 hops away from each gene.
}
