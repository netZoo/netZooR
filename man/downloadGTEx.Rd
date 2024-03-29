% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/YARN.R
\name{downloadGTEx}
\alias{downloadGTEx}
\title{Download GTEx files and turn them into ExpressionSet object}
\usage{
downloadGTEx(type = "genes", file = NULL, ...)
}
\arguments{
\item{type}{Type of counts to download - default genes.}

\item{file}{File path and name to automatically save the downloaded GTEx expression set. Saves as a RDS file.}

\item{...}{Does nothing currently.}
}
\value{
Organized ExpressionSet set.
}
\description{
Downloads the V6 GTEx release and turns it into an ExpressionSet object.
}
\examples{
# obj <- downloadGTEx(type='genes',file='~/Desktop/gtex.rds')
}
