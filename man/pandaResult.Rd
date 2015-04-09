\name{pandaResult}
\docType{data}
\alias{pandaResult}
\title{Analysis result from PANDA algorithm on toy data}
\description{
  This data panda object resulting from running the PANDA algorithm on the supplied toy dataset.
  data(pandaToyData)
  pandaResult <- panda(pandaToyData$motif, pandaToyData$expression,pandaToyData$ppi,hamming=.1,progress=TRUE)
}
\usage{pandaResult}
\format{A panda object}
\value{A panda object}
\source{}
\references{
Glass K, Huttenhower C, Quackenbush J, Yuan GC. Passing Messages Between Biological Networks to Refine Predicted Interactions. PLoS One. 2013 May 31;8(5):e64832.
}
\keyword{datasets}
