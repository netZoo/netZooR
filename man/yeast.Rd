\name{yeast}
\docType{data}
\alias{yeast}
\title{Toy data derived from three gene expression datasets and a mapping from transcription factors to genes.}
\description{
  This data is a list containing gene expression data from three separate yeast studies along with data mapping yeast transcription factors with genes based on the presence of a sequence binding motif for each transcription factor in the vicinity of each gene.  The motif data.frame, yeast$motif, describes a set of pairwise connections where a specific known sequence motif of a transcription factor was found upstream of the corresponding gene.   The expression data, yeast$exp.ko, yeast$exp.cc, and yeast$exp.sr, are three gene expression datasets measured in conditions of gene knockout, cell cycle, and stress response, respectively. 
}
\usage{
data(yeast)
}
\format{A list containing 4 data.frames}
\value{A list of length 4}
\references{
Glass K, Huttenhower C, Quackenbush J, Yuan GC. Passing Messages Between Biological Networks to Refine Predicted Interactions. PLoS One. 2013 May 31;8(5):e64832.
}
\keyword{datasets}
