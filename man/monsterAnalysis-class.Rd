\name{monsterAnalysis}
\docType{class}
\alias{monsterAnalysis-class}

\title{Class "monsterAnalysis" -- a list containing the transition matrix results from a monster analysis}
\description{This is the main class for MONSTER.}
\section{Extends}{
    Class \code{\link{list}}
}
\note{ 
    Note: This is a summary for reference. For an explanation of the actual usage, see the vignette.
    
    'monsterAnalysis' obj are the main class in use by the MONSTER package. The class extends \code{\link{list}} and provides the key results from the results of the MONSTER analysis.
    
    The transition matrix from the MONSTER analysis is contained in slot 'tm'.  This slot contains an m by m square matrix representing the main observed result from the \code{\link{monster}} method.
    
    For permutation testing, the slot 'nullTM' contains a list of matrices.  The list has length equal to the given number of permutations run with \code{\link{monster}}
    
    Data parameters, numGenes and numSamples indicate the number of genes and number of samples on which \code{\link{monster}} was originally run.
    

}

\section{Methods}{
    Class-specific methods.
    
    \code{\link{transitionNetworkPlot}} plots the transition matrix 
    using the igraph library
    
    \code{\link{transitionPCAPlot}} plots the transcription factors in the 
    transition matrix along the top principal components
    
    \code{\link{dTFIPlot}} plots the differential transcription factor 
    involvement for each transcription factor against the distribution 
    under the null 
    
    \code{\link{calculate.tm.p.values}} computes the significance of 
    each transcription factor involvement.
    
    \code{\link{hcl.heatmap.plot}} plots the transition matrix as a heatmap
    
}

\examples{
    # See vignette
}