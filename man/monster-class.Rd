\name{monster}
\docType{class}
\alias{monster-class}

\title{Class "monster" -- a list containing the transition matrix results from a monster analysis}
\description{This is the main class for MONSTER.}
\section{Objects from the Class}{
    Objects should be created with calls to \code{\link{newMRexperiment}}.
}
\section{Extends}{
    Class \code{list}
}
\note{ 
    Note: This is a summary for reference. For an explanation of the actual usage, see the vignette.
    
    MRexperiments are the main class in use by metagenomeSeq. The class extends eSet and provides additional slots which are populated during the analysis pipeline.
    
    MRexperiment dataset are created with calls to \code{\link{newMRexperiment}}.
    MRexperiment datasets contain raw count matrices (integers) accessible through \code{\link{MRcounts}}.
    Similarly, normalized count matrices can be accessed (following normalization) through \code{\link{MRcounts}} by calling norm=TRUE.
    Following an analysis, a matrix of posterior probabilities for counts is accessible through \code{\link{posteriorProbs}}.
    
    The normalization factors used in analysis can be recovered by \code{\link{normFactors}}, as can the library sizes of samples (depths of coverage), \code{\link{libSize}}.
    
    Similarly to other RNASeq bioconductor packages available, the rows of the matrix correspond to a feature (be it OTU, species, gene, etc.) and each column an experimental sample. Pertinent clinical information and potential confounding factors are stored in the phenoData slot (accessed via \code{pData}).
    
    To populate the various slots in an MRexperiment several functions are run.
    1) \code{\link{cumNormStat}} calculates the proper percentile to calculate normalization factors. The cumNormStat slot is populated.
2) \code{\link{cumNorm}} calculates the actual normalization factors using p = cumNormStat.

Other functions will place subsequent matrices (normalized counts (\code{\link{cumNormMat}}), posterior probabilities (\code{\link{posteriorProbs}}))

As mentioned above, \code{MRexperiment} is derived from the virtual class,\code{eSet} and thereby has a \code{phenoData} slot which allows for sample annotation. In the phenoData data frame factors are stored. The normalization factors and library size information is stored in a slot called expSummary that is an annotated data frame and is repopulated for subsetted data.

}

\section{Methods}{
    Class-specific methods.
    \describe{
        \item{\code{[}}{Subset operation, taking two arguments and
            indexing the sample and variable. Returns an \code{MRexperiment object}, including relevant
            metadata. Setting \code{drop=TRUE} generates an error. Subsetting the data, the experiment 
            summary slot is repopulated and pData is repopulated after calling factor (removing levels not present).}
    }
}

\examples{
    # See vignette
}