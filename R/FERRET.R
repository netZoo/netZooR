#'  FERRET_Results class
#'  This object type represents the list of results in a given directory.
#'  
#'  @name FERRET_Results-class
#'  @rdname FERRET_Results-class
#'  @exportClass FERRET_Results
#'  @slot results A named list of inferred networks, where each network is represented as adjacency list.
#'  Names correspond to the input file.
#'  @slot directory The name of the source directory.
#'  @slot interpretationOfNegative Whether to interpret a negative weight as a regulatory
#'  relationship with poor evidence or as an inhibitory relationship.
methods::setClass(
  Class="FERRET_Results",
  representation(results = "list",
                 directory = "character",
                 interpretationOfNegative = "character")
)

#'  FERRET_Comparisons class
#'  This object type represents the list of comparisons to make for the ingroup
#'  and the outgroup.
#'  
#'  @name FERRET_Comparisons-class
#'  @rdname FERRET_Comparisons-class
#'  @exportClass FERRET_Comparisons
#'  @slot ingroup A list of FERRET_Comparison objects
#'  @slot outgroup A list of FERRET_Comparison objects
methods::setClass(
  Class="FERRET_Comparisons",
  representation(ingroup = "list",
                 outgroup = "list")
)

#'  FERRET_Comparison class
#'  This object type represents a single comparison as a source-target pair.
#'  
#'  @name FERRET_Comparison-class
#'  @rdname FERRET_Comparison-class
#'  @exportClass FERRET_Comparison
#'  @slot source The name of a source network as represented in a FERRET_Results object.
#'  @slot target The name of a target network as represented in a FERRET_Results object.
methods::setClass(
  Class="FERRET_Comparison",
  representation(source = "character",
                 target = "character")
)

#'  FERRET_Similarities class
#'  This object type represents the list of similarities across cutoffs for
#'  the ingroup and outgroup.
#'  
#'  @name FERRET_Similarities-class
#'  @rdname FERRET_Similarities-class
#'  @exportClass FERRET_Similarities
#'  @slot ingroup A vector of similarity values across cutoffs for the ingroup.
#'  @slot outgroup A vector of similarity values across cutoffs for the outgroup.
methods::setClass(
  Class="FERRET_Similarities",
  representation(ingroup = "numeric",
                 outgroup = "numeric")
)

#'  FERRET_ROC_AUC
#'  This object type includes ROC and AUC scores.
#'  
#'  @name FERRET_ROC_AUC-class
#'  @rdname FERRET_ROC_AUC-class
#'  @exportClass FERRET_ROC_AUC
#'  @slot auc A scalar AUC score.
#'  @slot roc A named list of results, where the names are the types of metrics
#'  and each item in the list is an object of type FERRET_Similarities.
methods::setClass(
  Class="FERRET_ROC_AUC",
  representation(auc = "numeric",
                 roc = "list")
)


#' Download the data to use in FERRET. This method downloads data from an in-built URL.
#' Users must specify the data set of interest and a path.
#' @param datasetName One of "CPTAC_RCC" or "CPTAC_GBM".
#' @param destinationDir A directory where the data should be saved.
#' @export 
DownloadFERRETData <- function(datasetName, destinationDir){
  # First, make sure input is valid.
  validDatasetUrls <- c(CPTAC_RCC = "https://zenodo.org/api/records/10854694/files-archive", 
                        CPTAC_GBM = "https://zenodo.org/api/records/10850845/files-archive")

  # Ensure that users specify a valid dataset to download.
  # If users did not specify a valid dataset, stop.
  if(datasetName %in% names(validDatasetUrls)){
    
    # If the directory doesn't exist, create it.
    if(!file.exists(destinationDir)){
      tryCatch({
        dir.create(destinationDir)
      }, error = function(cond){
        stop(paste("Could not create directory ", destinationDir))
      })
    }
    
    # Download the data.
    zipLoc <- paste0(destinationDir, "/data.zip")
    unzipLoc <- paste0(destinationDir, "/data/")
    options(timeout = 6000)
    utils::download.file(validDatasetUrls[datasetName], zipLoc)
    if(!file.exists(unzipLoc)){
      dir.create(unzipLoc)
    }
    utils::unzip(zipfile = zipLoc, exdir = unzipLoc)
    
    message(paste("Your data are ready in", destinationDir))
  }else{
    stop(paste("Valid datasets include the following: ", paste(names(validDatasetUrls), collapse = ", ")))
  }
}

#' This function loads all results from a given directory and stores them in
#' a FERRET_Results object. Results are expected to be in an adjacency list
#' format, where the first column is the list transcription factors regulating genes,
#' the second column is the list of regulated genes, and the third column is the
#' score for each regulatory relationship.
#' @param resultDirectory The name of the directory.
#' @param interpretationOfNegative Negatives can refer to either poor associations ("poor")
#' or inhibitory ("inhibitory") regulation, depending on the meaning.
#' @param firstColumnIsRowname Whether to consider the first column as the list of
#' row names.
#' @returns An object of type FERRET_Results.
LoadResults <- function(resultDirectory, interpretationOfNegative = "poor",
                        firstColumnIsRowname = TRUE){
  
  # Check that the interpretation of negative is valid.
  if(interpretationOfNegative != "poor" && interpretationOfNegative != "inhibitory"){
    stop(paste("Valid values for interpretationOfNegative are 'poor' or 'inhibitory.'",
               "You entered:", interpretationOfNegative))
  }
  
  # Check that the directory exists.
  if(!file.exists(resultDirectory)){
    stop(paste("Invalid directory:", resultDirectory))
  }
  
  # Get the file names.
  fileNames <- list.files(resultDirectory)
  
  # Check that the directory is not empty.
  if(length(fileNames) == 0){
    stop(paste("The directory", resultDirectory, "is empty!"))
  }
  
  # If the interpretation of a negative is "poor", then simply read in all files.
  # Else, read in the inhibitory and activating networks separately.
  result <- NULL
  allFiles <- lapply(fileNames, function(file){
    results <- NULL
    # Read the file.
    tryCatch({
      if(firstColumnIsRowname == TRUE){
        results <- utils::read.csv(paste0(resultDirectory, "/", file),
                                   row.names = 1)
      }else{
        results <- utils::read.csv(paste0(resultDirectory, "/", file))
      }
      # Check that the columns are in the right order.
      if(ncol(results) != 3 || !is.numeric(results[,3])){
        stop("The third column must be numeric!")
      }else{
        # Add row names.
        rownames(results) <- paste(results[,1], results[,2], sep = "_")
        return(results)
      }
      
    }, error = function(cond){
      stop(paste("File has an invalid format:", file))
    })
    
  })
  names(allFiles) <- fileNames
  result <- methods::new("FERRET_Results",results = allFiles, directory = resultDirectory,
                         interpretationOfNegative = interpretationOfNegative)
  return(result)
}

#' Given the results, the source network, and the ingroup and outgroup, this function
#' builds a list of comparisons to make from the source network to the ingroup and
#' from the source network to the outgroup. If inhibitory or activating networks are specified,
#' then these comparisons are segmented out.
#' @param sourceNetwork The name of the file with the source network.
#' @param ingroupToCompare A list of files to compare with the source network (the in-group).
#' @param outgroupToCompare A list of files to compare with the target network (the out-group).
#' @param results An object of type FERRET_Results.
#' @returns A network in the format of an adjacency list
BuildComparisonObject <- function(sourceNetwork, ingroupToCompare, outgroupToCompare,
                                  results){
  # Check for invalid inputs.
  comparisons <- NULL
  if(!is(results, "FERRET_Results") || !is(sourceNetwork, "character") || 
     !is(ingroupToCompare, "character") || !is(outgroupToCompare, "character")){
    stop(paste("SourceNetwork must be a name (string). ingroupToCompare and outgroupToCompare must be",
               "a list of names (strings). results must be an object of type FERRET_Results."))
  }
  if(sourceNetwork %in% ingroupToCompare || sourceNetwork %in% outgroupToCompare){
    stop("The source network cannot be in the in-group or the out-group.")
  }
  if(length(intersect(ingroupToCompare, outgroupToCompare)) > 0){
    stop("The in-group and the out-group cannot overlap.")
  }
  
  # If including all edges, each list includes a comparison from the source to the group.
  ingroup <-lapply(ingroupToCompare, function(net){
    return(methods::new("FERRET_Comparison", source = sourceNetwork, target = net))
  })
  outgroup <-lapply(outgroupToCompare, function(net){
    return(methods::new("FERRET_Comparison", source = sourceNetwork, target = net))
  })
  comparisons <- methods::new("FERRET_Comparisons", ingroup = ingroup, outgroup = outgroup)
}

#' This function computes the AUC score for a given similarity index, where the
#' X axis is the index computed for the out-group and the Y axis is the index
#' computed on the in-group. Rather than the scores being between 0 and 1, we
#' set the minimum to be the minimum similarity across both in-group and out-group.
#' This is because all regulatory networks are expected to have some degree of
#' overlap including but not limited to basic cellular function.
#' @param results An object of type FERRET_Results
#' @param comparisons An object of type FERRET_Comparisons. Note that the 
#' in-group and out-group names in FERRET_Comparisons MUST match the names
#' in FERRET_Results.
#' @param metric The metric to use for computing similarity. One or more elements
#' in the following list:
#' - jaccard: Computes the Jaccard similarity (intersection / union) of edges.
#' - degree: Computes the overlap in degree of each node.
#' - modularity: Computes the differential modularity from a "baseline" network
#' to each "perturbed" network, then scales each modularity value between 0 and 1.
#' - subspace: Computes the geodesic distance between each pair of networks
#' on the Grassmann manifold for a given input k, then scales each
#' distance between 0 and 1.
#' Default is all elements.
#' @param k is the number of eigenvectors to include in the subspace representation
#' of each network. Default is 5.
#' @param numberOfCutoffs The number of cutoffs to include in the AUC/ROC curve.
#' These are determined using the range of values over all networks.
#' @param xlab The label for the X axis of the plot. Default is empty. This only
#' needs to be set when generating the plot.
#' @param ylab The label for the Y axis of the plot. Default is empty. This only
#' needs to be set when generating the plot.
#' @param plotCurve Whether or not to generate a plot. Default is TRUE.
#' @param mode "Percentile" or "Score".
#' @returns An object of type FERRET_ROC_AUC, also generates a  plot as a side effect.
#' @export
ComputeRobustnessAUC <- function(results, comparisons, metric = c("jaccard", "degree", "modularity", "subspace"), 
                        k = 5, numberOfCutoffs = 10, plotCurve = TRUE, xlab = "", ylab = "", mode = "Score"){
  
  # Check the inputs.
  if(!is(results, "FERRET_Results") || !is(comparisons, "FERRET_Comparisons")){
    stop(paste("Inputs must be of type FERRET_Results and FERRET_Comparisons. Construct these using",
                "LoadResults() and BuildComparisonObject(), respectively."))
  }
  
  # If we are plotting, reset the parameters accordingly to generate a panel.
  if(plotCurve == TRUE){
    if(length(metric) == 1){
      par(mfrow = c(1,1))
    }else if(length(metric) == 2){
      par(mfrow = c(1,2))
    }else if(length(metric) >= 3){
      par(mfrow = c(2,2))
    }
  }
  
  # If there are negatives in any of the networks and they do not represent
  # inhibitory edges, scale all of them.
  resultsScaled <- results
  minOfAll <- min(unlist(lapply(results@results, function(network){return(min(network[,3]))})))
  if(minOfAll < 0){
    # Check whether valid entry is included.
    if(length(results@interpretationOfNegative) == 0 || 
       (results@interpretationOfNegative != "inhibitory" && results@interpretationOfNegative != "poor")){
      stop(paste("If negatives are present, valid values for interpretationOfNegative in a FERRET_Results object are",
                 "'inhibitory' and 'poor'"))
    }
    
    # Scale network if appropriate.
    if(results@interpretationOfNegative == "poor" && minOfAll < 0){
      resultsScaled@results <- ScaleNetworks(networks = results@results, min = minOfAll)
    }
    
    # Scale network in both the positive and negative directions so that the scores
    # represent percentiles.
    if(mode == "Percentile"){
      resultsScaled@results <- ScaleNetworksByPercentile(networks = resultsScaled@results,
                                                         numberOfCutoffs = numberOfCutoffs)
    }
    
    # Combine AUC scores if inhibitory.
    auc <- list()
    if(results@interpretationOfNegative == "inhibitory"){
      
      # First compute all ROC scores.
      positives <- GetPositive(resultsScaled)
      negatives <- GetNegative(resultsScaled)
      if(mode == "Percentile"){
        positives@results <- ScaleNetworksByPercentile(networks = positives@results,
                                                           numberOfCutoffs = numberOfCutoffs)
        negatives@results <- ScaleNetworksByPercentile(networks = negatives@results,
                                                       numberOfCutoffs = numberOfCutoffs)
      }
      rocPositives <- ComputeRobustnessForOneEdgeType(results = positives, 
                                                      comparisons = comparisons,
                                                      numberOfCutoffs = numberOfCutoffs,
                                                      metric = metric, 
                                                      k = k, 
                                                      xlab = xlab,
                                                      ylab = ylab,
                                                      plotCurve = FALSE)
      rocNegatives <- ComputeRobustnessForOneEdgeType(results = negatives, 
                                                      comparisons = comparisons, 
                                                      numberOfCutoffs = numberOfCutoffs,
                                                      metric = metric, 
                                                      k = k, 
                                                      xlab = xlab,
                                                      ylab = ylab,
                                                      plotCurve = FALSE)

      # Next, combine them by concatenation.
      rocOverall <- lapply(names(rocPositives@roc), function(metric){
        pos <- rocPositives@roc[[metric]]
        neg <- rocNegatives@roc[[metric]]
        ingroup <- c(pos@ingroup, neg@ingroup)
        ingroupNames <- order(as.numeric(names(ingroup)))
        ingroup <- ingroup[ingroupNames]
        outgroup <- c(pos@outgroup, neg@outgroup)
        outgroupNames <- order(as.numeric(names(outgroup)))
        outgroup <- outgroup[outgroupNames]
        return(methods::new("FERRET_Similarities", ingroup = ingroup, outgroup = outgroup))
      })
      names(rocOverall) <- names(rocPositives@roc)

      # Finally, compute AUC.
      if("jaccard" %in% metric){
        auc[["Jaccard"]] <- AUCTrapezoid(unname(rocOverall$Jaccard@outgroup), unname(rocOverall$Jaccard@ingroup))
        if(plotCurve == TRUE){
          PlotROC(averageSims = rocOverall$Jaccard, auc = auc[["Jaccard"]], xlab = xlab, ylab = ylab,
                  main = "Jaccard Similarity")
        }
      }
      if("degree" %in% metric){
        auc[["Degree"]] <- AUCTrapezoid(unname(rocOverall$Degree@outgroup), unname(rocOverall$Degree@ingroup))
        if(plotCurve == TRUE){
          PlotROC(averageSims = rocOverall$Degree, auc = auc[["Degree"]], xlab = xlab, ylab = ylab,
                  main = "Degree Similarity")
        }
      }
      if("modularity" %in% metric){
        auc[["Modularity"]] <- AUCTrapezoid(unname(rocOverall$Modularity@outgroup), unname(rocOverall$Modularity@ingroup))
        if(plotCurve == TRUE){
          PlotROC(averageSims = rocOverall$Modularity, auc = auc[["Modularity"]], xlab = xlab, ylab = ylab,
                  main = "Modularity Similarity")
        }
      }
      if("subspace" %in% metric){
        auc[["Subspace"]] <- AUCTrapezoid(unname(rocOverall$Subspace@outgroup), unname(rocOverall$Subspace@ingroup))
        if(plotCurve == TRUE){
          PlotROC(averageSims = rocOverall$Subspace, auc = auc[["Subspace"]], xlab = xlab, ylab = ylab,
                  main = "Subspace Similarity")
        }
      }
      auc <-  methods::new("FERRET_ROC_AUC",auc = unlist(auc), roc = rocOverall)
    }
  }
  
  # Compute AUC on only the positive network (which will be the only network if
  # there are no inhibitory edges).
  if(results@interpretationOfNegative != "inhibitory" || minOfAll >= 0){
    if(mode == "Percentile"){
      resultsScaled@results <- ScaleNetworksByPercentile(networks = resultsScaled@results,
                                                         numberOfCutoffs = numberOfCutoffs)
    }
    auc <- ComputeRobustnessForOneEdgeType(results = resultsScaled, 
                                                   comparisons = comparisons, 
                                                   metric = metric, 
                                                   k = k, numberOfCutoffs = numberOfCutoffs,
                                                   plotCurve = plotCurve, xlab = xlab,
                                                   ylab = ylab)
  }
  
  # Reset par.
  par(mfrow = c(1,1))
 
  # Return the AUC scores.
  return(auc)
}

# This function scales networks in the case where a negative represents a poor
# relationship between two genes rather than an inhibitory one.
#' @title ScaleNetworks
#' @param networks A list of unscaled networks
#' @param min Minimum value to consider
#' @returns A list of scaled networks
ScaleNetworks <- function(networks, min){
  
  # Scale.
  allScaledNetworks <- lapply(networks, function(network){
    networkScaled <- network
    networkScaled[,3] <- network[,3] - min
    return(networkScaled)
  })
  return(allScaledNetworks)
}

# This function scales networks so that edges represent percentiles.
#' @title ScaleNetworksByPercentile
#' @param networks A list of unscaled networks
#' @param numberOfCutoffs Number of quantiles to consider
#' @returns A list of scaled networks
ScaleNetworksByPercentile <- function(networks, numberOfCutoffs){

  # Throw an error if the networks are not formatted appropriately.
  if(!is.list(networks) || is.data.frame(networks) || (ncol(networks[[1]]) != 3) || !is.numeric(networks[[1]][,3])){
    stop("networks parameter must be a list of networks with scores in the 3rd column!")
  }
  
  # Throw an error if the number of cutoffs is below 1 or is not an integer.
  if(!is.numeric(numberOfCutoffs) || numberOfCutoffs%%1 > 0 || numberOfCutoffs < 1){
    stop(paste("The number of cutoffs must be an integer and must be above 1. You entered:", numberOfCutoffs))
  }
  
  # Throw an error if the number of unique values in each network is not greater than the
  # number of cutoffs.
  for(network in networks){
    if(length(unique(network[which(network[,3] > 0),3])) > 0 &&
       length(unique(network[which(network[,3] >= 0),3])) < numberOfCutoffs + 1){
      stop("The number of unique positive values in each network must be greater than the number of cutoffs.")
    }
  }
  for(network in networks){
    if(length(unique(network[which(network[,3] < 0),3])) > 0 &&
       length(unique(network[which(network[,3] <= 0),3])) < numberOfCutoffs + 1){
      stop("The number of unique negative values in each network must be greater than the number of cutoffs.")
    }
  }
  
  # Compute where the percentiles should be set based on the number of cutoffs.
  percentilePiece <- 1 / (numberOfCutoffs + 1)
  percentiles <- percentilePiece * seq(0, numberOfCutoffs)
  
  # Scale.
  allScaledNetworks <- lapply(networks, function(network){
    
    # Compute quantiles for network if there are positive edges.
    networkScaled <- network
    if(length(unique(network[which(network[,3] > 0),3])) > 0){
      quantilesPos <- quantile(network[which(network[,3] >= 0),3], percentiles)
  
      # Match scores to quantiles.
      networkScaled[which(networkScaled[,3]>=0),3] <- unlist(lapply(networkScaled[which(networkScaled[,3]>=0),3], function(score){
        
        # Select the minimum quantile less than the score.
        quantileDiff <- quantilesPos - score
        quantilesWithNegativeDiff <- percentiles[which(quantileDiff <= 0)]
        quantileDiffsWhichNegative <- quantileDiff[which(quantileDiff <= 0)]
        matchingQuantile <- quantilesWithNegativeDiff[which.min(abs(quantileDiffsWhichNegative))]
        return(matchingQuantile)
      }))
    }
    
    return(networkScaled)
  })
  return(allScaledNetworks)
}

#' This function consolidates ROC curves and AUC scores to generate a single plot.
#' @param resultList A list of FERRET_ROC_AUC objects to consolidate.
#' @param xlab The label for the X axis of the plot. Default is empty. This only
#' needs to be set when generating the plot.
#' @param ylab The label for the Y axis of the plot. Default is empty. This only
#' needs to be set when generating the plot.
#' @returns The averaged, lowest, and highest values at each cutoff, also generates a  plot as a side effect.
#' @export
ConsolidateRobustness <- function(resultList, xlab, ylab){
  
  # Initialize mean, min, max.
  averageResults <- resultList[[1]]
  minResults <- resultList[[1]]
  maxResults <- resultList[[1]]
  
  # Find absolute minima and maxima.
  absoluteMin <- lapply(names(resultList[[1]]@roc), function(metric){
    allVals <- unlist(lapply(1:length(resultList), function(i){
      return(c(unname(resultList[[i]]@roc[[metric]]@outgroup),
               unname(resultList[[i]]@roc[[metric]]@ingroup)))
    }))
    return(min(allVals[which(!is.na(allVals))]))
  })
  names(absoluteMin) <- names(resultList[[1]]@roc)
  absoluteMax <- lapply(names(resultList[[1]]@roc), function(metric){
    allVals <- unlist(lapply(1:length(resultList), function(i){
      return(c(unname(resultList[[i]]@roc[[metric]]@outgroup),
               unname(resultList[[i]]@roc[[metric]]@ingroup)))
    }))
    return(max(allVals[which(!is.na(allVals))]))
  })
  names(absoluteMax) <- names(resultList[[1]]@roc)

  # Find the average.
  averageResultsRoc <- lapply(names(resultList[[1]]@roc), function(metric){
    averaged <- ConsolidateRobustnessHelper(resultList = resultList, 
                                                          functionToApply = "mean",
                                           metric = metric)
    return(methods::new("FERRET_Similarities", ingroup = averaged[["ingroup"]], 
                        outgroup = averaged[["outgroup"]]))
  })
  names(averageResultsRoc) <- names(resultList[[1]]@roc)
  averageResults@roc <- averageResultsRoc
  averageResultsAUC <- unlist(lapply(names(resultList[[1]]@roc), function(metric){
    return(AUCTrapezoid(x = unname(averageResults@roc[[metric]]@outgroup),
                        y = unname(averageResults@roc[[metric]]@ingroup),
                        absoluteMin  = absoluteMin[[metric]], absoluteMax = absoluteMax[[metric]]))
  }))
  names(averageResultsAUC) <- names(resultList[[1]]@roc)
  averageResults@auc <- averageResultsAUC
  
  # Find the min.
  minResultsRoc <- lapply(names(resultList[[1]]@roc), function(metric){
    minimum <- ConsolidateRobustnessHelper(resultList = resultList, 
                                               functionToApply = "min",
                                           metric = metric)
    return(methods::new("FERRET_Similarities", ingroup = minimum[["ingroup"]], 
                        outgroup = minimum[["outgroup"]]))
  })
  names(minResultsRoc) <- names(resultList[[1]]@roc)
  minResults@roc <- minResultsRoc
  minResultsAUC <- unlist(lapply(names(resultList[[1]]@roc), function(metric){
    return(AUCTrapezoid(unname(minResults@roc[[metric]]@outgroup),
                        unname(minResults@roc[[metric]]@ingroup),
                        absoluteMin  = absoluteMin[[metric]], absoluteMax = absoluteMax[[metric]]))
  }))
  names(minResultsAUC) <- names(resultList[[1]]@roc)
  minResults@auc <- minResultsAUC
  
  # Find the max.
  maxResultsRoc <- lapply(names(resultList[[1]]@roc), function(metric){
    maximum <- ConsolidateRobustnessHelper(resultList = resultList, 
                                           metric = metric,
                                           functionToApply = "max")
    return(methods::new("FERRET_Similarities", ingroup = maximum[["ingroup"]], 
                        outgroup = maximum[["outgroup"]]))
  })
  names(maxResultsRoc) <- names(resultList[[1]]@roc)
  maxResults@roc <- maxResultsRoc
  maxResultsAUC <- unlist(lapply(names(resultList[[1]]@roc), function(metric){
    return(AUCTrapezoid(unname(maxResults@roc[[metric]]@outgroup),
                        unname(maxResults@roc[[metric]]@ingroup),
                        absoluteMin  = absoluteMin[[metric]], absoluteMax = absoluteMax[[metric]]))
  }))
  names(maxResultsAUC) <- names(resultList[[1]]@roc)
  maxResults@auc <- maxResultsAUC
  
  # Plot the ROC curve with the boundaries.
  # Reset the parameters accordingly to generate a panel.
  if(length(resultList[[1]]@roc) == 1){
    par(mfrow = c(1,1))
  }else if(length(resultList[[1]]@roc) == 2){
    par(mfrow = c(1,2))
  }else if(length(resultList[[1]]@roc) >= 3){
    par(mfrow = c(2,2))
  }
  if("Jaccard" %in% names(resultList[[1]]@roc)){
    PlotROC(averageSims = averageResultsRoc$Jaccard, lowestSims = minResultsRoc$Jaccard,
            highestSims = maxResultsRoc$Jaccard, auc = averageResultsAUC[["Jaccard"]], xlab = xlab, ylab = ylab,
            main = "Jaccard Similarity", absoluteMin = absoluteMin[["Jaccard"]], absoluteMax = absoluteMax[["Jaccard"]])
  }
  if("Degree" %in% names(resultList[[1]]@roc)){
    PlotROC(averageSims = averageResultsRoc$Degree, lowestSims = minResultsRoc$Degree,
            highestSims = maxResultsRoc$Degree, auc = averageResultsAUC[["Degree"]], xlab = xlab, ylab = ylab,
            main = "Degree Similarity", absoluteMin = absoluteMin[["Degree"]], absoluteMax = absoluteMax[["Degree"]])
  }
  if("Modularity" %in% names(resultList[[1]]@roc)){
    PlotROC(averageSims = averageResultsRoc$Modularity, lowestSims = minResultsRoc$Modularity,
            highestSims = maxResultsRoc$Modularity, auc = averageResultsAUC[["Modularity"]], xlab = xlab, ylab = ylab,
            main = "Modularity Similarity", absoluteMin = absoluteMin[["Modularity"]], absoluteMax = absoluteMax[["Modularity"]])
  }
  if("Subspace" %in% names(resultList[[1]]@roc)){
    PlotROC(averageSims = averageResultsRoc$Subspace, lowestSims = minResultsRoc$Subspace,
            highestSims = maxResultsRoc$Subspace, auc = averageResultsAUC[["Subspace"]], xlab = xlab, ylab = ylab,
            main = "Subspace Similarity", absoluteMin = absoluteMin[["Subspace"]], absoluteMax = absoluteMax[["Subspace"]])
  }
}

#' This is a helper function for ConsolidateRobustness().
#' @param resultList A list of FERRET_ROC_AUC objects to consolidate.
#' @param functionToApply One of "mean", "min", or "max".
#' @param metric The metric to consider.
#' @returns a vector of new results.
ConsolidateRobustnessHelper <- function(resultList, functionToApply, metric){
  
  # Get range of values for X.
  xValues <- sort(unname(unlist(lapply(1:length(resultList), function(i){
    return(resultList[[i]]@roc[[metric]]@outgroup)
  }))))
  xBins <- seq(min(xValues), max(xValues), (max(xValues) - min(xValues)) / 10)
  
  # For each X bin, get all possible Y values.
  yValues <- lapply(1:length(xBins), function(binIdx){
    return(unlist(lapply(1:length(resultList), function(i){
      xPosInBin <- intersect(which(resultList[[i]]@roc[[metric]]@outgroup >= 0),
                             which(resultList[[i]]@roc[[metric]]@outgroup <= xBins[binIdx]))
      if(binIdx > 1){
        xPosInBin <- intersect(which(resultList[[i]]@roc[[metric]]@outgroup > xBins[binIdx-1]),
                               which(resultList[[i]]@roc[[metric]]@outgroup <= xBins[binIdx]))
      }
      yValsInBin <- unname(resultList[[i]]@roc[[metric]]@ingroup[xPosInBin])
      return(yValsInBin)
    })))
  })

  # Apply the function of interest for each Y bin.
  yConsolidated <- unlist(lapply(yValues, function(y){
    return(get(functionToApply)(y))
  }))
  
  # Remove NA, NAN, and Inf
  whichNotNA <- which(!is.na(yConsolidated))
  yConsolidated <- yConsolidated[whichNotNA]
  xBins <- xBins[whichNotNA]
  whichNotNAN <- which(!is.nan(yConsolidated))
  yConsolidated <- yConsolidated[whichNotNAN]
  xBins <- xBins[whichNotNAN]
  whichNotInf <- which(!is.infinite(yConsolidated))
  yConsolidated <- yConsolidated[whichNotInf]
  xBins <- xBins[whichNotInf]
  
  # Name the vectors.
  names(yConsolidated) <- as.character(1:length(xBins))
  names(xBins) <- as.character(1:length(xBins))
  
  # The results are a list.
  results <- list(ingroup = yConsolidated, outgroup = xBins)
  return(results)
}

#' This is a helper function for ConsolidateRobustness().
#' @param resultList A list of FERRET_ROC_AUC objects to consolidate.
#' @param functionToApply One of "mean", "min", or "max".
#' @param metric The metric to consider.
#' @param group One of "ingroup" or "outgroup"
#' @returns a vector of new results.
ConsolidateRobustnessHelperOld <- function(resultList, functionToApply, metric, group){
  
  # Get result for each cutoff.
  results <- unlist(lapply(1:length(resultList[[1]]@roc[[metric]]@ingroup), function(i){
    val <- NA
    
    # Consolidate across evaluations.
    allAtIndex <- unlist(lapply(resultList, function(results){
      atIndex <- NULL
      if(group == "ingroup"){
        atIndex <- results@roc[[metric]]@ingroup[i]
      }else{
        atIndex <- results@roc[[metric]]@outgroup[i]
      }
      return(atIndex)
    }))
    
    # For those results that are not NA, take the average.
    whichNotNA <- which(!is.na(allAtIndex))
    if(length(whichNotNA) > 0){
      val <- get(functionToApply)(allAtIndex[whichNotNA])
    }
    return(val)
  }))
  names(results) <- paste("cutoff", 1:length(resultList[[1]]@roc[[metric]]@ingroup))
  return(results)
}

#' This reads a file generated by WriteRobustnessAUC().
#' @param fileName A valid file name.
#' @returns An object of type FERRET_ROC_AUC
#' @export
ReadRobustnessAUC <- function(fileName){
  # Read file
  resultDf <- NULL
  tryCatch({
    resultDf <- utils::read.csv(fileName, row.names = 1)
  }, error = function(cond){
    stop(paste("Cannot read file", fileName))
  })
  
  # Extract AUC scores.
  metrics <- setdiff(colnames(resultDf), c("cutoff", "inOrOut"))
  auc <- unlist(resultDf[1,metrics])
  names(auc) <- metrics
  
  # Extract each metric's ROC scores.
  roc <- lapply(metrics, function(metric){
    # Extract the ingroup.
    ingroup <- unlist(resultDf[which(resultDf$inOrOut == "Ingroup"), metric])
    names(ingroup) <- unlist(resultDf[which(resultDf$inOrOut == "Ingroup"), "cutoff"])
    
    # Extract the outgroup.
    outgroup <- unlist(resultDf[which(resultDf$inOrOut == "Outgroup"), metric])
    names(outgroup) <- unlist(resultDf[which(resultDf$inOrOut == "Outgroup"), "cutoff"])
    
    # Create FERRET_Similarities object.
    similarities <- methods::new("FERRET_Similarities", ingroup = ingroup, outgroup = outgroup)
    return(similarities)
  })
  names(roc) <- metrics
  
  # Create a FERRET_ROC_AUC object.
  result <- methods::new("FERRET_ROC_AUC", auc = auc, roc = roc)
}

#' This function saves the results of ComputeRobustnessAUC to a file.
#' @param results An object of type FERRET_ROC_AUC
#' @param fileName A valid file name.
#' @returns NULL
#' @export
WriteRobustnessAUC <- function(results, fileName){
  # If arguments are of incorrect types, stop.
  if(!is(results, "FERRET_ROC_AUC")){
    stop("Error: Input must be an object of type FERRET_ROC_AUC")
  }
  # Add a row representing AUC.
  resultDfAuc <- as.data.frame(t(as.data.frame(results@auc)))
  colnames(resultDfAuc) <- names(results@auc)
  
  # Add a row for each of the ingroup ROC values.
  resultDfIngroup <- do.call(cbind, lapply(1:length(results@roc), function(i){
    df <- data.frame(V1 = results@roc[[i]]@ingroup)
  }))
  colnames(resultDfIngroup) <- names(results@roc)
  
  # Add a row for each of the outgroup ROC values.
  resultDfOutgroup <- do.call(cbind, lapply(1:length(results@roc), function(i){
    df <- data.frame(V1 = results@roc[[i]]@outgroup)
  }))
  colnames(resultDfOutgroup) <- names(results@roc)
  
  # Combine all values.
  resultDf <- do.call(rbind, list(resultDfAuc, resultDfIngroup, resultDfOutgroup))
  resultDf$cutoff <- c("NA", names(results@roc[[1]]@ingroup), names(results@roc[[1]]@outgroup))
  resultDf$inOrOut <- c("NA", rep("Ingroup", length(results@roc[[1]]@ingroup)),
                        rep("Outgroup", length(results@roc[[1]]@outgroup)))
  rownames(resultDf) <- c("AUC", paste0("IngroupROC_", names(results@roc[[1]]@ingroup)),
                          paste0("OutgroupROC_", names(results@roc[[1]]@outgroup)))
  
  # Save in file.
  tryCatch({
    utils::write.csv(resultDf, fileName)
  }, error = function(cond){
    stop(paste("Cannot write to file", fileName))
  })
  
}

#' This function finds the differential network between the common networks for
#' two different scenarios (e.g., the ingroup and outgroup). We note that, in this
#' scenario, we only report edges where the magnitude in the ingroup network is larger
#' than in the outgroup network and both edges are in the same direction. 
#' That is, negative means that a more significant inhibitory edge is present in the
#' ingroup network than in the outgroup network. It does NOT mean that
#' the magnitude is greater in the outgroup network than in the ingroup network.
#' @param results An object of type FERRET_Results
#' @param comparisons A list of FERRET_Comparisons object.
#' @returns A single network in the same format as the result networks.
#' @export
GetDifferentialNetworks <- function(results, comparisons){
  # First, compute the overlap between all in-group networks.
  network1 <- GetCommonNetworkAcrossComparisons(results = results, comparisons = comparisons@ingroup)

  # Then, compute the overlap between all out-group networks.
  network2 <- GetCommonNetworkAcrossComparisons(results = results, comparisons = comparisons@outgroup,
                                                useTargetsOnly = TRUE)

  # Expand the networks.
  network1Expanded <- network1
  network2Expanded <- network2
  network1EdgesOnly <- setdiff(rownames(network1), rownames(network2))
  network2EdgesOnly <- setdiff(rownames(network2), rownames(network1))
  network1Expanded[network2EdgesOnly,] <- network2[network2EdgesOnly,]
  network1Expanded[network2EdgesOnly, 3] <- 0
  network2Expanded[network1EdgesOnly,] <- network1[network1EdgesOnly,]
  network2Expanded[network1EdgesOnly, 3] <- 0
  network1Expanded <- network1Expanded[order(rownames(network1Expanded)),]
  network2Expanded <- network2Expanded[order(rownames(network2Expanded)),]

  # Find the differential network.
  bothNegativeAndIngroupMagnitudeGreater <- Reduce(intersect, list(which(network1Expanded[,3] < 0), 
                                                                    which(network2Expanded[,3] < 0),
                                                                   which(network1Expanded[,3] <= network2Expanded[,3])))
  bothPositiveAndIngroupMagnitudeGreater <- Reduce(intersect, list(which(network1Expanded[,3] > 0), 
                                                                   which(network2Expanded[,3] > 0),
                                                                   which(network1Expanded[,3] >= network2Expanded[,3])))
  whichTakeDifference <- c(bothNegativeAndIngroupMagnitudeGreater, bothPositiveAndIngroupMagnitudeGreater)
  ingroupNegativeOutgroupPositive <- intersect(which(network1Expanded[,3] < 0), which(network2Expanded[,3] >= 0))
  outgroupNegativeIngroupPositive <- intersect(which(network1Expanded[,3] > 0), which(network2Expanded[,3] <= 0))
  whichKeepOnlyIngroup <- c(ingroupNegativeOutgroupPositive, outgroupNegativeIngroupPositive)
  whichZero <- setdiff(1:nrow(network1Expanded), c(whichTakeDifference, whichKeepOnlyIngroup))
  
  network1Diff <- network1Expanded
  network1Diff[whichTakeDifference,3] <- network1Expanded[whichTakeDifference,3] - network2Expanded[whichTakeDifference,3]
  network1Diff[whichKeepOnlyIngroup,3] <- network1Expanded[whichKeepOnlyIngroup,3]
  network1Diff[whichZero,3] <- 0
  
  # Remove zero edges.
  network1Diff <- network1Diff[which(network1Diff[,3] != 0),]
  return(network1Diff)
}


#' This function finds the common network between the common networks for two
#' different scenarios (e.g., the ingroup and the outgroup).
#' @param results An object of type FERRET_Results
#' @param comparisons A FERRET_Comparisons object.
#' @returns A single network in the same format as the result networks.
#' @export
GetCommonNetwork <- function(results, comparisons){
  # First, compute the overlap between all in-group networks.
  network1 <- GetCommonNetworkAcrossComparisons(results = results, comparisons = comparisons@ingroup)
  
  # Then, compute the overlap between all out-group networks.
  network2 <- GetCommonNetworkAcrossComparisons(results = results, comparisons = comparisons@outgroup,
                                                useTargetsOnly = TRUE)
  
  # Intersect the networks.
  shared <- intersect(rownames(network1), rownames(network2))
  network1Shared <- network1[shared,]
  network2Shared <- network2[shared,]
  
  # Remove edges with opposite signs.
  sameSignEdges <- which(network1Shared[,3] * network2Shared[,3] > 0)
  network1Shared <- network1Shared[sameSignEdges,]
  network2Shared <- network2Shared[sameSignEdges,]
  
  # Add the score with the minimum magnitude.
  network1Shared[,3] <- pmin(abs(network1Shared[,3]), abs(network2Shared[,3]))
  network1Shared[which(network2Shared[,3] < 0),3] <- -1 * network1Shared[which(network2Shared[,3] < 0),3]
  
  return(network1Shared)
}

#' This function computes the percent weight of the input network that overlaps with the prior.
#' @param network The input network, where the transcription factors in the first column 
#' are gene symbols.
#' @param prior The prior
#' @returns The Jaccard index (scalar)
#' @export
GetOverlapWithPrior <- function(network, prior){
  # Check for empty prior.
  if(nrow(prior) == 0){
    stop("Prior cannot be empty!")
  }
  
  # Compute overlap.
  netEdges <- paste(network[,1], network[,2], sep = "_")
  rownames(network) <- netEdges
  priorEdges <- paste(prior[,1], prior[,2], sep = "_")
  sharedEdges <- intersect(netEdges, priorEdges)
  
  # Compute return value.
  sumShared <- sum(abs(network[sharedEdges, 3]))
  sumTotal <- sum(abs(network[,3]))
  retVal <- 0
  if(sumTotal > 0){
    retVal <- sumShared / sumTotal
  }
  return(retVal)
}

#' This function runs pathway analysis for the network targets.
#' @param network The input network, where the target genes in the second column
#' are gene symbols.
#' @param pathways A pathway object, generated using gmtPathways() from
#' the fgsea R package. The format of this object is a named list of 
#' vectors, where the elements of each vector are gene symbols.
#' @returns The pathway results formatted as a data frame.
#' @export
RunPathwayAnalysis <- function(network, pathways){
  
  # Obtain the targeting score for each gene.
  uniqueTargets <- unique(network[,2])
  differentialTargeting <- unlist(lapply(uniqueTargets, function(gene){
    return(sum(network[which(network[,2] == gene), 3]))
  }))
  names(differentialTargeting) <- uniqueTargets
  
  # Select only the genes with positive differential targeting.
  positiveDifferentialTargeting <- differentialTargeting[which(differentialTargeting > 0)]
  
  # Run pathway analysis.
  pathwayResult <- fgsea::fgsea(pathways = pathways, stats = positiveDifferentialTargeting,
                                scoreType = "pos")
  
  # Compile the "leading edge" vector into a list.
  leadingEdge <- unlist(lapply(1:length(pathwayResult$leadingEdge), function(i){
    return(paste(pathwayResult$leadingEdge[i][[1]], collapse = "; "))
  }))
  pathwayResultDf <- data.frame(pathway = pathwayResult$pathway,
                                pval = pathwayResult$pval,
                                padj = pathwayResult$padj,
                                pvalErrorSD = pathwayResult$log2err,
                                enrichmentScore = pathwayResult$ES,
                                normalizedEnrichmentScore = pathwayResult$NES,
                                remainingGeneCount = pathwayResult$size,
                                leadingGenesDrivingEnrichment = leadingEdge)
  
  # Return result.
  return(pathwayResultDf)
}

#' This function finds the network common to all folds and the original network.
#' the weight of each edge is the minimum magnitude weight over all networks.
#' @param results An object of type FERRET_Results
#' @param comparisons A list of FERRET_Comparisons objects. 
#' Note that the in-group and out-group names in FERRET_Comparisons MUST match the names
#' in FERRET_Results.
#' @param useTargetsOnly Whether to include the source network in the evaluation.
#' @returns A single network in the same format as the result networks.
GetCommonNetworkAcrossComparisons <- function(results, comparisons, useTargetsOnly = FALSE){
  # Flatten the comparisons list.
  flatComparisons <- unique(unlist(lapply(comparisons, function(comparison){
    retVal <- c(comparison@source, comparison@target)
    if(useTargetsOnly == TRUE){
      retVal <- comparison@target
    }
    return(retVal)
  })))

  # Extract the subset of networks in the comparisons list.
  networks <- results@results[which(names(results@results) %in% flatComparisons)]
  
  # If interpretation of negative is "poor", rescale the networks.
  minOfAll <- min(unlist(lapply(results@results, 
                                function(network){return(min(network[,3]))})))
  if(results@interpretationOfNegative == "poor" && minOfAll < 0){
    networks <- ScaleNetworks(networks = networks, min = minOfAll)
  }

  # Intersect all networks.
  # Initialize.
  runningNetwork <- networks[[1]]
  rownames(runningNetwork) <- paste(runningNetwork[,1], runningNetwork[,2], sep = "_")
  # Loop through all
  if(length(networks) > 1){
    for(network in networks[2:length(networks)]){
      # Extract the overlap between both the positives and negatives.
      
      # Intersect edges.
      rownames(network) <- paste(network[,1], network[,2], sep = "_")
      sharedEdges <- intersect(rownames(runningNetwork), rownames(network))
      runningNetwork <- runningNetwork[sharedEdges,]
      network <- network[sharedEdges,]
      
      # Remove edges with opposite signs.
      sameSignEdges <- which(network[,3] * runningNetwork[,3] > 0)
      runningNetwork <- runningNetwork[sameSignEdges,]
      network <- network[sameSignEdges,]
      
      # Add the score with the minimum magnitude.
      runningNetwork[,3] <- pmin(abs(network[,3]), abs(runningNetwork[,3]))
      runningNetwork[which(network[,3] < 0),3] <- -1 * runningNetwork[which(network[,3] < 0),3]
    }
  }
  
  # Return final network.
  return(runningNetwork)
}

#' This function takes in a network of associations between PCs and reconstructs
#' a network of associations between genes using the PC loadings.
#' @param pcNetwork The network of associations between PCs, formatted as a list with
#' source, target, and score info.
#' @param PC The principal components, formatted as a matrix with rows equal to the number
#' of genes and columns equal to the PCs.
#' @param zeroCutoff The zero cutoff.
#' @returns A single network in the same format as the input network.
#' @export
GetFullNetworkFromPCANetwork <- function(pcNetwork, PC, zeroCutoff = 0.001){
  
  # Create an adjacency matrix for the PC network.
  pcNetworkRenamed <- pcNetwork
  colnames(pcNetworkRenamed)[3] <- "Weight"
  pcsNotInNetwork <- setdiff(colnames(PC), c(pcNetworkRenamed[,1], pcNetworkRenamed[,2]))
  if(length(pcsNotInNetwork) > 0){
    expansion <- data.frame(source = pcsNotInNetwork,
                            target = rep(pcNetwork$target[1], length(pcsNotInNetwork)),
                            Weight = rep(0, length(pcsNotInNetwork)))
    rownames(expansion) <- paste(expansion$source, expansion$target, sep = "_")
    colnames(expansion) <- colnames(pcNetworkRenamed)
    pcNetworkRenamed <- rbind(pcNetworkRenamed, expansion)
  }
  
  pcAdjacency <- as.matrix(igraph::get.adjacency(igraph::graph.data.frame(pcNetworkRenamed), attr = "Weight"))
  pcAdjacency <- pcAdjacency[colnames(PC), colnames(PC)]
  
  # Use matrix multiplication to obtain an edge-edge network.
  sourceTargetMapping <- as.matrix(PC) %*% pcAdjacency %*% t(PC)

  # Melt the source-target mapping to an adjacency list.
  sourceTargetMappingAdj <- reshape2::melt(sourceTargetMapping)

  # Remove zeros and self-loops.
  sourceTargetMappingAdj <- sourceTargetMappingAdj[which(abs(sourceTargetMappingAdj[,3]) > zeroCutoff),]
  sourceTargetMappingAdj <- sourceTargetMappingAdj[which(sourceTargetMappingAdj[,1]
                                                         != sourceTargetMappingAdj[,2]),]

  # Set names
  rearrangedSourceTargetMapping <- sourceTargetMappingAdj
  colnames(rearrangedSourceTargetMapping) <- colnames(pcNetwork)
  rownames(rearrangedSourceTargetMapping) <- paste(rearrangedSourceTargetMapping[,1],
                                                   rearrangedSourceTargetMapping[,2], sep = "_")
  rearrangedSourceTargetMapping[,1] <- as.character(rearrangedSourceTargetMapping[,1])
  rearrangedSourceTargetMapping[,2] <- as.character(rearrangedSourceTargetMapping[,2])

  # Return
  return(rearrangedSourceTargetMapping)
}

#' Helper function that computes the AUC score for either all edges or
#' for activating and inhibitory edges separately.
#' @param results An object of type FERRET_Results
#' @param comparisons An object of type FERRET_Comparisons. Note that the 
#' in-group and out-group names in FERRET_Comparisons MUST match the names
#' in FERRET_Results.
#' @param metric The metric to use for computing similarity. One or more elements
#' in the following list:
#' - jaccard: Computes the Jaccard similarity (intersection / union) of edges.
#' - degree: Computes the overlap in degree of each node.
#' - modularity: Computes the differential modularity from a "baseline" network
#' to each "perturbed" network, then scales each modularity value between 0 and 1.
#' - subspace: Computes the geodesic distance between each pair of networks
#' on the Grassmann manifold for a given input k, then scales each
#' distance between 0 and 1.
#' Default is all elements.
#' @param k is the number of eigenvectors to include in the subspace representation
#' of each network. Default is 5.
#' @param xlab The label for the X axis of the plot. Default is empty. This only
#' needs to be set when generating the plot.
#' @param ylab The label for the Y axis of the plot. Default is empty. This only
#' needs to be set when generating the plot.
#' @param plotCurve Whether or not to generate a plot. Default is TRUE.
#' @param mode "Percentile" or "Score".
#' @param numberOfCutoffs Number of cutoffs to evaluate.
#' @returns An AUC score, also generates a  plot as a side effect.
#' @export
ComputeRobustnessForOneEdgeType <- function(results, comparisons, metric = c("jaccard", "degree", "modularity", "subspace"), 
                                 k = 5, plotCurve = TRUE, xlab = "", ylab = "", mode = "Score", numberOfCutoffs = 10){
  
  # Stop if no metrics were input.
  possibleMetrics <- c("jaccard", "modularity", "subspace", "degree")
  if(length(metric) == 0 || !is.character(metric) || length(setdiff(metric, possibleMetrics)) > 0){
    stop(paste("Metric must be one or more of the following:", paste(possibleMetrics, collapse = ", ")))
  }
  
  # Compute cutoffs for ROC/AUC using all results. Since these are already mapped
  # to quantiles, we only need to look at the unique values.
  cutoffs <- ObtainNetworkCutoffs(results, numberOfCutoffs)
  if(mode == "Percentile"){
    cutoffs <- seq(1, numberOfCutoffs, 1)
    uniqueVals <- unique(results@results[[1]][,3])
    if(length(uniqueVals) > 1){
      cutoffs <- sort(unique(results@results[[1]][,3]))
      cutoffs <- cutoffs[which(abs(cutoffs) > 0)]
    }
  }

  # For each metric, compute the similarities and the AUC scores.
  auc <- list()
  roc <- list()
  if("jaccard" %in% metric){
    jaccardSim <- ComputeSimilarities(results = results, comparisons = comparisons,
                                     cutoffs = cutoffs, metric = "jaccard")
    message("Finished all Jaccard similarity metrics")
    
    aucJaccard <- 0
    if(length(which(is.na(jaccardSim@ingroup))) < length(jaccardSim@ingroup) &&
       length(which(is.na(jaccardSim@outgroup))) < length(jaccardSim@outgroup)){
      aucJaccard <- AUCTrapezoid(unname(jaccardSim@outgroup), unname(jaccardSim@ingroup))
    }
    auc[["Jaccard"]] <- aucJaccard
    roc[["Jaccard"]] <- jaccardSim
    if(plotCurve == TRUE){
      PlotROC(averageSims = jaccardSim, auc = aucJaccard, xlab = xlab, ylab = ylab,
              main = "Jaccard Similarity")
    }
  }
  if("degree" %in% metric){
    degreeSim <- ComputeSimilarities(results = results, comparisons = comparisons,
                                   cutoffs = cutoffs, metric = "degree")
    message("Finished all Degree similarity metrics")
    aucDegree <- 0
    if(length(which(is.na(degreeSim@ingroup))) < length(degreeSim@ingroup) &&
       length(which(is.na(degreeSim@outgroup))) < length(degreeSim@outgroup)){
      aucDegree <- AUCTrapezoid(unname(degreeSim@outgroup), unname(degreeSim@ingroup))
    }
    auc[["Degree"]] <- aucDegree
    roc[["Degree"]] <- degreeSim
    if(plotCurve == TRUE){
      PlotROC(averageSims = degreeSim, auc = aucDegree, xlab = xlab, ylab = ylab,
              main = "Degree Similarity")
    }
  }
  if("modularity" %in% metric){
    modularitySim <- ComputeSimilarities(results = results, comparisons = comparisons,
                                           cutoffs = cutoffs, metric = "modularity")
    message("Finished all Modularity similarity metrics")
    aucModularity <- 0
    if(length(which(is.na(modularitySim@ingroup))) < length(modularitySim@ingroup) &&
       length(which(is.na(modularitySim@outgroup))) < length(modularitySim@outgroup)){
      aucModularity <- AUCTrapezoid(unname(modularitySim@outgroup), unname(modularitySim@ingroup))
    }
    auc[["Modularity"]] <- aucModularity
    roc[["Modularity"]] <- modularitySim
    if(plotCurve == TRUE){
      PlotROC(averageSims = modularitySim, auc = aucModularity, xlab = xlab, ylab = ylab,
              main = "Modularity Similarity")
    }
  }
  if("subspace" %in% metric){
    subspaceSim <- ComputeSimilarities(results = results, comparisons = comparisons, k = k,
                                       cutoffs = cutoffs, metric = "subspace")
    message("Finished all Subspace similarity metrics")
    aucSubspace <- 0
    if(length(which(is.na(subspaceSim@ingroup))) < length(subspaceSim@ingroup) &&
       length(which(is.na(subspaceSim@outgroup))) < length(subspaceSim@outgroup)){
      aucSubspace <- AUCTrapezoid(unname(subspaceSim@outgroup), unname(subspaceSim@ingroup))
    }
    auc[["Subspace"]] <- aucSubspace
    roc[["Subspace"]] <- subspaceSim
    if(plotCurve == TRUE){
      PlotROC(averageSims = subspaceSim, auc = aucSubspace, xlab = xlab, ylab = ylab,
              main = "Subspace Similarity")
    }
  }
  retVal <- methods::new("FERRET_ROC_AUC",auc = unlist(auc), roc = roc)

  # Return the AUC scores.
  return(retVal)
}

#' Given the values for axis X and axis Y, this function computes the AUC score
#' with linear interpolation between the points. This is a helper function for
#' ComputeAUC.
#' @param x The value on the X axis
#' @param y The value on the Y axis
#' @param absoluteMin The absolute minimum value. If NULL, it is the minimum
#' value across both X and Y.
#' @param absoluteMax The absolute maximum value. If NULL, it is the maximum
#' value across both X and Y.
#' @returns The AUC score
AUCTrapezoid <- function(x, y, absoluteMin = NULL, absoluteMax = NULL) {

  # Get minima and maxima.
  whichNotNA <- intersect(which(!is.na(x)), which(!is.na(y)))
  x <- x[whichNotNA]
  y <- y[whichNotNA]
  minX <- min(x)
  minY <- min(y)
  maxX <- max(x)
  maxY <- max(y)

  # Set minima and maxima.
  zeroCutoff <- 0.0000000001
  
  # If there are negatives, stop.
  if(minX < (-1 * zeroCutoff) || minY < (-1 * zeroCutoff)){
    stop("ERROR: You have input negative similarities.")
  }
  
  # If there is no variation in either x or y, return NA.
  auc <- NA
  if(minX != maxX || minY != maxY){
    
    # Sort x.
    orderX <- order(x)
    x <- x[orderX]
    y <- y[orderX]

    # When there are duplicate values of X, order according to Y within the duplicates.
    duplicates <- names(table(x))[which(table(x) > 1)]
    for(duplicate in duplicates){
      whichXDuplicate <- which(x == as.numeric(duplicate))
      y[whichXDuplicate] <- y[whichXDuplicate][order(y[whichXDuplicate])]
    }
    
    # Obtain global minima and maxima.
    min_total <- min(minY, minX)
    max_total <- max(maxY, maxX)
    
    # Set min and max.
    if(!is.null(absoluteMin)){
      min_total <- absoluteMin
    }
    if(!is.null(absoluteMax)){
      max_total <- absoluteMax
    }

    # Compute the AUC.
    auc <- 0
    for(i in 2:length(x)){
      # Add the rectangular portion.
      auc <- auc + (min(y[i], y[i-1]) - min_total) * (x[i] - x[i-1])

      # Add the triangular portion.
      auc <- auc + (max(y[i], y[i-1]) - min(y[i], y[i-1])) * (x[i] - x[i-1]) / 2
    }
    
    # Add the extrapolated portion.
    auc <- auc + (max_total - maxX) * (y[length(y)] - min_total)

    # Divide by the portion covered.
    auc <- auc / ((max_total - min_total) * (max_total - min_total))
  }
  return(auc)
}

#' This function plots the ROC curve for a given similarity index.
#' @param averageSims The average similarities for the in-groups and out-groups.
#' This is a FERRET_Similarities object.
#' @param lowestSims The lower boundary for the curve. This is a FERRET_Similarities object.
#' @param highestSims The upper boundary for the curve. This is a FERRET_Similarities object.
#' @param auc This is a total AUC score to include on the curve as text.
#' @param xlab This is the label to include on the X axis.
#' @param ylab This is the label to include on the Y axis.
#' @param main This is the title of the plot.
#' @param absoluteMin The absolute minimum similarity in the results.
#' @param absoluteMax The absolute maximum similarity in the results.
#' @returns Nothing, but generates a plot as a side effect.
#' @export
PlotROC <- function(averageSims, lowestSims = NULL, highestSims = NULL, auc, 
                    xlab, ylab, main, absoluteMin = NULL, absoluteMax = NULL){
  
  # Check for correct input.
  if(!is(averageSims, "FERRET_Similarities")){
    stop("The averageSims parameter must be of type FERRET_Similarities.")
  }
  
  # Check for NA values.
  if(length(which(is.na(averageSims@ingroup))) == length(averageSims@ingroup) ||
     length(which(is.na(averageSims@outgroup))) == length(averageSims@outgroup)){
    stop("Cannot plot ROC curve with all NA values.")
  }
  
  # Get X and Y.
  whichNotNA <- intersect(which(!is.na(averageSims@ingroup)), which(!is.na(averageSims@outgroup)))
  y <- averageSims@ingroup[whichNotNA]
  x <- averageSims@outgroup[whichNotNA]
  
  # Check for negatives.
  if(min(c(y, x)) < 0){
    stop("ERROR: You have input negative similarities.")
  }
  
  # Compute the minimum and maximum values.
  min_total <- min(c(y, x))
  max_total <- max(c(y, x))
  
  # Order the values on the X and Y axis according to X.
  orderX <- order(x)
  x <- x[orderX]
  y <- y[orderX]
  
  # When there are duplicate values of X, order according to Y within the duplicates.
  duplicates <- names(table(x))[which(table(x) > 1)]
  for(duplicate in duplicates){
    whichXDuplicate <- which(x == as.numeric(duplicate))
    y[whichXDuplicate] <- y[whichXDuplicate][order(y[whichXDuplicate])]
  }
  
  # If we have upper and lower boundaries, use these to set the min and max.
  limMin <- min_total
  limMax <- max_total
  if(!is.null(lowestSims)){
    whichNotNALow <- intersect(which(!is.na(lowestSims@ingroup)), which(!is.na(lowestSims@outgroup)))
    limMin <- min(c(lowestSims@ingroup[whichNotNALow], lowestSims@outgroup[whichNotNALow]))
  }
  if(!is.null(highestSims)){
    whichNotNAHigh <- intersect(which(!is.na(highestSims@ingroup)), which(!is.na(highestSims@outgroup)))
    limMax <- max(c(highestSims@ingroup[whichNotNAHigh], highestSims@outgroup[whichNotNAHigh]))
  }
  if(!is.null(absoluteMin)){
    limMin <- absoluteMin
  }
  if(!is.null(absoluteMax)){
    limMax <- absoluteMax
  }
  graphics::plot(x, y, xlab = xlab, ylab = ylab, type = "n", main = main,
       xlim = c(limMin, limMax), ylim = c(limMin, limMax))
  
  # Add the curve.
  graphics::lines(x, y)
  
  # Add the extrapolated curve.
  graphics::lines(x = c(x[length(x)], limMax), y = c(y[length(x)], y[length(x)]),
                 lty = "dotted")
  
  # Add the lower boundary.
  if(!is.null(lowestSims)){
    whichNotNALower <- intersect(which(!is.na(lowestSims@ingroup)), which(!is.na(lowestSims@outgroup)))
    yLowest <- lowestSims@ingroup[whichNotNA]
    xLowest <- lowestSims@outgroup[whichNotNA]
    orderXLowest <- order(xLowest)
    xLowest <- xLowest[orderXLowest]
    yLowest <- yLowest[orderXLowest]
    graphics::lines(xLowest, yLowest, lty = "dashed")
  }
  
  # Add the upper boundary.
  if(!is.null(highestSims)){
    whichNotNAHigher <- intersect(which(!is.na(highestSims@ingroup)), which(!is.na(highestSims@outgroup)))
    yHighest <- highestSims@ingroup[whichNotNA]
    xHighest <- highestSims@outgroup[whichNotNA]
    orderXHighest <- order(xHighest)
    xHighest <- xHighest[orderXHighest]
    yHighest <- yHighest[orderXHighest]
    graphics::lines(xHighest, yHighest, lty = "dashed")
  }
  
  # Include the text.
  graphics::text(paste(paste0("RAUC:"), format(round(auc, 2), nsmall = 2)), 
       x = 0.8 * (max_total - min_total) + min_total, 
       y = 0.1 * (max_total - min_total) + min_total)
}

#' This function obtains the cutoffs to use in computing the AUC scores.
#' @param results Object of type FERRET_Results
#' @param numberOfCutoffs The number of points to include in calculating the AUC score.
#' @returns A vector of cutoffs
ObtainNetworkCutoffs <- function(results, numberOfCutoffs){
  # Throw an error if the results are not of the correct type.
  if(!is(results, "FERRET_Results")){
    stop("Results must be of type FERRET_Results.")
  }
  
  # Throw an error if the number of cutoffs is below 1 or is not an integer.
  if(!is.numeric(numberOfCutoffs) || numberOfCutoffs%%1 > 0 || numberOfCutoffs < 1){
    stop(paste("The number of cutoffs must be an integer and must be above 1. You entered:", numberOfCutoffs))
  }
  
  # Obtain full range of values across all networks.
  allWeights <- unlist(lapply(results@results, function(network){return(as.numeric(network[,3]))}))
  minVal <- min(allWeights)
  maxVal <- max(allWeights)

  # Split this value into a given number of cutoffs.
  # If all networks are empty, set cutoffs to be 1-n. They will not be used.
  cutoffs <- seq(1, numberOfCutoffs, 1)
  
  # Throw an error if there are negative values.
  if(minVal < 0){
    stop("Negative scores must be adjusted before calculating cutoffs.")
  }
  
  if(length(allWeights) > 0){
    cutoffs <- seq(minVal, maxVal, by = ((maxVal - minVal) / numberOfCutoffs))
  }
  
  # Return the cutoffs.
  return(cutoffs)
}

#' This function returns each network filtered to only include those edges with
#' scores above a specified cutoff.
#' @param network The network to filter.
#' @param cutoff The cutoff to use when filtering the network.
#' @returns A network in the format of an adjacency list.
GetEdgesAboveCutoff <- function(network, cutoff){
  # Include a tolerance parameter in case of rounding errors.
  tolerance <- 0.0000000001
  return(network[which(network[,3] >= cutoff - tolerance),])
}

#' Helper function to compute the similarities between the baseline network,
#' the ingroup, and the outgroup.
#' @param results Object of type FERRET_Results
#' @param comparisons Object of type FERRET_Comparisons
#' @param cutoffs The list of cutoffs to use when filtering the networks
#' @param metric The metric to use for computing similarity. An element
#' in the following list:
#' - jaccard: Computes the Jaccard similarity (intersection / union) of edges.
#' - degree: Computes the overlap in degree of each node.
#' - modularity: Computes the differential modularity from a "baseline" network
#' to each "perturbed" network, then scales each modularity value between 0 and 1.
#' - subspace: Computes the geodesic distance between each pair of networks
#' on the Grassmann manifold for a given input k, then scales each
#' distance between 0 and 1.
#' Default is all elements.
#' @param k is the number of eigenvectors to include in the subspace representation
#' of each network. Default is 5.
#' @returns Object of type FERRET_Similarity
ComputeSimilarities <- function(results, comparisons, cutoffs, metric, k = 5){
  
  # Check that input types are correct.
  if(!is(results, "FERRET_Results") || !is(comparisons, "FERRET_Comparisons")){
    stop(paste("Inputs must be of type FERRET_Results and FERRET_Comparisons. Construct these using",
               "LoadResults() and BuildComparisonObject(), respectively."))
  }
  if(!is.numeric(cutoffs)){
    stop("You must specify at least one cutoff of numeric type!")
  }
  
  # Set a cutoff for zero.
  zeroCutoff <- 0.0000000001
   
  # For each cutoff, do all comparisons for in-group and out-group.
  sim <- lapply(cutoffs, function(cutoff){
    
    # Obtain a new FERRET_Results object where everything is filtered by cutoff.
    cutoffResults <- results
    cutoffResults@results <- lapply(results@results, function(result){
      return(GetEdgesAboveCutoff(result, cutoff))
    })

    # If we are doing subspace similarity, first compute the projections of all networks.
    projections <- NULL
    if(metric == "subspace"){
      # Compute the projection of the identity matrix.
      uniqueNodes <- sort(unique(unlist(lapply(cutoffResults@results, function(network){
        return(c(network[,1], network[,2]))
      }))))
      if(length(uniqueNodes) <= k){
        uniqueNodes <- c(uniqueNodes, as.character(seq(1, k+2, 1)))
      }
      
      # Compute the projections for the networks.
      projections <- lapply(cutoffResults@results, function(network){
        projection <- ComputeProjection(network = network, k, allNodes = uniqueNodes)
        projection <- do.call(cbind, lapply(colnames(projection), function(c){
          projectionC <- projection[,c]
          projectionC[which(abs(projectionC) < zeroCutoff)] <- 0
          return(projectionC)
        }))
      })
      names(projections) <- names(cutoffResults@results)
    }

    # Compute all ingroup similarities.
    ingroupSim <- ComputeSimilarityForGroup(results = cutoffResults, comparisons = comparisons@ingroup,
                                            metric = metric, k = k,
                                            projections = projections)
    ingroupSim[which(abs(ingroupSim) < zeroCutoff)] <- 0
    outgroupSim <- ComputeSimilarityForGroup(results = cutoffResults, comparisons = comparisons@outgroup,
                                             metric = metric, k = k, projections = projections)
    outgroupSim[which(abs(outgroupSim) < zeroCutoff)] <- 0

    # Compute the average similarity.
    return(list(ingroupSim = ingroupSim, outgroupSim = outgroupSim))
  })

  # Rearrange the returned object.
  simIngroup <- unlist(lapply(1:length(cutoffs), function(i){
    return(sim[[i]][["ingroupSim"]])
  }))
  simOutgroup <- unlist(lapply(1:length(cutoffs), function(i){
    return(sim[[i]][["outgroupSim"]])
  }))
  names(simIngroup) <- cutoffs
  names(simOutgroup) <- cutoffs

  # Return the object.
  return(methods::new("FERRET_Similarities",ingroup = simIngroup, outgroup = simOutgroup))
}

#' Helper function to return a subnetwork with only negative edges.
#' @param results A FERRET_Results object.
#' @returns A new adjacency list with only the negative scores, converted to positive.
GetNegative <- function(results){
  negResults <- results
  for(name in names(results@results)){
    network <- results@results[[name]]
    networkNeg <- network[which(network[,3] < 0),]
    networkNeg[,3] <- -1 * networkNeg[,3]
    negResults@results[[name]] <- networkNeg
  }
  
  return(negResults)
}

#' Helper function to return a subnetwork with only positive edges.
#' @param results A FERRET_Results object.
#' @returns A new adjacency list with only the positive scores.
GetPositive <- function(results){
  posResults <- results
  for(name in names(results@results)){
    network <- results@results[[name]]
    networkPos <- network[which(network[,3] >= 0),]
    posResults@results[[name]] <- networkPos
  }
  
  return(posResults)
}

#' Helper function to compute the similarities within the ingroup or the outgroup
#' @param results Object of type FERRET_Results
#' @param comparisons The ingroup or outgroup slot from an object of type
#' FERRET_Comparisons.
#' @param metric The metric to use for computing similarity. An element
#' in the following list:
#' - jaccard: Computes the Jaccard similarity (intersection / union) of edges.
#' - degree: Computes the overlap in degree of each node.
#' - modularity: Computes the differential modularity from a "baseline" network
#' to each "perturbed" network, then scales each modularity value between 0 and 1.
#' - subspace: Computes the geodesic distance between each pair of networks
#' on the Grassmann manifold for a given input k, then scales each
#' distance between 0 and 1.
#' Default is all elements.
#' @param k is the number of eigenvectors to include in the subspace representation
#' of each network. Default is 5.
#' @param projections A list of all projections for the results. Can be NULL if not using
#' subspace similarity or if inhibitory negative edges are present.
#' @param projectionIdentity The projection of the identity matrix. Can be NULL if not
#' using subspace similarity.
#' @returns A scalar similarity value.
ComputeSimilarityForGroup <- function(results, comparisons, metric, k = 5, projections = NULL,projectionIdentity = NULL){
  
  # Compute the similarity for the group.
  groupSimilarity <- unlist(lapply(comparisons, function(comparison){
    similarity <- NULL

    # Obtain both networks. If one is missing from the results, stop.
    network1 <- results@results[[comparison@source]]
    network2 <- results@results[[comparison@target]]
    
    if(is.null(network1)){
      stop(paste0("Invalid network '", comparison@source, "' in comparisons object!"))
    }
    if(is.null(network2)){
      stop(paste0("Invalid network '", comparison@target, "' in comparisons object!"))
    }
    
    # Perform analysis.
    cat(".")
    if(metric == "jaccard"){
      similarity <- JaccardSim(network1 = network1, network2 = network2)
    }else if(metric == "degree"){
      similarity <- DegreeSim(network1 = network1, network2 = network2)
    }else if(metric == "modularity"){
      similarity <- ModularitySim(network1 = network1, network2 = network2)
    }else if(metric == "subspace"){
      similarity <- SubspaceSim(projection1 = projections[[comparison@source]], 
                                projection2 = projections[[comparison@target]], 
                                k = k)
    }else{
      stop(paste("Invalid metric to use for similarity:", metric))
    }
    return(similarity)
  }))
  
  # Compute the average similarity.
  return(mean(groupSimilarity))
}

#' Helper function to compute the Jaccard similarities between two networks.
#' This is a weighted version of the Jaccard similarity; i.e., 
#' @param network1 A network for comparison.
#' @param network2 Another network for comparison.
#' @returns The similarity value.
JaccardSim <- function(network1, network2){

  # First, find the intersection and union of edges.
  sharedEdges <- intersect(rownames(network1), rownames(network2))
  allEdges <- union(rownames(network1),rownames(network2))
  
  # Sum together all intersecting edges. Use the minimum value.
  intersectionSum <- sum(unlist(lapply(sharedEdges, function(e){
    tryCatch({return(min(network1[e, 3], network2[e, 3]))},
             error = function(cond){
               message(cond)
             })
  })))
  
  # For all edges in either one of the networks, zero out the value
  # in the other network if not present, and then find the maximum value.
  unionNet1 <- network1[allEdges, 3]
  unionNet1[which(is.na(unionNet1))] <- 0
  unionNet2 <- network2[allEdges, 3]
  unionNet2[which(is.na(unionNet2))] <- 0
  unionSum <- sum(pmax(unionNet1, unionNet2))
  jaccard <- intersectionSum / unionSum
  
  # If one network is empty, the overlap is 0.
  if(length(sharedEdges) == 0 && length(allEdges) == 0){
    jaccard <- NA
  }
  return(jaccard)
}

#' Helper function to compute the Degree similarities between two networks.
#' @param network1 A network for comparison.
#' @param network2 Another network for comparison.
#' @returns The similarity value.
DegreeSim <- function(network1, network2){
  
  degreeOverlap <- NA
  
  # If both data frames are empty, return NA. Else, calculate.
  if(nrow(network1) == 0 && nrow(network2) == 0){
    degreeOverlap <- NA
  }
  else if(nrow(network1) == 0 || nrow(network2) == 0){
    degreeOverlap <- 0
  }else if(nrow(network1) > 0 && nrow(network2) > 0){
    
    # Expand the networks
    expandedNetworks <- ExpandNetworks(network1, network2)
    expandedNetwork1 <- expandedNetworks$network1
    expandedNetwork2 <- expandedNetworks$network2
    
    # Create adjacency matrices.
    colnames(expandedNetwork1) <- c("source", "target", "weight")
    colnames(expandedNetwork2) <- c("source", "target", "weight")
    net1Adj <- as.matrix(igraph::as_adjacency_matrix(igraph::graph_from_data_frame(expandedNetwork1, 
                                                                             directed = TRUE), 
                                               attr = "weight"))
    net2Adj <- as.matrix(igraph::as_adjacency_matrix(igraph::graph_from_data_frame(expandedNetwork2, 
                                                                             directed = TRUE), 
                                               attr = "weight"))
    
    # Rearrange matrices.
    allNodes1 <- c(network1[,1], network1[,2])
    allNodes2 <- c(network2[,1], network2[,2])
    nodeOrder <- sort(union(allNodes1, allNodes2))
    net1Order <- net1Adj[nodeOrder, nodeOrder]
    net2Order <- net2Adj[nodeOrder, nodeOrder]
    
    # Sum each row and each column.
    net1SumOut <- rowSums(net1Order)
    net2SumOut <- rowSums(net2Order)
    net1SumIn <- colSums(net1Order)
    net2SumIn <- colSums(net2Order)
    
    # Sum the edge counts.
    net1PercentageOut <- net1SumOut / max(sum(net1SumOut), .Machine$double.xmin) 
    net2PercentageOut <- net2SumOut / max(sum(net2SumOut), .Machine$double.xmin)
    net1PercentageIn <- net1SumIn / max(sum(net1SumIn), .Machine$double.xmin)
    net2PercentageIn <- net2SumIn / max(sum(net2SumIn), .Machine$double.xmin)
    
    # Compute the overall differences.
    netDiffOut <- abs(net1PercentageOut - net2PercentageOut) / 2
    netDiffIn <- abs(net1PercentageIn - net2PercentageIn) / 2
    
    # Convert to similarities.
    degreeOverlapOut <- 1 - sum(netDiffOut)
    degreeOverlapIn <- 1 - sum(netDiffIn)
    
    # Compute F1 score.
    degreeOverlap <- 0
    if(degreeOverlapIn != 0 || degreeOverlapOut != 0){
      degreeOverlap <- (2 * degreeOverlapIn * degreeOverlapOut) / (degreeOverlapIn + degreeOverlapOut)
    }
  }
  return(degreeOverlap)
}

#' Helper function to expand the networks to include zero-weighted edges.
#' @param network1 A network for comparison.
#' @param network2 Another network for comparison.
#' @returns A list of two expanded networks.
ExpandNetworks <- function(network1, network2){
  
  # Initialize a list to return.
  expanded <- list(network1 = network1, network2 = network2)
  
  # Expand each network to include zeros.
  # To add these nodes to the adjacency list, they only need to show up once.
  allNodes1 <- c(network1[,1], network1[,2])
  allNodes2 <- c(network2[,1], network2[,2])
  extraNodes1 <- setdiff(allNodes1, allNodes2)
  extraNodes2 <- setdiff(allNodes2, allNodes1)
  if(length(extraNodes1) > 0 || length(extraNodes2) > 0){
    # If the data frames are not empty, create the zeros using connections from
    # existing nodes to extra nodes. If they are empty, set everything in the other
    # network to 0 and add it.
    extraZeros1 <- NULL
    extraZeros2 <- NULL
    if(nrow(network1) > 0){
      extraZeros1 <- data.frame(source = rep(unique(network1[,1]), length(extraNodes2)),
                                target = rep(extraNodes2, each = length(unique(network1[,1]))), 
                                             score = rep(0, length(extraNodes2)))
      rownames(extraZeros1) <- paste(extraZeros1$source, extraZeros1$target, sep = "_")
      colnames(extraZeros1) <- colnames(network1[c(1, 2, 3)])
      extraZeros1 <- rbind(network1, extraZeros1)
    }else{
      extraZeros1 <- network2
      extraZeros1[,3] <- rep(0, nrow(network2))
    }
    if(nrow(network2) > 0){
      extraZeros2 <- data.frame(source = rep(unique(network2[,1]), length(extraNodes1)),
                                target = rep(extraNodes1, each = length(unique(network2[,1]))), 
                                             score = rep(0, length(extraNodes1)))
      rownames(extraZeros2) <- paste(extraZeros2$source, extraZeros2$target, sep = "_")
      colnames(extraZeros2) <- colnames(network1[c(1, 2, 3)])
      extraZeros2 <- rbind(network2, extraZeros2)
    }else{
      extraZeros2 <- network1
      extraZeros2[,3] <- rep(0, nrow(network1))
    }
    expanded <- list(network1 = extraZeros1, network2 = extraZeros2)
  }
  return(expanded)
}

#' Helper function to compute the Modularity similarities between two networks.
#' @param network1 A network for comparison.
#' @param network2 Another network for comparison.
#' @returns The similarity value.
ModularitySim <- function(network1, network2){
  
  similarity <- NA
  
  # If one data frame is empty, similarity is 0.
  if(nrow(network1) > 0 && nrow(network2) > 0){
    
    # Expand the networks
    expandedNetworks <- ExpandNetworks(network1, network2)
    network1 <- expandedNetworks$network1
    network2 <- expandedNetworks$network2
    
    # Compute the community structure of the first network using ALPACA functions.
    # If one of the networks is empty, ALPACA will throw an error, so just return 0.
    communities <- alpacaWBMlouvain(network1)[[1]]
    
    # Remove duplicates from community structure, which can happen in the case of
    # non-bipartite networks.
    uniqueNodes <- unique(names(communities))
    for(node in uniqueNodes){
      whichNode <- which(names(communities) == node)
      if(length(whichNode) > 1){
        communities <- communities[-whichNode[-1]]
      }
    }

    # Normalize the edge weights in the original network.
    colnames(network1)[3] <- "weight"
    colnames(network2)[3] <- "weight"
    w <- igraph::as_adjacency_matrix(igraph::graph_from_data_frame(network1, 
                                                                   directed = TRUE), 
                                     attr = "weight")
    A <- igraph::as_adjacency_matrix(igraph::graph_from_data_frame(network2, 
                                                                   directed = TRUE), 
                                     attr = "weight")
    w_tilde <- w * sum(A) / length(which(network1$weight > 0))

    # Compute the expected community structure of the baseline network Nij.
    # The sum of weights from i to all nodes in each community.
    nodeCommSrcSums <- as.matrix(do.call(cbind, lapply(sort(unique(communities)), function(community){
      return(data.frame(v1=unlist(lapply(1:nrow(w_tilde), function(gene){
        return(sum(w_tilde[gene, names(communities)[which(communities == community)]]))
      }))))
    })))
    
    # The sum of weights from all nodes in each community to j.
    nodeCommTargetSums <- as.matrix(do.call(cbind, lapply(1:nrow(w_tilde), function(gene){
      return(data.frame(v1 = unlist(lapply(sort(unique(communities)), function(community){
        return(sum(w_tilde[names(communities)[which(communities == community)], gene]))
      }))))
    })))
    
    # The sum of weights between all communities.
    nodeCrossCommSums <- as.matrix(do.call(cbind, lapply(sort(unique(communities)), function(community){
      return(data.frame(v1 = unlist(lapply(sort(unique(communities)), function(community2){
        return(sum(w_tilde[names(communities)[which(communities == community)], 
                           names(communities)[which(communities == community2)]]))
      }))))
    })))
    # Build first part of numerator (the scaled weight from node i to the community of to node j).
    numeratorTerm1 <- nodeCommSrcSums[,communities]
    
    # Build second part of numerator (the scaled weight from the community of node i to node j).
    numeratorTerm2 <- nodeCommTargetSums[communities,]
    
    # Build the denominator.
    denominatorStep1 <- nodeCrossCommSums[communities,]
    denominatorStep2 <- t(denominatorStep1[,communities])
    
    # Compute the expected community structure.
    N <- (numeratorTerm1 * numeratorTerm2) / denominatorStep2
    N[which(is.nan(N))] <- 0
    
    # Compute differential modularity per edge. In our implementation, we do not
    # maximize differemtial modularity with respect to M.
    D <- sum(abs(A - N))
    
    # Normalize.
    DNorm <- D / (sum(A) + sum(N))
    
    # Convert to similarity.
    similarity <- 1 - DNorm
  }else if(xor(nrow(network1) == 0,  nrow(network2) == 0)){
    similarity = 0
  }else{
    similarity <- NA
  }
  return(similarity)
}

#' Helper function to compute the Subspace similarities between two networks.
#' @param projection1 The projection of network 1.
#' @param projection2 The projection of network 2.
#' @param k Number of eigenvectors to include in the subspace projection. 
#' @returns The similarity value.
#' @import irlba
SubspaceSim <- function(projection1, projection2, k){
  
  subspaceSim <- NA
  
  if(max(abs(projection1)) > 0 || max(abs(projection2)) > 0){
    # Compute the projection distance between networks 1 and 2 as per Ding et al.
    projectionDist <- k - sum(diag(projection1 %*% t(projection1) %*% projection2 %*% t(projection2)))
    
    # The normalization factor should be k. Note that if we were to include the
    # zero matrix as projection1 or projection2, this is what we would get.
    projectionDistNormFactor <- k
    
    # Round down to 0 if the differences are very small.
    cutoff <- 0.00000000001
    if(projectionDistNormFactor < cutoff && projectionDist < cutoff){
      projectionDist <- 0
      projectionDistNormFactor <- 1
    }
    
    # Set similarity.
    subspaceSim <- 1 - (projectionDist / projectionDistNormFactor)
  }
  
  # Return normalized similarity.
  return(subspaceSim)
}

#' Helper function to compute the projection of a network.
#' @param network A network to project (N x N).
#' @param k Number of eigenvectors to include in the subspace projection. 
#' @param allNodes Additional nodes to add to the adjacency matrix.
#' @returns The projection (a matrix of N x k)
#' @import irlba
ComputeProjection <- function(network, k, allNodes){
  set.seed(1)
  
  # Check that k is a valid value.
  if(!is.numeric(k) || k - floor(k) != 0 || k <= 0){
    stop("Number of eigenvectors k must be an integer greater than 0!")
  }
  if(k >= length(unique(allNodes))){
    stop("Number of eigenvectors must be less than the dimensionality of the adjacency matrix!")
  }
  
  # Add all relationships with nodes that are not in the network.
  # To add these nodes to the adjacency list, they only need to show up once.
  colnames(network)[3] <- "weight"
  extraNodes <- setdiff(allNodes, c(network[,1], network[,2]))
  expandedNetwork <- network
  if(nrow(network) == 0){
    expandedNetwork <- data.frame(source = extraNodes,
                             target = extraNodes, weight = rep(0, length(extraNodes)))
  }else{
    extraZeros <- data.frame(source = extraNodes,
                             target = extraNodes, weight = rep(0, length(extraNodes)))
    rownames(extraZeros) <- paste(extraZeros$source, extraZeros$target, sep = "_")
    colnames(extraZeros) <- colnames(network[c(1,2,3)])
    expandedNetwork <- rbind(network, extraZeros)
  }

  # Compute the adjacency matrix.
  netAdj <- igraph::as_adjacency_matrix(igraph::graph_from_data_frame(expandedNetwork, 
                                          directed = TRUE), attr = "weight")

  # Compute the normalized graph Laplacian.
  degree <- matrix(data = rep(0, ncol(netAdj) * nrow(netAdj)), nrow = nrow(netAdj))
  laplacian <- degree - netAdj

  # Obtain the first k eigenvectors, of L, which corresponds to the
  # directions of highest variance of L. This is projection captures the
  # most information from L.
  projection <- NULL
  tryCatch({
    if(sum(laplacian) == 0){
      projection <- matrix(data = rep(0, nrow(laplacian) * k), nrow = nrow(laplacian))
      colnames(projection) <- paste0(rep("PC", k), as.character(seq(1, k, 1)))
    }else{
      projection <- as.matrix(irlba::prcomp_irlba(laplacian, n = k)$rotation)
    }
  }, error = function(cond){
    print(cond)
    stop("Could not project network onto a subspace!")
  })

  # Return the projection.
  return(projection)
}