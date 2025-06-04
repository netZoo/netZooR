#' Given a set of genes of interest, full bipartite networks with scores (one network for each sample), a significance
#' cutoff for statistical testing, and a hop constraint, BLOBFISH finds a subnetwork of
#' significant edges connecting the genes.
#' @param geneSet A character vector of genes comprising the targets of interest.
#' @param networks A list of bipartite (PANDA-like) networks, where each network is a data frame with the following format:
#' tf,gene,score
#' @param alpha The significance cutoff for the statistical test.
#' @param hopConstraint The maximum number of hops to be considered between gene pairs.
#' Must be an even number.
#' @param nullDistribution The null distribution, specified as a vector of values.
#' @param verbose Whether or not to print detailed information about the run.
#' @param topX Select the X lowest significant p-values for each gene. NULL by default.
#' @param doFDRAdjustment Whether or not to perform FDR adjustment.
#' @param pValueChunks The number of chunks to split when calculating the p-value. This
#' parameter allows the edges to be split into chunks to prevent memory errors.
#' @param pValueFile The file where the p-values should be saved. If NULL, they are not
#' saved and need to be recalculated.
#' @param loadPValues Whether p-values should be loaded from pValueFile or re-generated.
#' Default is FALSE.
#' @returns A bipartite subnetwork in the same format as the original networks.
#' @export
RunBLOBFISH <- function(geneSet, networks, alpha, hopConstraint, nullDistribution,
                        verbose = FALSE, topX = NULL, doFDRAdjustment = TRUE,
                        pValueChunks = 100, loadPValues = FALSE, pValueFile = "pvalues.RDS"){
  
  # Check for invalid inputs.
  if(!is.character(geneSet) || (!is.list(networks) && !is.data.frame(networks)) || !is.numeric(alpha) || !is.numeric(hopConstraint)){
    stop(paste("Wrong input type! geneSet must be a character vector. networks must be a list.",
               "alpha and hopConstraint must be scalar numeric values."))
  }else if(!is.data.frame(networks) && !is.data.frame(networks[[1]])){
    stop("Each network must be a data frame.")
  }else if(!is.data.frame(networks) && (ncol(networks[[1]]) != 3 || colnames(networks[[1]])[1] != "tf" ||
           colnames(networks[[1]])[2] != "gene" || colnames(networks[[1]])[3] != "score")){
    stop(paste("Each network must have transcription factors in the first column,",
               "target genes in the second column, and scores in the third column."))
  }else if(is.data.frame(networks) && (ncol(networks) < 3 || colnames(networks)[1] != "tf" ||
                                colnames(networks)[2] != "gene")){
    stop(paste("Each network must have transcription factors in the first column,",
               "target genes in the second column, and scores in the third column."))
  }else if(alpha > 1 || alpha <=0){
    stop("alpha must be between 0 and 1, not including 0.")
  }else if(hopConstraint < 2 || hopConstraint %% 2 != 0){
    stop("hopConstraint must be an even number of at least 2.")
  }else if(!is.numeric(nullDistribution)){
    stop("nullDistribution must be numeric.")
  }
  
  # Build the subnetwork.
  subnetwork <- BuildSubnetwork(geneSet = geneSet, 
                                networks = networks, 
                                alpha = alpha, 
                                hopConstraint = hopConstraint,
                                nullDistribution = nullDistribution,
                                verbose = verbose, topX = topX,
                                doFDRAdjustment = doFDRAdjustment,
                                pValueChunks = pValueChunks, 
                                loadPValues = loadPValues, 
                                pValueFile = pValueFile)
  
  # Return.
  return(subnetwork)
}

#' Generate a null distribution of edge scores for PANDA-like networks; that is,
#' the set of edges where (1) the TF does not have a binding motif in the gene region,
#' (2) the TF does not form a complex with any other TF that has a binding motif in
#' the gene region, and (3) the genes regulated by the TF are not coexpressed with the
#' gene in question. We obtain this by inputting an empty prior and an identity coexpression
#' matrix.
#' @param ppiFile The location of the protein-protein interaction network between transcription factors.
#' This should be a TSV file where the first two columns are the transcription
#' factors and the third is whether there is a PPI between them.
#' @param motifFile The location of the motif prior from genes to transcription factors. This should
#' be a TSV file where the first column is the transcription factors, the
#' second is the genes, and the third is whether the transcription factor's
#' binding motif is in the gene promoter region.
#' @param sampSize Number of samples to simulate
#' @param numberOfPandas Number of null PANDA networks to generate
#' @export
GenerateNullPANDADistribution <- function(ppiFile, motifFile, sampSize = 20,
                                          numberOfPandas = 10){
  
  # Set the seed.
  set.seed(1)
  
  # Read the motif.
  motif <- utils::read.table(motifFile, sep = "\t")
  ppi <- utils::read.table(ppiFile, sep = "\t")
  
  # Generate the null PANDA edges.
  nullPandas <- lapply(1:numberOfPandas, function(i){
    
    # Generate a random matrix of expression values, where any correlation that exists
    # is spurious.
    ngenes <- length(unique(motif[,2]))
    randomExpression <- matrix(data = stats::rnorm(ngenes * sampSize), nrow = ngenes)
    rownames(randomExpression) <- unique(motif[,2])

    # Save in a temporary file.
    fname <- paste0("tmp_", i, ".tsv")
    write.table(randomExpression, fname, row.names = TRUE, col.names = FALSE, sep = "\t",
                quote = FALSE)
    
    # Run PANDA.
    nullPanda <- pandaPy(expr_file = fname, motif_file=motifFile, ppi_file=ppiFile, save_tmp=FALSE,
                         save_memory=TRUE)$panda
    rownames(nullPanda) <- paste(nullPanda$TF, nullPanda$Gene, sep = "__")
    hist(nullPanda[,3])

    # Remove file.
    unlink(fname)
    
    # Return null results.
    rownames(motif) <- paste(motif[,1], motif[,2], sep = "__")
    rowsToExclude <- unlist(lapply(rownames(motif), function(edge){
      
      # Find other transcription factors that interact with this one.
      tf <- strsplit(edge, "__")[[1]][1]
      interactingTF <- unique(c(ppi[which(ppi[,2] == tf),1], ppi[which(ppi[,1] == tf),2]))
      
      # Create new edges including other transcription factors.
      newEdges <- paste(interactingTF, edge, sep = "__")

      # Return the old and new edges.
      return(c(edge, newEdges))
    }))
    
    # Return the scores for all tf-gene relationships not in the original motif.
    return(nullPanda[setdiff(rownames(nullPanda), rowsToExclude),"Score"])
  })

  # Return the values.
  nullPandasAll <- unlist(nullPandas)
  return(sample(nullPandasAll, n = length(nullPandasAll)))
}

#' Find the subnetwork of significant edges connecting the genes.
#' @param geneSet A character vector of genes comprising the targets of interest.
#' @param networks A list of bipartite (PANDA-like) networks, where each network is a data frame with the following format:
#' tf,gene,score
#' @param alpha The significance cutoff for the statistical test.
#' @param hopConstraint The maximum number of hops to be considered between gene pairs.
#' Must be an even number.
#' @param nullDistribution The null distribution, specified as a vector of values.
#' @param verbose Whether or not to print detailed information about the run.
#' @param topX Select the X lowest significant p-values for each gene. NULL by default.
#' @param doFDRAdjustment Whether or not to perform FDR adjustment.
#' @param pValueChunks The number of chunks to split when calculating the p-value. This
#' parameter allows the edges to be split into chunks to prevent memory errors.
#' @param pValueFile The file where the p-values should be saved. If NULL, they are not
#' saved and need to be recalculated.
#' @param loadPValues Whether p-values should be loaded from pValueFile or re-generated.
#' Default is FALSE.
#' @returns A bipartite subnetwork in the same format as the original networks.
BuildSubnetwork <- function(geneSet, networks, alpha, hopConstraint, nullDistribution,
                            verbose = FALSE, topX = NULL, doFDRAdjustment = TRUE,
                            pValueChunks = 100, loadPValues = FALSE, pValueFile = "pvalues.RDS"){
  
  # Name edges for each network.
  combinedNetwork <- networks
  if(!is.data.frame(combinedNetwork)){
    networksNamed <- lapply(networks, function(network){
      rownames(network) <- paste(network$tf, network$gene, sep = "__")
      return(network)
    })
    # Paste together the networks.
    combinedNetwork <- networksNamed[[1]]
    for(i in 2:length(networksNamed)){
      combinedNetwork[,2+i] <- networksNamed[[i]]$score
    }
  }
  
  # Find all significant edges in the network.
  pValues <- rep(NA, nrow(combinedNetwork))
  
  # If we are calculating p-values for the first time, calculate and save them.
  if(loadPValues == FALSE){
    pValues <- CalculatePValues(network = combinedNetwork,
                                pValueChunks = pValueChunks,
                                nullDistribution = nullDistribution,
                                doFDRAdjustment = doFDRAdjustment,
                                pValueFile = pValueFile)
  }else{
    # Read the saved p-values.
    pValues <- readRDS(pValueFile)
  }
  
  # Subset the network.
  whichSig <- which(pValues < alpha)
  significantEdges <- rownames(combinedNetwork)[whichSig]
  subnetwork <- combinedNetwork[significantEdges, c(1:2)]
  pValues <- pValues[significantEdges]
  genesWithNoSigEdges <- setdiff(geneSet, subnetwork$gene)
  geneSet <- intersect(geneSet, subnetwork$gene)
  if(verbose == TRUE){
    message(paste("The following genes had no significant edges:", paste(genesWithNoSigEdges, collapse = ",")))
    message(paste("Retained", length(significantEdges), "out of", length(rownames(combinedNetwork)), "edges"))
  }
  
  # For each gene, find the significant edges from each hop.
  significantSubnetworks <- FindSignificantEdgesForHop(geneSet = geneSet, pValues = pValues,
                                                       combinedNetwork = subnetwork,
                                                       hopConstraint = hopConstraint / 2,
                                                       verbose = verbose, topX = topX)
  
  # Find matches for each hop. Note that we do not need to consider cases where the 
  # number of hops are uneven. For instance, any two nodes that can be connected by
  # a node 3 hops away from node 1 and 1 hop away from node 2 are also connected by
  # a node 2 hops away from both nodes 1 and 2. The same goes for even numbers.
  # Even-odd pairs should not be considered. This can be proven.
  subnetwork <- FindConnectionsForAllHopCounts(subnetworks = significantSubnetworks,
                                               verbose = verbose)
  return(subnetwork)
}

#' Calculate p-values for all edges in the network using a Wilcoxon two-sample test
#' for each edge.
#' @param network A combination of PANDA-like networks, with the following format 
#' (e.g., 3 networks), provided as a data frame:
#' tf,gene,score1,score2,score3
#' @param nullDistribution The null distribution, specified as a vector of values.
#' @param pValueChunks The number of chunks to split when calculating the p-value. This
#' parameter allows the edges to be split into chunks to prevent memory errors.
#' @param doFDRAdjustment Whether or not to perform FDR adjustment.
#' @param pValueFile The file where the p-values should be saved. If NULL, they are not
#' saved and need to be recalculated.
#' @param verbose Whether or not to print detailed information about the run.
#' @returns A vector of p-values, one for each edge.
#' @export
CalculatePValues <- function(network, nullDistribution, pValueChunks = 100, 
                             doFDRAdjustment = TRUE, pValueFile = "pvalues.RDS",
                             verbose = FALSE){
  
  # Initialize p-values.
  pValues <- rep(NA, nrow(network))
  
  # Set the initial start and end indices.
  startIndex <- 1
  endIndex <- min(startIndex + ceiling(nrow(network) / pValueChunks),
                  nrow(network))
  i <- 1
  while(i < pValueChunks && startIndex <= endIndex){
    
    # Calculate p-values for this chunk.
    ourEdgeVals <- network[startIndex:endIndex, 3:ncol(network)]
    nullEdgeVals <- t(matrix(rep(nullDistribution, 
                                 nrow(ourEdgeVals)), ncol = nrow(ourEdgeVals)))
    pValues[startIndex:endIndex] <- matrixTests::row_wilcoxon_twosample(x = ourEdgeVals, 
                                                                        y = nullEdgeVals, 
                                                                        alternative = "greater")$pvalue
    
    # Print status.
    if(verbose == TRUE){
      message(paste("Completed p-values for chunk", i, "out of", pValueChunks))
    }
    
    # Update indices.
    startIndex <- endIndex + 1
    endIndex <- min(startIndex + ceiling(nrow(network) / pValueChunks),
                    nrow(network))
    i <- i + 1
  }
  
  # Adjust the p-values.
  if(doFDRAdjustment == TRUE){
    pValues <- stats::p.adjust(pValues, method = "fdr")
  }
  names(pValues) <- rownames(network)
  
  # Save the p-values.
  if(!is.null(pValueFile)){
    saveRDS(pValues, pValueFile)
  }
  
  # Return the p-values.
  return(pValues)
}
  
#' Find the subnetwork of significant edges n / 2 hops away from each gene.
#' @param geneSet A character vector of genes comprising the targets of interest.
#' @param combinedNetwork A concatenation of n PANDA-like networks with the following format:
#' tf,gene,score_net1, score_net2, ... , score_netn
#' @param pValues The p-values for all edges.
#' @param hopConstraint The maximum number of hops to be considered for a gene.
#' @param verbose Whether or not to print detailed information about the run.
#' @param topX Select the X lowest significant p-values for each gene. NULL by default.
#' @returns A bipartite subnetwork in the same format as the original networks.
FindSignificantEdgesForHop <- function(geneSet, combinedNetwork, hopConstraint, pValues,
                                       verbose = FALSE, topX = NULL){
  # Build the significant subnetwork for each gene, up to the hop constraint.
  uniqueGeneSet <- sort(unique(geneSet))
  geneSubnetworks <- lapply(uniqueGeneSet, function(gene){
    
    # Get all significant edges for a 1-hop subnetwork.
    if(verbose == TRUE){
      message(paste("Evaluating hop 1 for gene", gene))
    }
    subnetwork1Hop <- SignificantBreadthFirstSearch(networks = combinedNetwork, 
                                                    pValues = pValues, 
                                                    startingNodes = gene,
                                                    nodesToExclude = c(),
                                                    startFromTF = FALSE, 
                                                    verbose = verbose,
                                                    topX = topX)
    
    # Set the starting and excluded set for the next hop.
    startingNodes <- unique(subnetwork1Hop$tf)
    topXNew <- NULL
    if(!is.null(topX)){
      topXNew <- topX * length(startingNodes)
    }
    excludedSubset <- gene
    
    # Add to the list of all subnetworks.
    allSubnetworksForGene <- list(subnetwork1Hop)
    
    # Loop until we reach the maximum number of hops or there are no new edges
    # to traverse.
    hop <- 2
    while(hop <= hopConstraint && length(startingNodes) > 0){
      
      # If we are on an even hop, start from transcription factors.
      # If we are on an odd hop, start from genes.
      if(hop %% 2 == 0){
        
        # Find all significant edges in the next hop.
        if(verbose == TRUE){
          message(paste("Evaluating hop", hop, "for gene", gene))
        }
        subnetworkHops <- SignificantBreadthFirstSearch(networks = combinedNetwork, 
                                                        pValues = pValues, 
                                                        startingNodes = startingNodes,
                                                        nodesToExclude = excludedSubset,
                                                        startFromTF = TRUE, 
                                                        verbose = verbose,
                                                        topX = topXNew)
        
        # Set the starting and excluded set for the next hop.
        excludedSubset <- c(excludedSubset, startingNodes)
        startingNodes <- setdiff(unique(subnetworkHops$gene), excludedSubset)
        if(!is.null(topX)){
          topXNew <- topX * length(startingNodes)
        }
      }else{
        
        # Find all significant edges in the next hop.
        if(verbose == TRUE){
          message(paste("Evaluating hop", hop, "for gene", gene))
        }
        subnetworkHops <- SignificantBreadthFirstSearch(networks = combinedNetwork, 
                                                        pValues = pValues, 
                                                        startingNodes = startingNodes,
                                                        nodesToExclude = excludedSubset,
                                                        startFromTF = FALSE, 
                                                        verbose = verbose, 
                                                        topX = topXNew)
        
        # Set the starting and excluded set for the next hop.
        excludedSubset <- c(excludedSubset, startingNodes)
        startingNodes <- setdiff(unique(subnetworkHops$tf), excludedSubset)
        if(!is.null(topX)){
          topXNew <- topX * length(startingNodes)
        }
      }
      
      # Add to the list.
      allSubnetworksForGene[[length(allSubnetworksForGene) + 1]] <- subnetworkHops
      
      # Increment hops.
      hop <- hop + 1
    }
    return(allSubnetworksForGene)
  })
  
  # Add the names of the genes.
  names(geneSubnetworks) <- uniqueGeneSet
  return(geneSubnetworks)
}

#' Find all significant edges adjacent to the starting nodes, excluding the nodes
#' specified.
#' @param networks A concatenation of n PANDA-like networks with the following format:
#' tf,gene,score_net1, score_net2, ... , score_netn
#' Edges must be specified as "tf__gene".
#' @param pValues The p-values from the original network.
#' @param startingNodes The list of nodes from which to start.
#' @param nodesToExclude The list of nodes to exclude from the search.
#' @param startFromTF Whether to start from transcription factors (TRUE) or genes (FALSE).
#' @param verbose Whether or not to print detailed information about the run.
#' @param topX Select the X lowest significant p-values for each gene. NULL by default.
#' @returns A bipartite subnetwork in the same format as the original networks.
SignificantBreadthFirstSearch <- function(networks, pValues, startingNodes,
                                          nodesToExclude, startFromTF, 
                                          verbose = FALSE, topX = NULL){
  # Check that provided nodes overlap with the networks.
  if((length(setdiff(startingNodes, networks$tf)) > 0 && startFromTF == TRUE) ||
     (length(setdiff(startingNodes, networks$gene)) > 0 && startFromTF == FALSE)){
    stop("ERROR: Starting nodes do not overlap with network nodes")
  }
  if(length(setdiff(nodesToExclude, c(networks$tf, networks$gene))) > 0){
    stop("ERROR: List of nodes to exclude does not overlap with network nodes")
  }
  if(length(intersect(startingNodes, nodesToExclude)) > 0){
    stop("ERROR: Starting nodes cannot overlap with nodes to exclude")
  }
  
  # Identify genes and transcription factors to test, based on which of these we are
  # starting from.
  tfsToTest <- c()
  genesToTest <- c()
  if(startFromTF == TRUE){
    tfsToTest <- startingNodes
    genesToTest <- setdiff(unique(networks$gene), nodesToExclude)
  }else{
    genesToTest <- startingNodes
    tfsToTest <- setdiff(unique(networks$tf), nodesToExclude)
  }
  
  # Construct all edges to test based on the combination of these.
  geneLongList <- rep(genesToTest, length(tfsToTest))
  tfLongList <- unlist(lapply(tfsToTest, function(tf){
    return(rep(tf, length(genesToTest)))
  }))
  allPossibleEdges <- paste(tfLongList, geneLongList, sep = "__")
  allEdges <- intersect(allPossibleEdges, rownames(networks))
  
  # For each edge, measure its significance.
  subnetwork <- networks
  if(length(allEdges) > 0){
    
    # If topX is specified, filter again.
    significantEdges <- allEdges
    if(!is.null(topX) && length(allEdges) > topX){
      whichTopX <- order(pValues[allEdges])[1:topX]
      significantEdges <- allEdges[whichTopX]
    }
    
    # Return the edges meeting alpha.
    subnetwork <- networks[significantEdges, c(1:2)]
    if(verbose == TRUE){
      message(paste("Retained", length(significantEdges), "edges"))
    }
  }
  
  # Return the subnetwork.
  return(subnetwork)
}

#' For all hop counts up to the maximum, find subnetworks connecting each pair of
#' genes by exactly that number of hops. For instance, find each 
#' @param subnetworks A list of bipartite (PANDA-like) subnetworks for each gene, 
#' containing only the significant edges meeting the hop count criteria and
#' where each network is a data frame with the following format:
#' tf,gene
#' @param verbose Whether or not to print detailed information about the run.
#' @returns A bipartite subnetwork in the same format as the original networks.
FindConnectionsForAllHopCounts <- function(subnetworks, verbose = FALSE){
  
  # Find a subnetwork for each hop count.
  hopCountSubnetworks <- lapply(1:length(subnetworks[[1]]), function(hops){
    
    # For each pair of genes, find the subnetworks for this number of hops.
    geneSpecificHopCountSubnetwork <- lapply(1:(length(names(subnetworks))-1), function(i){
      genePairSpecificHopCountSubnetwork <- lapply((i+1):length(names(subnetworks)), function(j){
        
        # Get the subnetworks for the number of hops of interest.
        gene1 <- names(subnetworks)[i]
        gene2 <- names(subnetworks)[j]
        connectingSubnetwork <- data.frame(tf = NA, gene = NA)[0,]
        
        # If there were no edges at this hop count for one or both genes,
        # do not evaluate.
        if(length(subnetworks[[gene1]]) >= hops && length(subnetworks[[gene2]]) >= hops){
          subnetwork1 <- subnetworks[[gene1]][[hops]]
          subnetwork2 <- subnetworks[[gene2]][[hops]]
          
          # Initialize overlapping subnetwork.
          genesToRecurse1 <- c()
          genesToRecurse2 <- c()
          tfsToRecurse1 <- c()
          tfsToRecurse2 <- c()
          
          # If the number of hops is even, add edges from genes that overlap
          # If the number of hops is odd, add edges from transcription factors that overlap.
          if(hops %% 2 == 0){
            overlappingGenes <- intersect(subnetwork1$gene, subnetwork2$gene)
            if(verbose == TRUE){
              message(paste("Hop", hops, "-", length(overlappingGenes), "overlapped between", gene1, "and", gene2))
            }
            whichSubnet1Gene <- which(subnetwork1$gene %in% overlappingGenes)
            whichSubnet2Gene <- which(subnetwork2$gene %in% overlappingGenes)
            tfsToRecurse1 <- unique(subnetwork1[whichSubnet1Gene, "tf"])
            tfsToRecurse2 <- unique(subnetwork2[whichSubnet2Gene, "tf"])
            connectingSubnetwork <- rbind(connectingSubnetwork, subnetwork1[whichSubnet1Gene,],
                                          subnetwork2[whichSubnet2Gene,])
          }else{
            overlappingTF <- intersect(subnetwork1$tf, subnetwork2$tf)
            if(verbose == TRUE){
              message(paste("Hop", hops, "-", length(overlappingTF), "overlapped between", gene1, "and", gene2))
            }
            whichSubnet1TF <- which(subnetwork1$tf %in% overlappingTF)
            whichSubnet2TF <- which(subnetwork2$tf %in% overlappingTF)
            genesToRecurse1 <- unique(subnetwork1[whichSubnet1TF, "gene"])
            genesToRecurse2 <- unique(subnetwork2[whichSubnet2TF, "gene"])
            connectingSubnetwork <- rbind(connectingSubnetwork, subnetwork1[whichSubnet1TF,],
                                          subnetwork2[whichSubnet2TF,])
          }
          
          # Recurse back over the number of hops.
          if(hops-1 >= 1){
            for(hop in (hops-1):1){
              subnetwork1 <- subnetworks[[gene1]][[hop]]
              subnetwork2 <- subnetworks[[gene2]][[hop]]
              
              # If the current number of hops is even, add edges from TFs connected to genes of interest.
              # If the current number of hops is odd, add edges from genes connected to TFs of interest.
              if(hop %% 2 == 0){
                whichTFConnectedToGene1 <- which(subnetwork1$gene %in% genesToRecurse1)
                whichTFConnectedToGene2 <- which(subnetwork2$gene %in% genesToRecurse2)
                tfsToRecurse1 <- unique(subnetwork1[whichTFConnectedToGene1, "tf"])
                tfsToRecurse2 <- unique(subnetwork2[whichTFConnectedToGene2, "tf"])
                connectingSubnetwork <- rbind(connectingSubnetwork, subnetwork1[whichTFConnectedToGene1,],
                                              subnetwork2[whichTFConnectedToGene2,])
              }else{
                whichGeneConnectedToTF1 <- which(subnetwork1$tf %in% tfsToRecurse1)
                whichGeneConnectedToTF2 <- which(subnetwork2$tf %in% tfsToRecurse2)
                genesToRecurse1 <- unique(subnetwork1[whichGeneConnectedToTF1, "gene"])
                genesToRecurse2 <- unique(subnetwork2[whichGeneConnectedToTF2, "gene"])
                connectingSubnetwork <- rbind(connectingSubnetwork, subnetwork1[whichGeneConnectedToTF1,],
                                              subnetwork2[whichGeneConnectedToTF2,])
              }
            }
          }
        }
        
        # Return the subnetwork, which should now contain all of the edges connecting the
        # gene pair at the prespecified number of hops.
        return(connectingSubnetwork)
      })
      # Bind together the subnetwork for each gene pair.
      connectingSubnetworkAll <- do.call(rbind, genePairSpecificHopCountSubnetwork)
      return(connectingSubnetworkAll)
    })
    
    # Bind together the subnetworks for each gene.
    return(do.call(rbind, geneSpecificHopCountSubnetwork))
  })
  
  # Bind together the subnetworks for each hop count.
  compositeSubnetwork <- do.call(rbind, hopCountSubnetworks)
  compositeSubnetworkEdges <- paste(compositeSubnetwork$tf, compositeSubnetwork$gene, sep = "__")
  uniqueEdges <- sort(unique(compositeSubnetworkEdges))
  compositeSubnetworkDedup <- do.call(rbind, lapply(uniqueEdges, function(edge){
    whichFirstEdge <- which(compositeSubnetworkEdges == edge)[1]
    return(compositeSubnetwork[whichFirstEdge,])
  }))
  rownames(compositeSubnetworkDedup) <- uniqueEdges
  
  # Remove all genes connected to a single transcription factor. These genes were
  # added because they are regulated by a transcription factor that co-regulates
  # two seed genes. Similarly, remove all transcription factors connected to a 
  # single gene.
  geneCounts <- table(compositeSubnetworkDedup$gene)
  tfCounts <- table(compositeSubnetworkDedup$tf)
  genesToRemove <- names(geneCounts)[which(geneCounts == 1)]
  genesToRemove <- setdiff(genesToRemove, names(subnetworks))
  tfsToRemove <- names(tfCounts)[which(tfCounts == 1)]
  compositeSubnetworkDedup <- compositeSubnetworkDedup[which(compositeSubnetworkDedup$gene %in% setdiff(names(geneCounts), 
                                                                                                        genesToRemove)),]
  compositeSubnetworkDedup <- compositeSubnetworkDedup[which(compositeSubnetworkDedup$tf %in% setdiff(names(tfCounts), 
                                                                                                      tfsToRemove)),]
  return(compositeSubnetworkDedup)
}

#' Plot the networks, using different colors for transcription factors, genes of interest,
#' and additional genes.
#' @param network A data frame with the following format:
#' tf,gene
#' @param genesOfInterest Which genes of interest to highlight
#' @param tfColor Color for the transcription factors
#' @param geneColorMapping Color mapping from a set of genes to a color. The
#' nodes and edges connected to them will be this color. If NULL, all genes and
#' their edges will be gray. The format is a data frame, where the first column ("gene")
#' is the name of the gene and the second ("color") is the color.
#' @param nodeSize Size of node
#' @param edgeWidth Width of edges
#' @param vertexLabels Which vertex labels to include. By default, none are included.
#' @param layoutBipartite Whether or not to layout as a bipartite graph.
#' @param vertexLabelSize The size of label to use for the vertex, as a fraction of the default.
#' @param vertexLabelOffset Number of pixels in the offset when plotting labels.
#' Default is TRUE.
#' @returns A bipartite plot of the network
#' @export
PlotNetwork <- function(network, genesOfInterest,
                        tfColor = "blue", nodeSize = 1,
                        edgeWidth = 0.5, vertexLabels = NA, vertexLabelSize = 0.7,
                        vertexLabelOffset = 0.5, layoutBipartite = TRUE, geneColorMapping = NULL){
  
  # Convert from factor to character.
  network$tf <- as.character(network$tf)
  network$gene <- as.character(network$gene)
  
  # Set the node attributes.
  uniqueNodes <- unique(c(network$tf, network$gene))
  nodeAttrs <- data.frame(node = uniqueNodes,
                          color = rep("gray", length(uniqueNodes)),
                          size = rep(nodeSize, length(uniqueNodes)),
                          frame.width = rep(0, length(uniqueNodes)),
                          label.color = "black", label.cex = vertexLabelSize,
                          label.dist = vertexLabelOffset)
  rownames(nodeAttrs) <- uniqueNodes
  
  # Add TF colors.
  nodeAttrs[which(uniqueNodes %in% network$tf), "color"] <- tfColor
  
  # Add gene colors.
  rownames(geneColorMapping) <- geneColorMapping$gene
  geneColorMapping <- geneColorMapping[intersect(rownames(geneColorMapping), uniqueNodes),]
  if(!is.null(geneColorMapping)){
    nodeAttrs[rownames(geneColorMapping), "color"] <- geneColorMapping$color
  }

  # Add edge attributes.
  if(!is.null(geneColorMapping)){
    for(gene in rownames(geneColorMapping)){
      network[which(network$gene == gene), "color"] <- geneColorMapping[gene, "color"]
    }
  }
  network$width <- edgeWidth

  # Create a graph object.
  graph <- igraph::graph_from_data_frame(network, vertices = nodeAttrs, directed = FALSE)
  V(graph)$type <- V(graph)$name %in% network$tf
  
  # Plot.
  labels <- V(graph)$name
  whichEmpty <- which(labels %in% setdiff(labels, vertexLabels))
  labels[whichEmpty] <- rep(NA, length(whichEmpty))

  if(layoutBipartite == TRUE){
    LO <- layout_as_bipartite(graph)
    LO <- LO[,c(2,1)]
    igraph::plot.igraph(graph, layout = LO, vertex.label = labels)
  }else{
    igraph::plot.igraph(graph, vertex.label = labels)
  }
 }
