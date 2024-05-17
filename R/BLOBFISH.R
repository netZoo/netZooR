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
#' @returns A bipartite subnetwork in the same format as the original networks.
#' @export
RunBLOBFISH <- function(geneSet, networks, alpha, hopConstraint, nullDistribution,
                        verbose = FALSE, topX = NULL){
  
  # Check for invalid inputs.
  if(!is.character(geneSet) || !is.list(networks) || !is.numeric(alpha) || !is.numeric(hopConstraint)){
    stop(paste("Wrong input type! geneSet must be a character vector. networks must be a list.",
               "alpha and hopConstraint must be scalar numeric values."))
  }else if(!is.data.frame(networks[[1]])){
    stop("Each network must be a data frame.")
  }else if(ncol(networks[[1]]) != 3 || colnames(networks[[1]])[1] != "tf" ||
           colnames(networks[[1]])[2] != "gene" || colnames(networks[[1]])[3] != "score"){
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
                                verbose = verbose, topX = topX)
  
  # Return.
  return(subnetwork)
}

#' Generate a null distribution of edge scores for PANDA-like networks; that is,
#' the set of edges where (1) the TF does not have a binding motif in the gene region,
#' (2) the TF does not form a complex with any other TF that has a binding motif in
#' the gene region, and (3) the genes regulated by the TF are not coexpressed with the
#' gene in question. We obtain this by inputting an empty prior and an identity coexpression
#' matrix.
#' @param ngenes Number of genes to simulate
#' @param nTranscriptionFactors Number of transcription factors to simulate
#' @param sampSize Number of samples to simulate
#' @param numberOfPandas Number of null PANDA networks to generate
#' @export
GenerateNullPANDADistribution <- function(ngenes = 20000, nTranscriptionFactors = 600, sampSize = 20,
                                          numberOfPandas = 10){
  
  # Set the seed.
  set.seed(1)
  
  # Generate a motif file where only TFs 1-3 are connected to genes.
  motiffname <- paste0("tmpMotif.tsv")
  motifs <- data.frame(tf = rep(paste0("tf", 1:nTranscriptionFactors), ngenes),
                       gene = unlist(lapply(1:ngenes, function(j){
                         return(rep(paste0("gene", j), nTranscriptionFactors))
                       })), score = rep(0, ngenes * nTranscriptionFactors))
  motifs[intersect(which(motifs$tf == "tf1"), which(motifs$gene == "gene1")), "score"] <- 1
  motifs[intersect(which(motifs$tf == "tf2"), which(motifs$gene == "gene2")), "score"] <- 1
  motifs[intersect(which(motifs$tf == "tf3"), which(motifs$gene == "gene3")), "score"] <- 1
  write.table(motifs, motiffname, quote = FALSE, sep = "\t",
              row.names = FALSE, col.names = FALSE)
  
  # Generate a PPI where TFs 1 is connected to TFs 4-5.
  ppifname <- paste0("tmpPPI.tsv")
  ppi <- data.frame(tf = c(paste0("tf", 1:nTranscriptionFactors), "tf1", "tf1"),
                    gene = c(paste0("tf", 1:nTranscriptionFactors), "tf4", "tf5"), 
                    score = rep(1, nTranscriptionFactors + 2))
  write.table(ppi, ppifname, quote = FALSE, sep = "\t",
              row.names = FALSE, col.names = FALSE)
  
  # Generate the null PANDA edges.
  nullPandas <- lapply(1:numberOfPandas, function(i){
    
    # Generate a random matrix of expression values, where any correlation that exists
    # is spurious.
    randomExpression <- matrix(data = stats::rnorm(ngenes * sampSize), nrow = ngenes)
    rownames(randomExpression) <- paste0("gene", 1:ngenes)

    # Save in a temporary file.
    fname <- paste0("tmp_", i, ".tsv")
    write.table(randomExpression, fname, row.names = TRUE, col.names = FALSE, sep = "\t",
                quote = FALSE)
    
    # Run PANDA.
    nullPanda <- pandaPy(expr_file = fname, motif_file=motiffname, ppi_file=ppifname, save_tmp=FALSE)$panda
    rownames(nullPanda) <- paste(nullPanda$TF, nullPanda$Gene, sep = "__")

    # Remove file.
    unlink(fname)
    
    # Return null results.
    tfListToExclude <- c("tf1", "tf2", "tf3", "tf4", "tf5")
    geneListToExclude <- c("gene1", "gene2", "gene3", "gene1", "gene1")
    rowsToExclude <- paste(tfListToExclude, geneListToExclude, sep = "__")
    return(nullPanda[setdiff(rownames(nullPanda), rowsToExclude),"Score"])
  })
  
  # Remove files.
  unlink(motiffname)
  unlink(ppifname)

  # Return the values.
  return(unlist(nullPandas))
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
#' @returns A bipartite subnetwork in the same format as the original networks.
BuildSubnetwork <- function(geneSet, networks, alpha, hopConstraint, nullDistribution,
                            verbose = FALSE, topX = NULL){
  
  # Name edges for each network.
  networksNamed <- lapply(networks, function(network){
    rownames(network) <- paste(network$tf, network$gene, sep = "__")
    return(network)
  })

  # Paste together the networks.
  combinedNetwork <- networksNamed[[1]]
  for(i in 2:length(networksNamed)){
    combinedNetwork[,2+i] <- networksNamed[[i]]$score
  }
  
  # For each gene, find the significant edges from each hop.
  significantSubnetworks <- FindSignificantEdgesForHop(geneSet = geneSet,
                                                       combinedNetwork = combinedNetwork,
                                                       alpha = alpha,
                                                       hopConstraint = hopConstraint / 2,
                                                       nullDistribution = nullDistribution,
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
  
#' Find the subnetwork of significant edges n / 2 hops away from each gene.
#' @param geneSet A character vector of genes comprising the targets of interest.
#' @param combinedNetwork A concatenation of n PANDA-like networks with the following format:
#' tf,gene,score_net1, score_net2, ... , score_netn
#' @param alpha The significance cutoff for the statistical test.
#' @param hopConstraint The maximum number of hops to be considered for a gene.
#' @param nullDistribution The null distribution, specified as a vector of values.
#' @param verbose Whether or not to print detailed information about the run.
#' @param topX Select the X lowest significant p-values for each gene. NULL by default.
#' @returns A bipartite subnetwork in the same format as the original networks.
FindSignificantEdgesForHop <- function(geneSet, combinedNetwork, alpha, hopConstraint, nullDistribution,
                                       verbose = FALSE, topX = NULL){
  # Build the significant subnetwork for each gene, up to the hop constraint.
  uniqueGeneSet <- sort(unique(geneSet))
  geneSubnetworks <- lapply(uniqueGeneSet, function(gene){
    
    # Get all significant edges for a 1-hop subnetwork.
    if(verbose == TRUE){
      message(paste("Evaluating hop 1 for gene", gene))
    }
    subnetwork1Hop <- SignificantBreadthFirstSearch(networks = combinedNetwork, 
                                                    alpha = alpha, 
                                                    startingNodes = gene,
                                                    nodesToExclude = c(),
                                                    startFromTF = FALSE, doFDRAdjustment = FALSE,
                                                    nullDistribution, verbose = verbose,
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
                                                        alpha = alpha, 
                                                        startingNodes = startingNodes,
                                                        nodesToExclude = excludedSubset,
                                                        startFromTF = TRUE, doFDRAdjustment = FALSE,
                                                        nullDistribution = nullDistribution,
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
                                                        alpha = alpha, 
                                                        startingNodes = startingNodes,
                                                        nodesToExclude = excludedSubset,
                                                        startFromTF = FALSE, doFDRAdjustment = FALSE,
                                                        nullDistribution = nullDistribution,
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
#' @param alpha The significance cutoff for the statistical test.
#' @param startingNodes The list of nodes from which to start.
#' @param nodesToExclude The list of nodes to exclude from the search.
#' @param startFromTF Whether to start from transcription factors (TRUE) or genes (FALSE).
#' @param doFDRAdjustment Whether or not to adjust the p-values using FDR.
#' @param nullDistribution The null distribution, specified as a vector of values.
#' @param verbose Whether or not to print detailed information about the run.
#' @param topX Select the X lowest significant p-values for each gene. NULL by default.
#' @returns A bipartite subnetwork in the same format as the original networks.
SignificantBreadthFirstSearch <- function(networks, alpha, startingNodes,
                                          nodesToExclude, startFromTF, doFDRAdjustment = FALSE,
                                          nullDistribution, verbose = FALSE, topX = NULL){
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
  allEdges <- paste(tfLongList, geneLongList, sep = "__")

  # For each edge, measure its significance.
  subnetwork <- networks[c(), 3:ncol(networks)]
  if(length(allEdges) > 0){
    ourEdgeVals <- networks[allEdges, 3:ncol(networks)]
    #nullDistribution <- sample(nullDistribution, size = floor(length(nullDistribution) / length(allEdges)) * length(allEdges))
    nullEdgeVals <- t(matrix(rep(nullDistribution, length(allEdges)), ncol = length(allEdges)))
    pValues <- matrixTests::row_wilcoxon_twosample(x = ourEdgeVals, y = nullEdgeVals, alternative = "greater")$pvalue

    # Adjust the p-values.
    if(doFDRAdjustment == TRUE){
      pValues <- stats::p.adjust(pValues, method = "fdr")
    }
    whichSig <- which(pValues < alpha)

    # If topX is specified, filter again.
    if(!is.null(topX) && length(whichSig) > topX){
      whichTopX <- order(pValues[whichSig])[1:topX]
      whichSig <- whichSig[whichTopX]
    }
    
    # Return the edges meeting alpha.
    significantEdges <- allEdges[whichSig]
    subnetwork <- networks[significantEdges, c(1:2)]
    if(verbose == TRUE){
      message(paste("Retained", length(significantEdges), "out of", length(allEdges), "edges"))
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
      return(do.call(rbind, genePairSpecificHopCountSubnetwork))
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
  return(compositeSubnetworkDedup)
}

#' Plot the networks, using different colors for transcription factors, genes of interest,
#' and additional genes.
#' @param network A data frame with the following format:
#' tf,gene
#' @param genesOfInterest Which genes of interest to highlight
#' @param geneOfInterestColor Color for the genes of interest
#' @param tfColor Color for the transcription factors
#' @param otherGenesColor Color for the other genes
#' @param nodeSize Size of node
#' @param edgeWidth Width of edges
#' @param vertexLabels Which vertex labels to include. By default, none are included.
#' @param layoutBipartite Whether or not to layout as a bipartite graph.
#' @param vertexLabelSize The size of label to use for the vertex, as a fraction of the default.
#' @param vertexLabelOffset Number of pixels in the offset when plotting labels.
#' Default is TRUE.
#' @returns A bipartite plot of the network
PlotNetwork <- function(network, genesOfInterest, geneOfInterestColor = "red",
                        tfColor = "blue", otherGenesColor = "gray", nodeSize = 1,
                        edgeWidth = 0.5, vertexLabels = NA, vertexLabelSize = 0.7,
                        vertexLabelOffset = 0.5, layoutBipartite = TRUE){
  
  # Set the node attributes.
  uniqueNodes <- unique(c(network$tf, network$gene))
  nodeAttrs <- data.frame(node = uniqueNodes,
                          color = rep(otherGenesColor, length(uniqueNodes)),
                          size = rep(nodeSize, length(uniqueNodes)),
                          frame.width = rep(0, length(uniqueNodes)),
                          label.color = "black", label.cex = vertexLabelSize,
                          label.dist = vertexLabelOffset)
  rownames(nodeAttrs) <- uniqueNodes
  nodeAttrs[which(uniqueNodes %in% network$tf), "color"] <- tfColor
  nodeAttrs[which(uniqueNodes %in% genesOfInterest), "color"] <- geneOfInterestColor
  
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
    igraph::plot.igraph(graph, layout = LO, vertex.label = labels, edge.width = edgeWidth)
  }else{
    igraph::plot.igraph(graph, vertex.label = labels, edge.width = edgeWidth)
  }
 }
