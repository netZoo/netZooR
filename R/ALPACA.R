#' Main ALPACA function
#'
#' This function compares two networks and finds the sets of nodes that best characterize the change in modular structure
#' @param net.table A table of edges, with the first column representing the TFs ("from" nodes) and the second column representing the targets ("to" nodes). The third column contains the edge weights corresponding to the control or healthy network, and the fourth column contains the edge weights for the disease network or network of interest.
#' @param file.stem The folder location and title under which all results will be stored.
#' @param verbose Indicates whether the full differential modularity matrix should also be written to a file. Defaults to FALSE.
#'  modularity 
#' @return List where first element is the membership vector and second element is the contribution score of each node to its module's total differential modularity
#' @import igraph
#' @import Matrix
#' @importFrom utils write.table
#' @rawNamespace import(GOstats, except= makeGOGraph)
#' @import org.Hs.eg.db
#' @export
#' 

alpaca <- function(net.table,file.stem,verbose=F)
{
  net.table[,1] <- paste(net.table[,1],"A",sep="_")
  net.table[,2] <- paste(net.table[,2],"B",sep="_")
  
  print("Detecting communities in control network...")
  ctrl.pos <- net.table[net.table[,3]>=0,1:3]
  
  ctrl.elist <- data.frame(red=ctrl.pos[,2],blue=ctrl.pos[,1],weights=ctrl.pos[,3])
  ctrl.condor <- create.condor.object(ctrl.elist)
  ctrl.condor <- condor.cluster(ctrl.condor,project=F)
  ctrl.memb <- c(ctrl.condor$red.memb[,2],ctrl.condor$blue.memb[,2])
  names(ctrl.memb) <- c(as.character(ctrl.condor$red.memb[,1]),as.character(ctrl.condor$blue.memb[,1]))
  if (!(is.null(file.stem))) write.table(ctrl.memb, paste(c(file.stem,"_ALPACA_ctrl_memb.txt"),collapse=""),row.names=T,col.names=F,quote=F,sep="\t")
  
  pos.table <- net.table[intersect(which(net.table[,3]>=0),which(net.table[,4]>=0)),]
  pos.graph <- graph.edgelist(as.matrix(pos.table[,1:2]),directed=T)
  
  if (length(setdiff(V(pos.graph)$name,names(ctrl.memb)))>0)
  {
    uncounted <- setdiff(V(pos.graph)$name,names(ctrl.memb))
    unc.memb <- sample(1:max(ctrl.memb),length(uncounted),replace=T)
    names(unc.memb) <- uncounted
    ctrl.memb <- c(ctrl.memb,unc.memb)
  }
  
  print("Computing differential modularity matrix...")
  dwbm <- alpaca.computeDWBMmat.mscale(pos.table,ctrl.memb[V(pos.graph)$name])
  if (verbose) 
  {write.table(dwbm, paste(c(file.stem,"_DWBM.txt"),collapse=""),row.names=T,col.names=T,quote=F,sep="\t")
    write.table(rownames(dwbm),paste(c(file.stem,"_DWBM_rownames.txt"),collapse=""),row.names=T,col.names=T,quote=F,sep="\t")
    write.table(colnames(dwbm),paste(c(file.stem,"_DWBM_colnames.txt"),collapse=""),row.names=T,col.names=T,quote=F,sep="\t")
  }
  
  ntfs <- nrow(dwbm)
  ngenes <- ncol(dwbm)
  this.B <- array(0,dim=c(ntfs+ngenes,ntfs+ngenes))
  this.B[1:ntfs,(ntfs+1):(ntfs+ngenes)] <- dwbm
  this.B[(ntfs+1):(ntfs+ngenes),1:ntfs] <- t(dwbm)
  
  print("Computing differential modules...")
  louv.memb <- alpaca.genlouvain(this.B)
  names(louv.memb) <- c(rownames(dwbm),colnames(dwbm))
  
  print("Computing node scores...")
  louv.Ascores <- NULL
  louv.Bscores <- NULL
  for (i in 1:max(louv.memb))
  {
    this.comm <- names(louv.memb)[louv.memb==i]
    this.tfs <- this.comm[grep("_A$",this.comm)]
    this.genes <- this.comm[grep("_B$",this.comm)]
    if (length(this.tfs)>=1){
      if (length(this.tfs)>1){
        tf.sums <- apply(dwbm[this.tfs,this.genes],1,sum)
        gene.sums <- apply(dwbm[this.tfs,this.genes],2,sum)
        } else if (length(this.tfs)==1) {
            tf.sums <- sum(dwbm[this.tfs,this.genes])
            gene.sums <- dwbm[this.tfs,this.genes]
        } 
      this.denom <- sum(dwbm[this.tfs,this.genes])
      louv.Ascores <- c(louv.Ascores,tf.sums/this.denom)
      louv.Bscores <- c(louv.Bscores,gene.sums/this.denom)
    }
  }
  
  louv.scores <- c(louv.Ascores,louv.Bscores)
  
  if (!is.null(file.stem)) {
    write.table(cbind(names(louv.memb), as.vector(louv.memb)),paste(c(file.stem,"_ALPACA_final_memb.txt"),collapse=""),col.names=F,row.names=F,quote=F,sep="\t")
    write.table(cbind(names(louv.scores),louv.scores), paste(c(file.stem,"_ALPACA_scores.txt"),collapse=""),row.names=F,col.names=F,quote=F,sep="\t")
  }
  
  list(louv.memb,louv.scores)
}

#' Extract core target genes in differential modules
#'
#' This function outputs the top target genes in each module, ranked by their contribution to the differential modularity of the particular module in which they belong.
#' @param module.result A table of edges, with the first column representing the TFs ("from" nodes) and the second column representing the targets ("to" nodes). The third column contains the edge weights corresponding to the control or healthy network, and the fourth column contains the edge weights for the disease network or network of interest.
#' @param set.lengths The desired lengths of the top gene lists.
#'  
#' @return List with two elements. First element is a list of the top target genes in each cluster. Second element is a vector with the names of the gene sets. The names are in the format "number_length", where number is the module number label and length is the length of the gene set.
#' @import igraph
#' @import Matrix
#' @rawNamespace import(GOstats, except= makeGOGraph)
#' @import org.Hs.eg.db
#' @export
#' 

alpaca.ExtractTopGenes <- function(module.result,set.lengths)
{
  mod.memb <- module.result[[1]]
  mod.scores <- module.result[[2]]
  mod.ord <- names(mod.scores)[order(mod.scores,decreasing=T)]
  mod.Bord <- mod.ord[grep("_B$",mod.ord)]
  
  mod.top <- NULL
  mod.top.names <- NULL
  count <- 0
  for (i in 1:max(mod.memb))
  {
    mod.top.names <- c(mod.top.names,paste(i,set.lengths,sep="_"))
    this.comm <- names(mod.memb)[mod.memb==i]
    this.comm.ord <- mod.Bord[mod.Bord %in% this.comm]
    for (j in 1:length(set.lengths))
    {
      count <- count+1
      if (length(this.comm.ord)<set.lengths[j]) mod.top[[count]] <-sapply(this.comm.ord,alpaca.node.to.gene) else mod.top[[count]] <- sapply(this.comm.ord[1:set.lengths[j]],alpaca.node.to.gene)
    }
  }
  list(mod.top,mod.top.names)
}

#' The top GO term associated genes in each module
#'
#' Get all the genes in the top-scoring lists which are annotated with the enriched GO terms. Only GO terms with at least 3 genes in the overlap are included.
#' @param go.result The result of the GO term analysis (alpaca.list.to.go)
#' @param dm.top The result of extracting the top genes of the differential modules (dm.top)
#'  
#' @return A vector with strings representing gene lists, each element of the vector has the genes in that GO term and community pasted together with spaces in between.
#' @import igraph
#' @import Matrix
#' @rawNamespace import(GOstats, except= makeGOGraph)
#' @import org.Hs.eg.db
#' @import GO.db
#' @export
#' 

alpaca.GOtab.to.genes <- function(go.result,dm.top)
{
  names(dm.top[[1]]) <- dm.top[[2]]
  dm.go.genes <- NULL
  for (i in 1:nrow(go.result))
  {
    if (i %% 10 ==0) print(i)
    if (go.result[i,7]>2){
      this.go.id <- as.character(go.result[i,3])
      this.genes <- alpaca.go.to.genes(this.go.id)
      this.label <- as.character(go.result[i,1])
      comm.top <- dm.top[[1]][[this.label]]
      dm.go.genes <- c(dm.go.genes,paste(intersect(comm.top,this.genes),collapse=" "))
    } else dm.go.genes <- c(dm.go.genes,"")
  }
  dm.go.genes
}


#' Translating gene identifiers to gene symbols
#'
#' Takes a list of gene sets named using gene identifiers and converts them to a list of symbols given a user-defined annotation table.
#' @param mod.top A list of gene sets (gene identifiers)
#' @param annot.vec A vector of gene symbols with gene identifiers as the names of the vector, that defines the translation between annotations.
#'  
#' @return A list of sets of gene symbols.
#' @import igraph
#' @import Matrix
#' @rawNamespace import(GOstats, except= makeGOGraph)
#' @import org.Hs.eg.db
#' @export
#' 

alpaca.TopEnsembl.to.TopSym <- function(mod.top,annot.vec)
{
  res.sym <- NULL
  for (i in 1:length(mod.top))
  {
    x <- mod.top[[i]]
    res.sym[[i]] <- annot.vec[x[x %in% names(annot.vec)]]
  }
  res.sym
}

#' Edge subtraction method (CONDOR optimizaton)
#'
#' Takes two networks, subtracts edges individually, and then clusters the subtracted network using CONDOR. 
#' @param net.table A table of edges, with the first column representing the TFs ("from" nodes) and the second column representing the targets ("to" nodes). The third column contains the edge weights corresponding to the control or healthy network, and the fourth column contains the edge weights for the disease network or network of interest.
#' @param file.stem The folder location and title under which all results will be stored.
#'  
#' @return List where first element is the membership vector and second element is the contribution score of each node to its community's modularity in the final edge-subtracted network
#' @import igraph
#' @import Matrix
#' @importFrom utils write.table
#' @rawNamespace import(GOstats, except= makeGOGraph)
#' @import org.Hs.eg.db
#' @export
#' 

alpaca.DeltaZAnalysis <- function(net.table,file.stem)
{
  net.table[,1] <- paste(net.table[,1],"A",sep="_")
  net.table[,2] <- paste(net.table[,2],"B",sep="_")
  pos.wts <- unique(c(which(net.table[,3]>=0),which(net.table[,4]>=0)))
  net.table <- net.table[pos.wts,]
  
  net.delz <- as.numeric(net.table[,4])-as.numeric(net.table[,3])
  delz.tab <- cbind(net.table[net.delz>0,1:2],net.delz[net.delz>0])
  
  delz.elist <- data.frame(red=delz.tab[,2],blue=delz.tab[,1],weights=delz.tab[,3])
  delz.condor <- create.condor.object(delz.elist)
  delz.condor <- condor.cluster(delz.condor,project=F)
  delz.memb <- c(delz.condor$red.memb[,2],delz.condor$blue.memb[,2])
  names(delz.memb) <- c(as.character(delz.condor$red.memb[,1]),as.character(delz.condor$blue.memb[,1])) 
  
  delz.core <- condor.qscore(delz.condor)
  delz.A <- delz.core$qscores$blue.qscore
  delz.B <- delz.core$qscores$red.qscore
  delz.scores <- c(delz.A[,3],delz.B[,3])
  names(delz.scores) <- c(as.character(delz.A[,1]),as.character(delz.B[,1]))
  
  if (!is.null(file.stem)) {
    write.table(delz.memb, paste(c(file.stem,"_DelZ_memb.txt"),collapse=""),row.names=T,col.names=F,quote=F,sep="\t")
    write.table(delz.scores, paste(c(file.stem,"_DelZ_scores.txt"),collapse=""),row.names=T,col.names=F,quote=F,sep="\t")
  }
  
  list(delz.memb,delz.scores)
}

#' Edge subtraction method (Louvain optimizaton)
#'
#' Takes two networks, subtracts edges individually, and then clusters the subtracted network using Louvain method. 
#' @param net.table A table of edges, with the first column representing the TFs ("from" nodes) and the second column representing the targets ("to" nodes). The third column contains the edge weights corresponding to the control or healthy network, and the fourth column contains the edge weights for the disease network or network of interest.
#' @param file.stem The folder location and title under which all results will be stored.
#'  
#' @return List where first element is the membership vector and second element is the contribution score of each node to its community's modularity in the final edge-subtracted network
#' @import igraph
#' @import Matrix
#' @importFrom utils write.table
#' @rawNamespace import(GOstats, except= makeGOGraph)
#' @import org.Hs.eg.db
#' @export
#' 

alpaca.DeltaZAnalysis.Louvain <- function(net.table,file.stem)
{
  net.table[,1] <- paste(net.table[,1],"A",sep="_")
  net.table[,2] <- paste(net.table[,2],"B",sep="_")
  pos.wts <- unique(c(which(net.table[,3]>=0),which(net.table[,4]>=0)))
  net.table <- net.table[pos.wts,]
  
  net.delz <- as.numeric(net.table[,4])-as.numeric(net.table[,3])
  delz.tab <- data.frame(TF=net.table[net.delz>0,1],Targ=net.table[net.delz>0,2],wts=net.delz[net.delz>0])
  
  delz.res <- alpaca.WBM.louvain(delz.tab)
  delz.memb <- delz.res[[1]]
  
  delz.scores <- delz.res[[2]]
  
  if (!is.null(file.stem)) {
    write.table(delz.memb, paste(c(file.stem,"_DelZ_memb.txt"),collapse=""),row.names=T,col.names=F,quote=F,sep="\t")
    write.table(delz.scores, paste(c(file.stem,"_DelZ_scores.txt"),collapse=""),row.names=T,col.names=F,quote=F,sep="\t")
  }
  
  list(delz.memb,delz.scores)
}

#' Community comparison method (CONDOR optimizaton)
#'
#' Takes two networks, finds community structure of each one individually using CONDOR, and then ranks the nodes that show the biggest difference in their community membership. 
#' @param net.table A table of edges, with the first column representing the TFs ("from" nodes) and the second column representing the targets ("to" nodes). The third column contains the edge weights corresponding to the control or healthy network, and the fourth column contains the edge weights for the disease network or network of interest.
#'  
#' @return Vector of nodes ordered by how much they change their community membership between the two networks.
#' @import igraph
#' @import Matrix
#' @rawNamespace import(GOstats, except= makeGOGraph)
#' @import org.Hs.eg.db
#' @export
#' 

alpaca.RotationAnalysis <- function(net.table)
{
  
  net.table[,1] <- paste(net.table[,1],"A",sep="_")
  net.table[,2] <- paste(net.table[,2],"B",sep="_")
  
  net1 <- net.table[net.table[,3]>0,1:3]
  net2 <- net.table[net.table[,4]>0,c(1:2,4)]
  
  net1.elist <- data.frame(red=net1[,2],blue=net1[,1],weights=net1[,3])
  net1.condor <- create.condor.object(net1.elist)
  net1.condor <- condor.cluster(net1.condor,project=F)
  net1.memb <- c(net1.condor$red.memb[,2],net1.condor$blue.memb[,2])
  names(net1.memb) <- c(as.character(net1.condor$red.memb[,1]),as.character(net1.condor$blue.memb[,1])) 
  
  net2.elist <- data.frame(red=net2[,2],blue=net2[,1],weights=net2[,3])
  net2.condor <- create.condor.object(net2.elist)
  net2.condor <- condor.cluster(net2.condor,project=F)
  net2.memb <- c(net2.condor$red.memb[,2],net2.condor$blue.memb[,2])
  names(net2.memb) <- c(as.character(net2.condor$red.memb[,1]),as.character(net2.condor$blue.memb[,1])) 
  
  alpaca.CommunityStructureRotation(net1.memb,net2.memb)
}

#' Community comparison method (CONDOR optimizaton)
#'
#' Takes two networks, finds community structure of each one individually using a generalization of the Louvain method, and then ranks the nodes that show the biggest difference in their community membership. 
#' @param net.table A table of edges, with the first column representing the TFs ("from" nodes) and the second column representing the targets ("to" nodes). The third column contains the edge weights corresponding to the control or healthy network, and the fourth column contains the edge weights for the disease network or network of interest.
#'  
#' @return Vector of nodes ordered by how much they change their community membership between the two networks.
#' @import igraph
#' @import Matrix
#' @rawNamespace import(GOstats, except= makeGOGraph)
#' @import org.Hs.eg.db
#' @export
#' 

alpaca.RotationAnalysis.Louvain <- function(net.table)
{
  
  net.table[,1] <- paste(net.table[,1],"A",sep="_")
  net.table[,2] <- paste(net.table[,2],"B",sep="_")
  
  net1 <- net.table[net.table[,3]>0,1:3]
  net2 <- net.table[net.table[,4]>0,c(1:2,4)]
  
  net1.res <- alpaca.WBM.louvain(net1)
  net1.memb <- net1.res[[1]]
  
  net2.res <- alpaca.WBM.louvain(net2)
  net2.memb <- net2.res[[1]]
  
  alpaca.CommunityStructureRotation(net1.memb,net2.memb)
}

#' Generalized Louvain method for bipartite networks
#'
#' This function implements a generalized form of the Louvain method for weighted bipartite networks. 
#' @param net.frame A table of edges, with the first column representing the TFs ("from" nodes) and the second column representing the targets ("to" nodes). The third column contains the edge weights corresponding to the network of interest.
#'  
#' @return List where first element is the community membership vector and second element is the contribution score of each node to its community's portion of the bipartite modularity.
#' @import igraph
#' @import Matrix
#' @rawNamespace import(GOstats, except= makeGOGraph)
#' @import org.Hs.eg.db
#' @export
#' 

alpaca.WBM.louvain <- function(net.frame)
{
  wbm.mat <- alpaca.computeWBMmat(net.frame)
  ntfs <- nrow(wbm.mat)
  ngenes <- ncol(wbm.mat)
  this.B <- array(0,dim=c(ntfs+ngenes,ntfs+ngenes))
  this.B[1:ntfs,(ntfs+1):(ntfs+ngenes)] <- wbm.mat
  this.B[(ntfs+1):(ntfs+ngenes),1:ntfs] <- t(wbm.mat)
  
  print("Running alpaca.genlouvain...")
  louv.memb1 <- alpaca.genlouvain(this.B)
  louv.memb <- as.vector(louv.memb1)
  names(louv.memb) <- c(rownames(wbm.mat),colnames(wbm.mat))
  
  print("Computing node scores...")
  louv.Ascores <- NULL
  louv.Bscores <- NULL
  for (i in 1:max(louv.memb))
  {
    this.comm <- names(louv.memb)[louv.memb==i]
    this.tfs <- this.comm[grep("_A$",this.comm)]
    this.genes <- this.comm[grep("_B$",this.comm)]
    if (length(this.tfs)>1){
      tf.sums <- apply(wbm.mat[this.tfs,this.genes],1,sum)
      gene.sums <- apply(wbm.mat[this.tfs,this.genes],2,sum)} else {
        tf.sums <- sum(wbm.mat[this.tfs,this.genes])
        gene.sums <- wbm.mat[this.tfs,this.genes]
      }
    this.denom <- sum(wbm.mat[this.tfs,this.genes])
    tf.norm <- tf.sums/this.denom
    names(tf.norm) <- this.tfs
    gene.norm <- gene.sums/this.denom
    names(gene.norm) <- this.genes
    louv.Ascores <- c(louv.Ascores,tf.norm)
    louv.Bscores <- c(louv.Bscores,gene.norm)
  }
  
  louv.scores <- c(louv.Ascores,louv.Bscores)
  list(louv.memb,louv.scores)
  
}

#' Compute modularity matrix for weighted bipartite network
#'
#' This function computes the modularity matrix for a weighted bipartite network.
#' @param edge.mat A table of edges, with the first column representing the TFs ("from" nodes) and the second column representing the targets ("to" nodes). The third column contains the edge weights corresponding to the network of interest.
#'  
#' @return Modularity matrix with rows representing TFs ("from" nodes) and columns repesenting targets ("to" nodes)
#' @import igraph
#' @import Matrix
#' @rawNamespace import(GOstats, except= makeGOGraph)
#' @import org.Hs.eg.db
#' @export
#' 

alpaca.computeWBMmat <- function(edge.mat)
{
  red.nodes <- unique(edge.mat[,1])
  blue.nodes <- unique(edge.mat[,2])
  
  A = array(0,dim=c(length(red.nodes),length(blue.nodes)))
  rownames(A) <- red.nodes
  colnames(A) <- blue.nodes
  A[as.matrix(edge.mat[,1:2])] <- as.numeric(edge.mat[,3])
  ki = apply(A,1,sum)
  dj = apply(A,2,sum)
  m = sum(ki)
  
  subt.mat <- outer(ki,dj)/m
  (A-subt.mat)/m
}

#' Remove tags from gene names
#'
#' In gene regulatory networks, transcription factors can act as both "from" nodes (regulators) and "to" nodes (target genes), so the network analysis functions automatically tag the two columns to differentiate them. This function removes those tags from the gene identifiers. 
#' @param x Tagged node identifier
#'  
#' @return Untagged node name
#' @import igraph
#' @import Matrix
#' @rawNamespace import(GOstats, except= makeGOGraph)
#' @import org.Hs.eg.db
#' @export
#' 

alpaca.node.to.gene <- function(x){strsplit(x,split="_")[[1]][1]}

#' GO term enrichment for a list of gene sets
#'
#' GO term enrichment is run using the GOstats package, and corrected for multiple testing using the Benjamini-Hochberg method. 
#' @param gene.list A list consisting of vectors of genes; genes must be identified by their official gene symbols.
#' @param univ.vec A vector of  all gene symbols that were present in the original network. This set is used as the universe for running the hypergeometric test in GOstats.
#' @param comm.nums A vector of names for the gene sets in the input parameter "gene.list". These are used to create the table of final results.
#'  
#' @return A table whose rows represent enriched GO terms (p_adj<0.05) and columns describe useful properties, like the name of the GO term, the label of the gene set which is enriched in that GO term, the adjusted p-value and Odds Ratio.
#' @import igraph
#' @import Matrix
#' @rawNamespace import(GOstats, except= makeGOGraph)
#' @import org.Hs.eg.db
#' @rawNamespace import(AnnotationDbi, except= select)
#' @export
#' 

alpaca.list.to.go <- function(gene.list,univ.vec,comm.nums){
  
  xx <- AnnotationDbi::as.list(org.Hs.egSYMBOL)
  base.sym <- unlist(xx)
  base.entrez <- names(xx)
  
  univ.entrez <- base.entrez[base.sym %in% univ.vec]
  comm.entrez <- base.entrez[base.sym %in% gene.list[[1]]]
  
  #params<- new("GOHyperGParams",geneIds = comm.entrez, universeGeneIds = univ.entrez, ontology= "BP", pvalueCutoff = 1, conditional = T, testDirection="over", annotation="org.Hs.eg.db")
  
  go.tab <- NULL
  for (i in 1:length(gene.list))
  {
    print(i)
    comm.entrez <- base.entrez[base.sym %in% gene.list[[i]]]
    if (length(comm.entrez)>3){
      params<- new("GOHyperGParams",geneIds = comm.entrez, universeGeneIds = univ.entrez, ontology= "BP", pvalueCutoff = 1, conditional = T, testDirection="over", annotation="org.Hs.eg.db")
      go.res <- hyperGTest(params)
      this.df <- summary(go.res)
      this.tab <- this.df[this.df[,5]>1,]
      this.df.p <- p.adjust(this.tab[,2],method="BH")
      if (sum(this.df.p<0.05)>0)
      {
        this.tab1 <- this.tab[this.df.p<0.05,]
        this.df.p1 <- this.df.p[this.df.p<0.05]
        go.tab <- rbind(go.tab,cbind(rep(comm.nums[i],length(this.df.p1)),this.df.p1,this.tab1))
      }
    }
  }
  go.tab		
}

#' Enrichment in ranked list
#'
#' This function computes the enrichment of selected nodes in a ranked list, using Wilcoxon, Kolmogorov-Smirnov, and Fisher exact tests.  
#' @param node.ordered An ordered list of nodes (high-scoring to low-scoring).
#' @param true.pos The selected set of nodes being tested for enrichment among the ranked list.
#'  
#' @return A vector of 4 values. 1) Wilcoxon p-value, 2) KS p-value, 3) Fisher p-value, 4) Fisher odds ratio.
#' @import igraph
#' @import Matrix
#' @rawNamespace import(GOstats, except= makeGOGraph)
#' @import org.Hs.eg.db
#' @export
#' 

alpaca.TestNodeRank <- function(node.ordered,true.pos)
{
  node.ind <- rev(1:length(node.ordered))
  names(node.ind) <- node.ordered
  
  node.neg <- node.ordered[!(node.ordered %in% true.pos)]
  node.pos <- intersect(node.ordered,true.pos)
  
  wil.p <- wilcox.test(node.ind[node.pos], node.ind[node.neg],exact=F,alternative="greater",conf.int=T)$p.value 
  ks.p <- ks.test(node.ind[node.pos],node.ind[node.neg],exact=F,alternative="less")$p.value
  
  top10perc <- node.ordered[1:floor(0.1*length(node.ordered))]
  a <- length(intersect(true.pos,top10perc))
  b <- length(top10perc)-a
  c <- length(intersect(true.pos,node.ordered))-a
  d <- length(node.ordered)-a-b-c
  fish.p <- fisher.test(array(c(a,b,c,d),dim=c(2,2)),alternative="greater")$p.value
  fish.or <- (a*d)/(b*c)
  
  c(wil.p,ks.p,fish.p,fish.or)
}

#' Comparing node community membership between two networks
#'
#' This function uses the pseudo-inverse to find the optimal linear transformation mapping the community structures of two networks, then ranks nodes in the network by how much they deviate from the linear mapping.
#' @param net1.memb The community membership for Network 1.
#' @param net2.memb The community membership for Network 2.
#'  
#' @return A ranked list of nodes.
#' @import igraph
#' @import Matrix
#' @rawNamespace import(GOstats, except= makeGOGraph)
#' @examples 
#' a <- 1 #place holder
#' @import org.Hs.eg.db
#' @export
#' 

alpaca.CommunityStructureRotation <- function(net1.memb,net2.memb){
  
  net1.mat <- array(0,dim=c(length(net1.memb),max(net1.memb)))
  rownames(net1.mat) <- names(net1.memb)
  net1.mat[cbind(1:length(net1.memb),net1.memb)] <- 1
  
  net2.mat <- array(0,dim=c(length(net2.memb),max(net2.memb)))
  rownames(net2.mat) <- names(net2.memb)
  net2.mat[cbind(1:length(net2.memb),net2.memb)] <- 1
  
  common.nodes <- intersect(rownames(net1.mat),rownames(net2.mat))
  net1.comm <- net1.mat[common.nodes,]
  net2.comm <- net2.mat[common.nodes,]
  
  net1.svd <- svd(net1.comm)
  net1.inv <- net1.svd$v %*% diag(net1.svd$d^(-1)) %*% t(net1.svd$u)
  rot.mat <- as.matrix(net1.inv) %*% as.matrix(net2.comm)
  
  node.scores <- apply(cbind(net1.comm,net2.comm),1,function(x){as.vector(as.numeric(x[1:ncol(net1.comm)])) %*% rot.mat %*% as.numeric(x[(ncol(net1.comm)+1):(ncol(net1.comm)+ncol(net2.comm))])})
  
  node.scores1 <- (max(node.scores)-node.scores)/(max(node.scores)-min(node.scores))
  
  names(node.scores1)[order(node.scores1,decreasing=T)]
}

#' Differential modularity matrix
#'
#' This function computes the differential modularity matrix for weighted bipartite networks. The community structure of the healthy network is rescaled by the ratio of m (the total edge weight) of each network.
#' @param edge.mat A table of edges, with the first column representing the TFs ("from" nodes) and the second column representing the targets ("to" nodes). The third column contains the edge weights corresponding to the control or healthy network, and the fourth column contains the edge weights for the disease network or network of interest.
#' @param ctrl.memb The community membership for the control (healthy) network.
#'  
#' @return The differential modularity matrix, with rows representing "from" nodes and columns representing "to" nodes.
#' @import igraph
#' @import Matrix
#' @rawNamespace import(GOstats, except= makeGOGraph)
#' @import org.Hs.eg.db
#' @export
#' 

alpaca.computeDWBMmat.mscale <- function(edge.mat,ctrl.memb){
  
  graph.nodes <- names(ctrl.memb)
  all.A <- intersect(graph.nodes,unique(edge.mat[,1]))
  all.B <- intersect(graph.nodes,unique(edge.mat[,2]))
  
  A.cond <-  array(0,dim=c(length(all.A),length(all.B)))
  rownames(A.cond) <- all.A
  colnames(A.cond) <- all.B
  A.ctrl <- A.cond
  A.ctrl[as.matrix(edge.mat[,1:2])] <- as.numeric(edge.mat[,3])
  A.cond[as.matrix(edge.mat[,1:2])] <- as.numeric(edge.mat[,4])
  m.ctrl <- sum(A.ctrl)
  m.cond <- sum(A.cond)
  
  #cond.deg <- c(apply(A.cond,1,sum),apply(A.cond,2,sum))
  #ctrl.deg <-  c(apply(A.ctrl,1,sum),apply(A.ctrl,2,sum))
  #deg.ratio <- apply(cbind(cond.deg,ctrl.deg),1,function(x){
  #	if (x[2]>0) x[1]/x[2] else 0})
  deg.ratio <- m.cond/m.ctrl
  
  
  Dij <- A.cond
  
  for (i in 1:(max(ctrl.memb)-1))
  {
    #print(i)
    mod1.A <- intersect(names(ctrl.memb)[ctrl.memb==i],all.A)
    mod1.B <- intersect(names(ctrl.memb)[ctrl.memb==i],all.B)
    
    for (j in (i+1):max(ctrl.memb))
    {
      mod2.A <- intersect(names(ctrl.memb)[ctrl.memb==j],all.A)
      mod2.B <- intersect(names(ctrl.memb)[ctrl.memb==j],all.B)
      
      if (length(mod1.A)==1) {
        this.di1 <- sum(A.ctrl[mod1.A,mod2.B] * deg.ratio)
        this.kj1 <- A.ctrl[mod1.A,mod2.B] * deg.ratio} else if (length(mod2.B)==1) {
          this.di1 <- A.ctrl[mod1.A,mod2.B] * deg.ratio
          this.kj1 <- sum(A.ctrl[mod1.A,mod2.B] * deg.ratio)} else {
            this.di1 <- apply(A.ctrl[mod1.A,mod2.B],1,sum) * deg.ratio
            this.kj1 <- apply(A.ctrl[mod1.A,mod2.B],2,sum) * deg.ratio}
      
      if (length(mod2.A)==1) {
        this.di2 <- sum(A.ctrl[mod2.A,mod1.B] * deg.ratio)
        this.kj2 <- A.ctrl[mod2.A,mod1.B] * deg.ratio} else if (length(mod1.B)==1){
          this.di2 <- A.ctrl[mod2.A,mod1.B] * deg.ratio
          this.kj2 <- sum(A.ctrl[mod2.A,mod1.B] * deg.ratio)} else {
            this.di2 <- apply(A.ctrl[mod2.A,mod1.B],1,sum) * deg.ratio
            this.kj2 <- apply(A.ctrl[mod2.A,mod1.B],2,sum) * deg.ratio}
      
      didj1 <- outer(this.di1,this.kj1)
      didj2 <- outer(this.di2,this.kj2)
      
      this.m1 <- (sum(this.di1)+sum(this.kj1))/2+0.01
      this.m2 <- (sum(this.di2)+sum(this.kj2))/2+0.01
      
      Dij[mod1.A,mod2.B] <- Dij[mod1.A,mod2.B]-didj1/this.m1
      Dij[mod2.A,mod1.B] <- Dij[mod2.A,mod1.B]-didj2/this.m2
    }	
  }
  
  for (i in 1:max(ctrl.memb))
  {
    mod1.A <- intersect(names(ctrl.memb)[ctrl.memb==i],all.A)
    mod1.B <- intersect(names(ctrl.memb)[ctrl.memb==i],all.B)
    
    if (length(mod1.A)==1) {
      this.di1 <- sum(A.ctrl[mod1.A,mod1.B] * deg.ratio)
      this.kj1 <- A.ctrl[mod1.A,mod1.B] * deg.ratio} else if (length(mod1.B)==1){
        this.di1 <- A.ctrl[mod1.A,mod1.B] * deg.ratio
        this.kj1 <- sum(A.ctrl[mod1.A,mod1.B] * deg.ratio)} else {
          this.di1 <- apply(A.ctrl[mod1.A,mod1.B],1,sum) * deg.ratio
          this.kj1 <- apply(A.ctrl[mod1.A,mod1.B],2,sum) * deg.ratio}
    
    didj1 <- outer(this.di1,this.kj1)		
    this.m1 <- (sum(this.di1)+sum(this.kj1))/2+0.01
    
    Dij[mod1.A,mod1.B] <- Dij[mod1.A,mod1.B]-didj1/this.m1	
  }
  
  Dij/m.cond
}


#' Create alpaca.metanetwork for Louvain optimization
#'
#' Computes the "effective" adjacency matrix of a alpaca.metanetwork whose nodes represent communities in the larger input matrix.
#' @param J The modularity matrix
#' @param S The community membership vector from the previous round of agglomeration.
#'  
#' @return The differential modularity matrix, with rows representing "from" nodes and columns representing "to" nodes.
#' @import igraph
#' @import Matrix
#' @rawNamespace import(GOstats, except= makeGOGraph)
#' @import org.Hs.eg.db
#' @export
#' 

alpaca.metanetwork <- function(J,S)
{
  PP <- sparseMatrix(i=1:length(S),j=S,x=1)
  t(PP) %*% J %*% PP 
}

#' Renumbering community membership vector
#'
#' This is a helper function alpaca.genlouvain. It re-numbers the communities so that they run from 1 to N increasing through the vector.
#' @param S The community membership vector derived from the previous round of agglomeration.
#'  
#' @return The renumbered membership vector.
#' @import igraph
#' @import Matrix
#' @rawNamespace import(GOstats, except= makeGOGraph)
#' @import org.Hs.eg.db
#' @export
#' 

alpaca.tidyconfig <- function(S)
{
  T <- rep(0,length(S))
  for (i in 1:length(S))
  {
    if (T[i]==0) T[S==S[i]] <- max(T)+1
  }
  T
}

#' Generalized Louvain optimization
#'
#' This function implements the Louvain optimization scheme on a general symmetric matrix. First, nodes are all placed in separate communities, and merged iteratively according to which merge moves result in the greatest increase in the modularity sum. Note that nodes are iterated in the order of the input matrix (not randomly) so that all results are reproducible. Second, the final community membership is used to form a alpaca.metanetwork whose nodes represent communities from the prevous step, and which are connected by effective edge weights. The merging process is then repeated on the alpaca.metanetwork. These two steps are repeated until the modularity sum does not increase more than a very small tolerance factor. New
#' @param B Symmetric modularity matrix
#'  
#' @return The community membership vector
#' @import igraph
#' @import Matrix
#' @rawNamespace import(GOstats, except= makeGOGraph)
#' @import org.Hs.eg.db
#' @export
#' 

alpaca.genlouvain <- function(B)
{
  eps <- 2e-16
  if (sum(B-t(B))>0) B <- (B+t(B))/2  #force symmetric matrix
  M <- B
  
  n <- nrow(B)
  S <- t(1:n)
  
  dtot <- 0
  
  S2 <- 1:n
  Sb <- NULL
  
  n.outer <- 0
  
  while(!identical(Sb,S2))
  {
    n.outer <- n.outer+1
    if (n.outer>50) {
      print("Reached greater than 50 outer iterations.")
      return(0)
    }
    
    y <- unique(S2)
    y <- y[order(y,decreasing=F)]
    Sb <- S2
    print(paste(c("Merging",length(y),"communities"),collapse=" "))
    
    yb <- NULL
    
    G <- sparseMatrix(i=1:length(y),j=y,x=1)
    dstep <- 1
    nsteps <- 0
    
    while((!identical(yb,y)) && (dstep/dtot>2*eps))
    {
      yb <- y
      dstep <- 0
      nsteps <- nsteps+1
      print(nsteps)
      if (nsteps>50) {
        print("Reached greater than 50 inner iterations.")
        return(0)
      }
      
      #ord.i <- sample(1:nrow(M),nrow(M),replace=F)
      ord.i <- 1:nrow(M)
      
      for (i in ord.i)  
      {
        u <- unique(c(y[i],y[M[,i]>0]))
        u <- u[order(u,decreasing=F)]
        dH <- t(M[,i]) %*% G[,u]
        
        yi <- which(u==y[i])
        dH[yi] <- dH[yi] - M[i,i]
        k <- max.col(dH)
        #if (length(k)>1) k <- sample(k,1)
        
        if (dH[k]>dH[yi]){
          dtot <- dtot+dH[k]-dH[yi]
          dstep <- dstep+dH[k]-dH[yi]
          G[i,y[i]] <- 0
          G[i,u[k]] <- 1
          y[i] <- u[k]
        }
      }
    }
    
    y <- alpaca.tidyconfig(y)
    for (i in 1:length(y))
    {
      S[S==i] <- y[i]
      S2[S2==i] <- y[i]
    }
    
    if (identical(Sb,S2))
    {
      return(S)
    }
    
    M <- alpaca.metanetwork(B,S2)		
  }
}

#' Simulated networks
#'
#' This function creates a pair of networks given user-defined parameters for the modular structure of the first (healthy) network and the type of added module in the second (disease) network.
#' @param comm.sizes A two-column matrix indicating the number of "from" nodes (left column) and number of "to" nodes (right column) in each community (row).
#' @param edge.mat A matrix indicating the number of edges from the TFs in community i (rows) to target genes in community j (columns).
#' @param num.module The number of modules that will be added to simulate the disease network.
#' @param size.module A two-column matrix indicating the number of "from" and "to" nodes in each new module (row) that will be added to simulate the disease network.
#' @param dens.module A vector of length num.module, indicating the edge density of each added module.
#'  
#' @return A list with two elements. The first element is a four-column edge table of the same form that is input into the differential modularity function. The second element is a list of all the new nodes in the modules that were added to create the disease network.
#' @import igraph
#' @import Matrix
#' @rawNamespace import(GOstats, except= makeGOGraph)
#' @import org.Hs.eg.db
#' @export
#' 

alpaca.simulateNetwork <- function(comm.sizes,edge.mat,num.module,size.module,dens.module)
{
  if (length(size.module)==2) size.module <- rbind(size.module,c(0,0))
  
  num.comm <- nrow(comm.sizes)
  net1.Anodes <- paste("A",1:sum(comm.sizes[,1]),sep="")
  net1.Bnodes <- paste("B",1:sum(comm.sizes[,2]),sep="")
  A.comm <- B.comm <- NULL
  for (i in 1:num.comm)
  {
    A.comm[[i]] <- net1.Anodes[(sum(comm.sizes[1:i,1])+1-comm.sizes[i,1]):sum(comm.sizes[1:i,1])]
    B.comm[[i]] <- net1.Bnodes[(sum(comm.sizes[1:i,2])+1-comm.sizes[i,2]):sum(comm.sizes[1:i,2])]
  }
  
  A.memb <- B.memb <- NULL
  for (i in 1:length(A.comm))
  {
    A.memb <- rbind(A.memb,cbind(A.comm[[i]],rep(i,length(A.comm[[i]]))))
    B.memb <- rbind(B.memb,cbind(B.comm[[i]],rep(i,length(B.comm[[i]]))))
  }
  
  all.edges <- NULL
  weights.mat <- NULL
  for (i in 1:nrow(comm.sizes))
  {
    this.A.comm <- A.comm[[i]]
    for (j in 1:nrow(comm.sizes))
    {
      this.B.comm <- B.comm[[j]]
      this.edges <- NULL
      for (k in 1:length(this.A.comm))
        this.edges <- rbind(this.edges,cbind(rep(this.A.comm[k],length(this.B.comm)),this.B.comm))
      this.weights <- array(0,dim=c(nrow(this.edges),2))
      this.weights[sample(1:nrow(this.edges),edge.mat[i,j],replace=F),1] <- 1
      this.weights[sample(1:nrow(this.edges),edge.mat[i,j],replace=F),2] <- 1
      
      all.edges <- rbind(all.edges,this.edges)
      weights.mat <- rbind(weights.mat,this.weights)
    }
  }
  
  new.module <- NULL
  for (i in 1:num.module)
  {
    this.sizes <- size.module[i,]
    this.dens <- dens.module[i]
    new.tfs <- sample(net1.Anodes,this.sizes[1],replace=F)
    new.genes <- sample(net1.Bnodes,this.sizes[2],replace=F)
    new.edges <- intersect(which(all.edges[,1] %in% new.tfs),which(all.edges[,2] %in% new.genes))
    num.edges.add <- floor(length(new.edges)*this.dens)
    edges.to.add <- sample(new.edges,num.edges.add,replace=F)
    weights.mat[edges.to.add,2] <- 1
    new.module[[i]] <- c(new.tfs,new.genes)
  }
  
  not.both.zero <- unique(c(which(weights.mat[,1]>0),which(weights.mat[,2]>0)))
  
  sim.edge.mat <- cbind(data.frame(all.edges[not.both.zero,]),weights.mat[not.both.zero,])
  
  list(sim.edge.mat,new.module)
}

#' Map GO terms to gene symbols
#'
#' This function extracts all the gene symbols associated with a GO term and its descendants. (v1)
#' @param go.term The GO Biological Process ID (string).
#'  
#' @return A vector of all gene symbols associated with the GO term.
#' @import igraph
#' @import Matrix
#' @rawNamespace import(GOstats, except= makeGOGraph) 
#' @import org.Hs.eg.db
#' @import GO.db
#' @rawNamespace import(AnnotationDbi, except= select)
#' @export
#' 

alpaca.go.to.genes <- function(go.term)
{
  GOoffspring <- AnnotationDbi::as.list(GOBPOFFSPRING)
  kids <- GOoffspring[[go.term]]
  GOentrez <- AnnotationDbi::as.list(org.Hs.egGO2EG)
  
  go.entrez <- GOentrez[[go.term]]
  kids.entrez <- unique(unlist(GOentrez[kids]))
  all.entrez <- unique(c(go.entrez,kids.entrez))
  
  hs.sym <- org.Hs.egSYMBOL
  mapped_genes <- mappedkeys(hs.sym)
  hs.sym1 <- AnnotationDbi::as.list(hs.sym[mapped_genes])
  base.sym <- unlist(hs.sym1)
  base.entrez <- names(hs.sym1)
  
  base.sym[all.entrez]  
}




