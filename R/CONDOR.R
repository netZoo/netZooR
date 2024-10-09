#' Main clustering function for condor.
#' 
#' This function performs community structure clustering using
#' the bipartite modularity described in
#' \code{\link{condorModularityMax}}. This function uses a standard
#' (non-bipartite) community structure clustering of the uni-partite,
#' weighted projection of the original bipartite graph as an initial
#' guess for the bipartite modularity.
#' @param condor.object Output of make.condor.object. This function uses 
#' \code{condor.object$edges}
#' @param cs.method is a string to specify which unipartite community 
#' structure algorithm should be used for the seed clustering. 
#' Options are \code{LCS} (\code{\link[igraph]{multilevel.community}}), 
#' \code{LEC} (\code{\link[igraph]{leading.eigenvector.community}}), 
#' \code{FG} (\code{\link[igraph]{fastgreedy.community}}).
#' @param project Provides options for initial seeding of the bipartite 
#' modularity maximization.
#' If TRUE, the nodes in the first column of \code{condor.object$edges}
#' are projected and clustered using \code{cs.method}. If FALSE, the 
#' complete bipartite network is clustered using the unipartite clustering 
#' methods listed in \code{cs.method}.
#' @param low.memory If TRUE, uses \code{\link{condorModularityMax}}
#' instead of \code{\link{condorMatrixModularity}}. This is a slower
#' implementation of the modularity maximization, which does not store any
#' matrices in memory. Useful on a machine with low RAM. However, runtimes
#' are (much) longer.
#' @param deltaQmin convergence parameter determining the minimum required increase
#' in the modularity for each iteration. Default is min(10^{-4},1/(number of edges)),
#' with number of edges determined by \code{nrow(condor.object$edges)}. User can
#' set this parameter by passing a numeric value to deltaQmin.
#' @return \code{condor.object} with \code{\link{condorModularityMax}} output
#' included.
#' @examples 
#' r = c(1,1,1,2,2,2,3,3,3,4,4);
#' b = c(1,2,3,1,2,4,2,3,4,3,4);
#' reds <- c("Alice","Sue","Janine","Mary")
#' blues <- c("Bob","John","Ed","Hank")
#' elist <- data.frame(red=reds[r],blue=blues[b])
#' condor.object <- createCondorObject(elist)
#' condor.object <- condorCluster(condor.object)
#' @import igraph
#' @import Matrix
#' @export
#' 
condorCluster <- function(condor.object,cs.method="LCS",project=TRUE,low.memory=FALSE,deltaQmin="default"){
  
  elist <- condor.object$edges
  
  #throw exception if subgraphs are disjoint
  if (max(table(elist[, 1])) == 1) {
    stop("Unable to cluster. Subgraphs are disjoint.")
  }
  
  #extract weights, if any
  if(ncol(elist) > 2){
    weights <- elist[, 3]
  }
  else {
    weights <- rep(1, nrow(elist))
  }
  
  #make sure there's only one connected component
  g.component.test <- graph_from_data_frame(elist,directed=FALSE)
  if(!is.connected(g.component.test)){
    stop("More than one connected component detected,
         method requires only one connected component")
  }  
  
  G <- graph_from_data_frame(elist,directed=FALSE)
  
  project.weights <- weights
  #Use unipartite community structure method for first pass
  #project network into gene space to obtain initial community assignment for genes
  if(project){
    #SNP row indices
    reds = as.integer(factor(elist[,1]))
    red.names = levels(factor(elist[,1]))
    #Gene column indices
    blues = as.integer(factor(elist[,2]))
    blue.names = levels(factor(elist[,2]))
    N = max(blues)
    
    #Sparse matrix with the upper right block of the true Adjacency matrix. Notice the dimension is reds x blues
    sM = sparseMatrix(i=reds,j=blues,x=weights,dims=c(length(unique(reds)),length(unique(blues))),index1=TRUE);
    #Project into gene space, projected adjacency matrix has dim = genes x genes
    gM = t(sM) %*% sM;
    rm(sM)
    gc()
    colnames(gM) <- blue.names
    rownames(gM) <- blue.names
    G1 = graph.adjacency(gM,mode="undirected",weighted=TRUE,diag=FALSE);
    #if(clusters(G1)$no > 1){print("Warning more than one component! May cause indexing error")}
    #V(G1)$name <- sort(unique(as.vector(esub[,2])))
    #remove loops and multiple edges
    gcc.initialize = simplify(max.component(G1))
    project.weights <- edge.attributes(gcc.initialize)$weight
  }
  
  #option to treat the bipartite network as if it is unipartite
  #for community initialization only
  if(!project) {
    gcc.initialize <- G 
    blue.names = levels(factor(elist[,2]))
    #blue.indx <- V(G)$name %in% blue.names
  }
  
  if(cs.method=="LCS"){cs0 = multilevel.community(gcc.initialize, weights=project.weights)}
  if(cs.method=="LEC"){cs0 = leading.eigenvector.community(gcc.initialize, weights=project.weights)}
  if(cs.method=="FG"){cs0 = fastgreedy.community(gcc.initialize, weights=project.weights)}
  print(paste("modularity of projected graph",max(cs0$modularity)))
  
  #initial condition for blues' community membership
  if(project){
    T0 <- data.frame(as.integer(factor(blue.names)),as.vector(membership(cs0)))
  }
  else {
    blue.membs <- membership(cs0)[blue.names]
    T0 <- data.frame(as.integer(factor(blue.names)),blue.membs)
  }
  #run bipartite modularity maximization using initial assignments T0
  if(low.memory){
    condor.object <- condorModularityMax(condor.object,T0=T0,weights=weights,deltaQmin=deltaQmin)
  }
  if(!low.memory){
    condor.object <- condorMatrixModularity(condor.object,T0=T0,weights=weights,deltaQmin=deltaQmin)
  }
  
  
  return(condor.object)
  }


#' Compare qscore distribution of a subset of nodes to all other nodes.
#' 
#' Compute one-sided KS and wilcox tests to determine if a subset of nodes
#' has a stochastically larger qscore distribution.
#' @param test_nodes is a list containing the subset of nodes (of one node class
#' --blue or red--only) to be tested
#' @param q is a two column data frame containing the node names in the 
#' first column and the q-scores in the second column.
#' @param perm if TRUE, run permutation tests. Else, run 
#' \code{\link[stats]{ks.test}} and \code{\link[stats]{wilcox.test}} only.
#' @param plot.hist if TRUE, produces two histograms of test statistics 
#' from permutation tests, one for KS and one for 
#' wilcoxon and a red dot for true labeling. Only works if perm=TRUE.
#' @param nsamp Number of permutation tests to run 
#' @return if \code{perm=FALSE}, the analytical p-values from 
#' \code{\link[stats]{ks.test}} and \code{\link[stats]{wilcox.test}}
#' @return if \code{perm=TRUE}, the permutation p-values are provided in 
#' addition to the analytical values.
#' @note \code{\link[stats]{ks.test}} and \code{\link[stats]{wilcox.test}}
#' will throw warnings due to the presence of ties, so the p-values will be 
#' approximate. See those functions for further details.
#' @rawNamespace import(stats, except= c(cov2cor,decompose,toeplitz,lowess,update,spectrum))
#' @examples
#' r = c(1,1,1,2,2,2,3,3,3,4,4);
#' b = c(1,2,3,1,2,4,2,3,4,3,4);
#' reds <- c("Alice","Sue","Janine","Mary")
#' blues <- c("Bob","John","Ed","Hank")
#' elist <- data.frame(red=reds[r],blue=blues[b])
#' condor.object <- createCondorObject(elist)
#' condor.object <- condorCluster(condor.object)
#' condor.object <- condorQscore(condor.object) 
#' q_in <- condor.object$qscores$red.qscore
#' out <- condorCoreEnrich(c("Alice","Mary"),q=q_in,perm=TRUE,plot.hist=TRUE)
#' @export
#' 
condorCoreEnrich = function(test_nodes,q,perm=FALSE,plot.hist=FALSE,nsamp=1000){
  qtest <- q[q[,1] %in% test_nodes,3]
  #
  qall <- q[,3]
  qnot_test <- q[!(q[,1] %in% test_nodes),3]
  ks_out <- ks.test(qtest,qnot_test,exact=FALSE,alternative="less")
  
  w_out <- wilcox.test(qtest,qnot_test,exact=FALSE,alternative="greater")
  if(perm){
    qnull <- q[!(q[,1] %in% test_nodes),3]
    ks_true <- ks.test(qtest,qnull,exact=FALSE,alternative="less")$statistic
    ks_rand <- ks.permute(qtest,qnull,nsamp=nsamp)
    pv_permuted <- perm.pval(ks_true,ks_rand)
    #now for the wilcox test
    w_true <- wilcox.test(qtest,qnull,exact=FALSE,alternative="greater")$statistic
    w_rand <- wilcox.permute(qtest,qnull,nsamp=nsamp)
    pvw_permuted <- perm.pval(w_true,w_rand)
    analytical.pvals = cbind.data.frame(ks.pvalue=ks_out$p.value,wilcox.pvalue=w_out$p.value)
    perm.pvals = data.frame(ks.perm=pv_permuted,wilcox.perm=pvw_permuted)
    perm.scores <- data.frame(ks.rand=ks_rand,wilcox.rand=w_rand)
    true.scores <- data.frame(ks.true=ks_true,wilcox.true=w_true)
    out <- list(analytical.pvals=analytical.pvals,perm.pvals=perm.pvals,
                perm.scores=perm.scores,true.scores=true.scores)
  }
  if(!perm){out=cbind.data.frame(ks.pvalue=ks_out$p.value,wilcox.pvalue=w_out$p.value)}
  if(perm && plot.hist){plot.enrich.hist(qik_enrich_out=out)}
  
  return(out)
}

ks.permute = function(A,B,nsamp=1000){
  n = length(A)
  ks_out = vector()
  all <- c(A,B)
  for(i in seq_len(nsamp)){
    #randomly choose n values from all, assign the rest to B_rand
    aind = sample(seq_len(length(all)),n)
    A_rand = all[aind]
    B_rand = all[-aind]
    #note: ks.test will throw a warning because of ties. This is expected.
    ks_temp = ks.test(A_rand,B_rand,exact=FALSE,alternative="less")$statistic
    ks_out <- c(ks_out,ks_temp)
  }
  return(ks_out)
}

wilcox.permute= function(A,B,nsamp=1000){
  n = length(A)
  w_out = vector()
  all <- c(A,B)
  for(i in seq_len(nsamp)){
    #randomly choose n values from all, assign the rest to B_rand
    aind = sample(seq_len(length(all)),n)
    A_rand = all[aind]
    B_rand = all[-aind]
    #note: wilcoxon.test will throw a warning because of ties. This is expected.
    w_temp = wilcox.test(A_rand,B_rand,exact=FALSE,alternative="greater")$statistic
    w_out <- c(w_out,w_temp)
  }
  return(w_out)
}

perm.pval = function(stat_true,stat_random){
  #definition of permutation pvalue from
  #T. Knijnenburg et. al. Bioinformatics 2009 vol. 25, eqn (2)
  x0 = stat_true
  y_n = stat_random
  N = length(stat_random)
  pval = (1 + sum(y_n >= x0))/N
  return(pval)
}

plot.enrich.hist = function(qik_enrich_out,ks=TRUE,wilcoxon=TRUE,...){
  if(ks & wilcoxon){
    par(mfrow=c(1,2))
    par(mar=c(3,4,3,1),oma=c(1,1,1,1))
  }
  if(wilcoxon){
    w_rand <- qik_enrich_out$perm.scores$wilcox.rand
    w_true <- qik_enrich_out$true.scores$wilcox.true
    pvw_permuted <- qik_enrich_out$perm.pvals$wilcox.perm
    h1 <- hist(w_rand,plot=FALSE)[c("breaks","counts")]
    hist(w_rand,col="steelblue",xlim=c(min(h1$breaks),1.1*max(c(w_rand,w_true))),
         xlab="",ylab="",xaxt="n",yaxt="n",main="")
    rect(min(c(h1$breaks,w_true)),0,max(c(h1$breaks,w_true))+0.1,1.03*max(h1$counts),lwd=2)
    axis(1, at = pretty(w_rand), pos = 0,font=2,lwd.ticks=2)
    axis(2, at = pretty(h1$counts), pos = min(h1$breaks),lwd=1,lwd.ticks=2,font=2)
    points(w_true,5,cex=2,col="red",pch=19)
    mtext("Wilcoxon Test Score",side=1,line=1.75,font=2)
    mtext("Counts",side=2,line=1.5,font=2)
    mtext(paste0("Perm. P-value = ",
                 as.expression(signif(pvw_permuted,digits=3))),side=3,font=2)
  }
  if(ks){
    ks_rand <- qik_enrich_out$perm.scores$ks.rand
    ks_true <- qik_enrich_out$true.scores$ks.true
    pv_permuted <- qik_enrich_out$perm.pvals$ks.perm
    h1 <- hist(ks_rand,plot=FALSE)[c("breaks","counts")]
    hist(ks_rand,col="steelblue",xlim=c(0,max(c(ks_rand,ks_true))+0.1),
         xlab="",ylab="",xaxt="n",yaxt="n",main="")
    rect(0,0,max(c(ks_rand,ks_true))+0.1,1.03*max(h1$counts),lwd=2)
    axis(1, at = pretty(c(ks_rand,ks_true)),pos = 0,font=2,lwd.ticks=2)
    axis(2, at = pretty(h1$counts), pos = 0,lwd=0,lwd.ticks=2,font=2)
    points(ks_true,5,cex=2,col="red",pch=19)
    mtext("KS Test Score",side=1,line=1.75,font=2)
    mtext("Counts",side=2,line=1.5,font=2)
    mtext(paste0("Perm. P-value = ",
                 as.expression(signif(pv_permuted,digits=3))),side=3,font=2)
  }
}


#' Iteratively maximize bipartite modularity.
#' 
#' This function is based on the bipartite modularity
#' as defined in "Modularity and community detection in bipartite networks"
#' by Michael J. Barber, Phys. Rev. E 76, 066102 (2007)
#' This function uses a slightly different implementation from the paper. It 
#' does not use the "adaptive BRIM" method for identifying the number of 
#' modules. Rather, it simply continues to iterate until the difference in 
#' modularity between iterations is less that 10^-4. Starting from a random 
#' initial condition, this could take some time. Use 
#' \code{\link{condorCluster}} for quicker runtimes and likely better 
#' clustering, it initializes the blue 
#' node memberships by projecting the blue nodes into a unipartite "blue" 
#' network and then identify communities in that network using a standard 
#' unipartite community detection algorithm run on the projected network.
#' See \code{\link{condorCluster}} for more details on that.
#' This function loads the entire adjacency matrix in memory, so if your
#' network has more than ~50,000 nodes, you may want to use
#' \code{\link{condorModularityMax}}, which is slower, but does not store
#' the matrices in memory. Or, of course, you could move to a larger machine.
#' @param condor.object is a list created by 
#' \code{\link{createCondorObject}}. \code{condor.object$edges} must 
#' contain the edges in the giant connected component of a bipartite network 
#' @param T0 is a two column data.frame with the initial community 
#' assignment for each "blue" node, assuming there are more reds than blues, 
#' though this is not strictly necessary. The first column contains the 
#' node name, the second column the community assignment.
#' @param weights edgeweights for each edge in \code{edgelist}.
#' @param deltaQmin convergence parameter determining the minimum required increase
#' in the modularity for each iteration. Default is min(10^-4,1/(number of edges)),
#' with number of edges determined by \code{nrow(condor.object$edges)}. User can
#' set this parameter by passing a numeric value to deltaQmin.
#' @return Qcoms data.frame with modularity of each community.
#' @return modularity modularity value after each iteration.
#' @return red.memb community membership of the red nodes
#' @return blue.memb community membership of the blue.nodes
#' @examples 
#' r = c(1,1,1,2,2,2,3,3,3,4,4);
#' b = c(1,2,3,1,2,4,2,3,4,3,4);
#' reds <- c("Alice","Sue","Janine","Mary")
#' blues <- c("Bob","John","Ed","Hank")
#' elist <- data.frame(red=reds[r],blue=blues[b])
#' condor.object <- createCondorObject(elist)
#' #randomly assign blues to their own community
#' T0 <- data.frame(nodes=blues,coms=seq_len(4))
#' condor.object <- condorMatrixModularity(condor.object,T0=T0)
#' @import nnet
#' @export
#' 
condorMatrixModularity = function(condor.object,T0=cbind(seq_len(q),rep(1,q)),weights=1,deltaQmin="default"){
  
  #assign convergence parameter
  if(deltaQmin == "default"){
    #number of edges
    m = nrow(condor.object$edges)
    deltaQmin <- min(10^-4,1/m)
  }
  if(deltaQmin != "default" & !is.numeric(deltaQmin)){
    stop("deltaQmin must be either 'default' or a numeric value")
  }
  
  
  #define a local function
  machine.which.max = function(X){
    thresh <- .Machine$double.eps
    X[abs(X) <= thresh] <- 0
    true_max <- which.max(X)
    return(true_max)
  }
  machine.max = function(X){
    thresh <- .Machine$double.eps
    X[abs(X) <= thresh] <- 0
    true_max <- max(X)
    return(true_max)
  }
  machine.which.is.max = function(X){
    thresh <- .Machine$double.eps
    X[abs(X) <= thresh] <- 0
    true_max <- which.is.max(X)
    return(true_max)
  }
  #Convert the edgelist to a sparseMatrix object
  esub <- condor.object$edges
  #make sure there's only one connected component
  g.component.test <- graph_from_data_frame(esub,directed=FALSE)
  if(!is.connected(g.component.test)){
    stop("More than one connected component,
         method requires only one connected component")
  }
  reds <- as.integer(factor(esub[,1]))
  red.names <- levels(factor(esub[,1]))
  blues <- as.integer(factor(esub[,2]))
  blue.names <- levels(factor(esub[,2]))
  #ensure that nrows > ncols
  if(length(red.names) < length(blue.names)){
    stop("Adjacency matrix dimension mismatch: This code requires nrows > ncols")
  }
  
  #The upper right block of the true Adjacency matrix. notices the dimension is reds x blues
  #A = sparseMatrix(i=reds,j=blues,x=weights,dims=c(length(unique(reds)),length(unique(blues))),index1=TRUE);
  A = matrix(0,nrow=length(unique(reds)),ncol=length(unique(blues)))
  edges = cbind(reds,blues)
  A[edges] <- weights
  rownames(A) <- red.names
  colnames(A) <- blue.names
  
  p = nrow(A)
  q = ncol(A)
  N = p+q
  
  #initialize community assignments for blue nodes.
  T0[,1] = as.integer(factor(T0[,1]))
  T0 = T0[order(T0[,1]),]
  Tind <- data.matrix(T0)
  Rind = data.matrix(cbind(seq_len(p),rep(0,length=p)))
  cs = sort(unique(Tind[,2]))
  
  #variables to track modularity changes after each loop
  Qhist <- vector();
  Qnow  <- 0
  deltaQ <- 1
  
  
  #Sparse matrix with the upper right block of the true Adjacency matrix. Notice the dimension is reds x blues
  ki = rowSums(A)
  dj = colSums(A)
  m = sum(ki) # m = sum(dj) too
  
  #Make B tilde, Btilde = Aij - (ki*dj)/m
  Btilde = -(ki %o% dj)/m ###
  Btilde[edges] = weights + Btilde[edges]
  #initialize Tm
  #Tm = sparseMatrix(i=Tind[,1],j=Tind[,2],x=1,index1=TRUE)
  Tm = matrix(0,nrow=q,ncol=max(Tind[,2]))
  Tm[Tind] <- 1 
  #Begin iterations?
  
  #______________________________________
  iter=1
  while(deltaQ > deltaQmin){
    
    ### Step 1, assign red nodes
    # T tilde = Btilde %*% T
    Ttilde = Btilde %*% Tm
    
    #Find first max for all rows, update R
    Rind[,2] <- unname(apply(Ttilde, 1, machine.which.max))
    
    #Special condition if all nodes are stuck in one large community,
    #if TRUE, randomly assign two nodes to new communities.
    if(iter > 1 && length(unique(c(Rind[,2],Tind[,2]))) == 1){
      random_nodes <- sample(seq_len(nrow(Rind)),2)
      Rind[random_nodes,2] <- max(c(Rind[,2],Tind[,2])) + seq_len(2)
    }
    
    #Check to see if new communities should be made
    negative_contribution <- unname(apply(Ttilde, 1,machine.max)) < 0
    if( any( negative_contribution )){
      #add new communities
      num_new_com <- sum(negative_contribution)
      cs <- length(unique(c(Tind[,2],Rind[,2])))
      Rind[negative_contribution,2] <- (cs + 1):(cs + num_new_com)
    }
    
    #Rm = sparseMatrix(i=Rind[,1],j=Rind[,2],x=1,index1=TRUE)
    Rm = matrix(0,nrow=p,ncol=max(Rind[,2]))
    Rm[data.matrix(Rind)] <- 1
    ### Step 2, assign blue nodes
    
    # R tilde = transpose(Btilde) %*% R
    Rtilde = crossprod(Btilde,Rm)
    
    #Find first max for all rows, update T
    Tind[,2] <- unname(apply(Rtilde, 1, machine.which.max))
    
    #Check to see if new communities should be made
    negative_contribution <- unname(apply(Rtilde, 1,machine.max)) < 0
    if( any( negative_contribution )){
      #add new communities
      num_new_com <- sum(negative_contribution)
      cs <- length(unique(c(Tind[,2],Rind[,2])))
      Tind[negative_contribution,2] <- (cs + 1):(cs + num_new_com)
    }
    #Tm dimensions, note the extra empty community.
    #Tm = sparseMatrix(i=Tind[,1],j=Tind[,2],x=1,index1=TRUE) 
    Tm = matrix(0,nrow=q,ncol=max(Tind[,2]))
    Tm[data.matrix(Tind)] <- 1 
    
    Qthen <- Qnow
    #replace with diag(crossprod(T,BTR))/m
    Qcom <- diag(crossprod(Rm,Btilde %*% Tm))/m
    Qnow <- sum(Qcom)
    if(abs(Qnow) < .Machine$double.eps){Qnow <- 0}
    Qhist = c(Qhist,Qnow)
    
    print(paste("Q =",Qnow,sep=" "))
    if(Qnow != 0){
      deltaQ = Qnow - Qthen
    }
    iter=iter+1
  }
  
  #__________end_while_____________
  
  cs <- seq_len(ncol(Tm))
  if(any(sort(unique(Tind[,2])) != sort(unique(Rind[,2])))){stop("number of red and blue communities unequal")}
  #drop empty communities
  qcom_temp <- cbind(Qcom,cs)
  qcom_out <- qcom_temp[abs(Qcom) > .Machine$double.eps,]
  
  #if communities were dropped, relabel so community labels can function
  #as row/column indices in the modularity matrix, B_ij.
  if(nrow(qcom_out) < nrow(qcom_temp)){
    qcom_out[,2] <- as.integer(factor(qcom_out[,2]))
    Rind[,2] <- as.integer(factor(Rind[,2]))
    Tind[,2] <- as.integer(factor(Tind[,2]))
  }
  colnames(qcom_out) <- c("Qcom","community")
  condor.object$Qcoms = qcom_out
  condor.object$modularity=Qhist
  condor.object$red.memb=data.frame(red.names=red.names[Rind[,1]],com=Rind[,2])
  condor.object$blue.memb=data.frame(blue.names=blue.names[Tind[,1]],com=Tind[,2])
  
  return(condor.object)
  }


#' Iteratively maximize bipartite modularity.
#' 
#' This function is based on the bipartite modularity
#' as defined in "Modularity and community detection in bipartite networks"
#' by Michael J. Barber, Phys. Rev. E 76, 066102 (2007)
#' This function uses a slightly different implementation from the paper. It 
#' does not use the "adaptive BRIM" method for identifying the number of 
#' modules. Rather, it simply continues to iterate until the difference in 
#' modularity between iterations is less that 10^-4. Starting from a random 
#' initial condition, this could take some time. Use 
#' \code{\link{condorCluster}} for quicker runtimes and likely better 
#' clustering, it initializes the blue 
#' node memberships by projecting the blue nodes into a unipartite "blue" 
#' network and then identify communities in that network using a standard 
#' unipartite community detection algorithm run on the projected network.
#' See \code{\link{condorCluster}} for more details that.
#' @param condor.object is a list created by 
#' \code{\link{createCondorObject}}. \code{condor.object$edges} must 
#' contain the edges in the giant connected component of a bipartite network 
#' @param T0 is a two column data.frame with the initial community 
#' assignment for each "blue" node, assuming there are more reds than blues, 
#' though this is not strictly necessary. The first column contains the 
#' node name, the second column the community assignment.
#' @param weights edgeweights for each edge in \code{edgelist}.
#' @param deltaQmin convergence parameter determining the minimum required increase
#' in the modularity for each iteration. Default is min(10^-4,1/(number of edges)),
#' with number of edges determined by \code{nrow(condor.object$edges)}. User can
#' set this parameter by passing a numeric value to deltaQmin.
#' @return Qcoms data.frame with modularity of each community. 
#' @return modularity modularity value after each iteration.
#' @return red.memb community membership of the red nodes
#' @return blue.memb community membership of the blue.nodes
#' @examples 
#' r = c(1,1,1,2,2,2,3,3,3,4,4);
#' b = c(1,2,3,1,2,4,2,3,4,3,4);
#' reds <- c("Alice","Sue","Janine","Mary")
#' blues <- c("Bob","John","Ed","Hank")
#' elist <- data.frame(red=reds[r],blue=blues[b])
#' condor.object <- createCondorObject(elist)
#' #randomly assign blues to their own community
#' T0 <- data.frame(nodes=blues,coms=1)
#' condor.object <- condorModularityMax(condor.object,T0=T0)
#' @import Matrix
#' @import nnet
#' @export
#' 
condorModularityMax = function(condor.object,T0=cbind(seq_len(q),rep(1,q)),weights=1,deltaQmin="default"){
  
  #assign convergence parameter
  if(deltaQmin == "default"){
    #number of edges
    m = nrow(condor.object$edges)
    deltaQmin <- min(10^-4,1/m)
  }
  if(deltaQmin != "default" & !is.numeric(deltaQmin)){
    stop("deltaQmin must be either 'default' or a numeric value")
  }
  
  
  
  #Convert the edgelist to a sparseMatrix object
  esub <- condor.object$edges
  #make sure there's only one connected component
  g.component.test <- graph_from_data_frame(esub,directed=FALSE)
  if(!is.connected(g.component.test)){
    stop("More than one connected component,
         method requires only one connected component")
  }
  reds <- as.integer(factor(esub[,1]))
  red.names <- levels(factor(esub[,1]))
  blues <- as.integer(factor(esub[,2]))
  blue.names <- levels(factor(esub[,2]))
  #ensure that nrows > ncols
  if(length(red.names) < length(blue.names)){
    stop("Adjacency matrix dimension mismatch: This code requires nrows > ncols")
  }
  
  #Sparse matrix with the upper right block of the true Adjacency matrix. notices the dimension is reds x blues
  A = sparseMatrix(i=reds,j=blues,x=weights,dims=c(length(unique(reds)),length(unique(blues))),index1=TRUE);
  rownames(A) <- red.names
  colnames(A) <- blue.names
  
  p = nrow(A)
  q = ncol(A)
  N = p+q
  
  #Sparse matrix with the upper right block of the true Adjacency matrix. Notice the dimension is reds x blues
  ki = rowSums(A)
  dj = colSums(A)
  m = sum(ki) # m = sum(dj) too
  
  #initialize community assignments for red and blue nodes.
  T0[,1] = as.integer(factor(T0[,1]))
  Tmat <- T0
  R = cbind(seq_len(p),rep(0,length=p))
  cs = sort(unique(Tmat[,2]))
  #variables to track modularity changes after each loop
  Qhist <- vector();
  Qnow  <- 0
  deltaQ <- 1
  #______________________________________
  iter=1
  while(deltaQ > deltaQmin){
    btr <- BTR <- bt <- BT <- vector();
    #calculate T tilde
    for(i in seq_len(p)){
      if(i %% 2500 == 0){print(sprintf("%s%% through iteration %s",round(i/p*100,digits=1),iter))}
      #find the optimal community for node i
      bt <- rep(0,length(cs))
      for(k in cs){
        ind <- Tmat[,2] == k
        if(any(ind)){
          #bt[k] = sum(A[i,ind] - (ki[i]*dj[ind])/m)
          bt[k] = sum((A[i,] - (ki[i]*dj)/m)[ind])
        }
      }
      #note that which.max returns the FIRST max if more than one max
      h = which.max(bt)
      if(bt[h] < 0){
        print("making new comm")
        R[i,2] <- max(Tmat[,2])+1
        cs <- sort(c(cs,R[i,2]))
      }
      #if(length(h) > 1){h <- sample(h,1)}
      if(bt[h] >= 0){
        R[i,2] <- h # assign blue vertex i to comm k such that Q is maximized
        bt[-h] <- 0 # BTR is zero if i is not in k (see definition of Q)
      }
      
      #BT <- rbind(BT,bt)
    }
    #calculate R tilde, i.e., B_transpose * R
    for(j in seq_len(q))
    {
      #initialize jth row of (B_transpose) * R
      btr = rep(0,length(cs))
      #calculate Q for assigning j to different communities
      for(k in cs)
      {
        #if node j is in community k, else BTR[j,k] = 0
        ind <- R[,2] == k
        if(any(ind)){
          #btr[k] = sum(A[ind,j]-(ki[ind]*dj[j])/m)
          btr[k] = sum((A[,j]-(ki*dj[j])/m)[ind])
        }
      }
      
      g = which.max(btr)
      #if there is no comm. assignment to increase modularity, make
      #a new community with that node only
      if(btr[g] < 0)
      {
        print("making new comm")
        Tmat[j,2] <- max(R[,2])+1
        btr <- rep(0,length(cs)+1)
        cs <- c(cs,Tmat[j,2])
        btr[length(cs)] <- sum(A[,j]-(ki*dj[j])/m)
        print(btr)
      }
      if(btr[g] >= 0)
      {
        Tmat[j,2] <- g
        btr[-g] <- 0
      }
      
      #add another column to BTR if a new community is added
      if( !is.vector(BTR) && dim(BTR)[2] < length(btr)){BTR <- cbind(BTR,0) }
      BTR <- rbind(BTR,btr)
    }
    
    Tt =  t(sparseMatrix(i=Tmat[,1],j=Tmat[,2],x=1,dims=c(q,length(cs)),index1=TRUE))
    Qthen <- Qnow
    Qcom <- diag(Tt %*% BTR)/m
    Qnow <- sum(Qcom)
    Qhist = c(Qhist,Qnow)
    
    print(paste("Q =",Qnow,sep=" "))
    if(Qnow != 0){
      deltaQ = Qnow - Qthen
    }
    iter=iter+1
  }
  #__________end_while_____________
  qcom_temp <- cbind(Qcom,sort(unique(cs)))
  #drop empty communities
  qcom_out <- qcom_temp[qcom_temp[,1] > 0,]
  #if communities were dropped, relabel so community labels can function
  #as row/column indices in the modularity matrix, B_ij.
  if(nrow(qcom_out) < nrow(qcom_temp)){
    qcom_out[,2] <- as.integer(factor(qcom_out[,2]))
    R[,2] <- as.integer(factor(R[,2]))
    Tmat[,2] <- as.integer(factor(Tmat[,2]))
  }
  colnames(qcom_out) <- c("Qcom","community")
  condor.object$Qcoms = qcom_out
  condor.object$modularity=Qhist
  condor.object$red.memb=data.frame(red.names=red.names[R[,1]],com=R[,2])
  condor.object$blue.memb=data.frame(blue.names=blue.names[Tmat[,1]],com=Tmat[,2])
  
  return(condor.object)
  }


#' Plot adjacency matrix with links grouped and colored by community
#' 
#' This function will generate the network link 'heatmap' with colored dots
#' representing within-community links and black dots between-community 
#' links
#' @param condor.object output of either \code{\link{condorCluster}} or 
#' \code{\link{condorModularityMax}}
#' @param color_list vector of colors accepted by \code{col} inside the 
#' \code{\link[graphics]{plot}} function. There must be as many colors as 
#' communities.
#' @param point.size passed to \code{cex} in the 
#' \code{\link[graphics]{plot}}
#' @param xlab x axis label
#' @param ylab y axis label
#' @return produces a \code{\link[graphics]{plot}} output.
#' @references \url{http://tools.medialab.sciences-po.fr/iwanthue/} for 
#'  a nice color generator at 
#' @note For the condor paper \url{http://arxiv.org/abs/1509.02816}, I used
#'   35 colors from the "Tarnish" palette with "hard" clustering
#' @examples
#' r = c(1,1,1,2,2,2,3,3,3,4,4);
#' b = c(1,2,3,1,2,4,2,3,4,3,4);
#' reds <- c("Alice","Sue","Janine","Mary")
#' blues <- c("Bob","John","Ed","Hank")
#' elist <- data.frame(red=reds[r],blue=blues[b])
#' condor.object <- createCondorObject(elist)
#' condor.object <- condorCluster(condor.object)
#' condorPlotCommunities(condor.object,
#' color_list=c("darkgreen","darkorange"),point.size=2,
#' xlab="Women",ylab="Men")
#' @rawNamespace import(data.table, except= c(dcast,melt))
#' @importFrom graphics axis box hist mtext par plot points rect
#' @export
#'  
condorPlotCommunities = function(condor.object,color_list,point.size=0.01,
                                   xlab="SNP",ylab="Gene"){
  
  dt0 <- data.table(condor.object$edges)
  setnames(dt0,seq_len(2),c("SNP","gene"))
  dt1 <- data.table(condor.object$red.memb)
  setnames(dt1,c("SNP","red.memb"))
  dt2 <- data.table(condor.object$blue.memb)
  setnames(dt2,c("gene","blue.memb"))
  dt3 <- merge(dt0,dt1,by="SNP",all.x=TRUE)    
  eqtl_object <- merge(dt3,dt2,by="gene",all.x=TRUE)
  setkey(eqtl_object,"SNP")
  eqtl_all <- data.table(eqtl_object[!is.na(SNP)])
  #this groups red and blue nodes in the same community. very important
  eqtl_block <- eqtl_all[blue.memb==red.memb]
  # coerce non-factor inpout
  eqtl_block$gene = as.factor(eqtl_block$gene)
  eqtl_block$SNP = as.factor(eqtl_block$SNP)
  eqtl_all$gene = as.factor(eqtl_all$gene)
  eqtl_all$SNP = as.factor(eqtl_all$SNP)
  
  #setkeyv(eqtl1,c("SNP","blue.memb","gene","red.memb"))
  if(nlevels(eqtl_block$SNP) != length(unique(eqtl_block$SNP))){
    print("warning: empty levels in SNP column. This may cause silent issues with plotting.")}
  #select all links that connect nodes in the same community
  setkey(eqtl_block,"blue.memb","red.memb")
  #make new index for each node that will correspond to it's row/col number
  red_tmp <- data.table(rindx=seq_len(nlevels(eqtl_block$SNP)),SNP=unique(eqtl_block$SNP))
  red_indx <- merge(red_tmp,unique(eqtl_block,by="SNP")[,c("SNP","red.memb"),with=FALSE],by="SNP")
  red_indx[,red.com.size:=length(unique(SNP)),by=red.memb]
  red_indx[red.com.size > 1,rindx:=sample(x=rindx),by=red.memb][,red.memb:=NULL,]
  setkey(red_indx,"SNP")
  blue_tmp <- data.table(bindx=seq_len(nlevels(eqtl_block$gene)),gene=unique(eqtl_block$gene))
  blue_indx <- merge(blue_tmp,unique(eqtl_block,by="gene")[,c("gene","blue.memb"),with=FALSE],by="gene")
  #shuffle nodes within each community to make density homogeneous
  blue_indx[,blue.com.size:=length(unique(gene)),by=blue.memb]
  blue_indx[blue.com.size > 1,bindx:=sample(x=bindx),by=blue.memb][,blue.memb:=NULL,]
  setkey(blue_indx,"gene")
  
  if(dim(red_indx)[1] != nlevels(eqtl_all$SNP) && dim(blue_indx)[1] != nlevels(eqtl_all$gene)){
    print("Warning! not all nodes in block!")
  }
  
  #in the unlikely event a node is only connected to nodes in OTHER comms
  #if(nlevels(eqtl_all$SNP) != nlevels(eqtl_all$SNP)){
  #  tmp = setdiff(levels(eqtl_all$SNP),levels(eqtl_all$SNP))
  m1 <- merge(eqtl_all,red_indx,by="SNP",all=TRUE)
  #setkey(m1,"gene")
  m2 <- merge(m1,blue_indx,by="gene",all=TRUE)
  
  #pdf("Community_structure_matrix.pdf",height=7,width=12)
  #setEPS()
  #postscript(paste0(figure_out,".eps"),height=7,width=15)#,width=720,height=480,res=300,pointsize=3)
  par(mar=c(3,3,3,0.5)+0.1)
  #plot links that connect nodes in different communities
  m2[red.memb != blue.memb][plot(rindx,bindx,cex=point.size,xaxt="n",yaxt="n",
                                 xaxs="i",yaxs="i",ylim=c(0,max(m2$bindx)+1),
                                 xlim=c(0,max(m2$rindx)+1),xlab="",ylab="",pch=19)]
  #plot links that connect nodes in same communities
  m2[red.memb==blue.memb][points(rindx,bindx,cex=point.size,pch=19,
                                 col=color_list[red.memb])]  
  box(lwd=2)
  mtext(xlab,side=3,font=2,cex=2.5,padj=-0.25)
  mtext(ylab,side=2,font=2,cex=2.5,padj=-0.5)
  
  ## Add community labels to top 
  cs <- cumsum(rle(sort(m2[!duplicated(SNP)]$red.memb))$lengths)
  lens <- rle(sort(m2[!duplicated(SNP)]$red.memb))$lengths
  lpts <- cs - lens/2
  axis(1,at=lpts,labels=seq_len(length(color_list)),lwd.ticks=-0.1,cex.axis=1.25,padj=0.25,font=2)
  #dev.off()
}



#' Plot weighted adjacency matrix with links grouped by community
#' 
#' This function will generate the network link 'heatmap' for a weighted network
#' @param condor.object output of either \code{\link{condorCluster}} or 
#' \code{\link{condorModularityMax}}
#' @param main plot title
#' @param xlab x axis label
#' @param ylab y axis label
#' @return produces a \code{\link[graphics]{plot}} output.
#' @examples
#' data(small1976)
#' condor.object <- createCondorObject(small1976)
#' condor.object <- condorCluster(condor.object, project=FALSE)
#' condorPlotHeatmap(condor.object)
#' @import gplots
#' @export
#'  
condorPlotHeatmap = function(condor.object, main="", xlab="blues", ylab="reds"){
  bo <- condor.object
  # convert edge lists to adjacency matrices (n reds x m blues)
  adj = get.adjacency(bo$G, attr="weight", sparse=FALSE)
  # reorder reds according to community membership
  reds = as.character(bo$red.memb[order(bo$red.memb[,2]),1])
  adj = adj[reds,]
  # reorder blues according to community membership
  blues = as.character(bo$blue.memb[order(bo$blue.memb[,2]),1])
  adj = adj[,blues]
  rowsep = cumsum(as.vector(table(bo$red.memb[,2])))
  colsep = cumsum(as.vector(table(bo$blue.memb[,2])))
  labCol <- as.character(sort(bo$blue.memb[,2]))
  labCol[duplicated(labCol)] <- ""
  labRow <- as.character(sort(bo$red.memb[,2]))
  labRow[duplicated(labRow)] <- ""
  heatmap.2(adj, Rowv=FALSE, Colv=FALSE, dendrogram="none", keysize=1.25,
            col=colorpanel(10, "white", "black"), scale="none",
            key=TRUE, symkey=FALSE, density.info="none", trace="none",
            main=main, sepcolor ="#DDDDDD", colsep=colsep, rowsep=rowsep,
            sepwidth = c(0.025, 0.025), ylab=ylab, xlab=xlab, margins=c(3,3),
            labCol=labCol, labRow=labRow, offsetRow=0, offsetCol=0,
            breaks=sort(c(0.1,seq(0, max(adj),length.out=10))))
}


#'Calculate Qscore for all nodes 
#'
#'Qscore is designed to calculate the fraction of the modularity 
#'contributed by each node to its community's modularity 
#' @param condor.object output of \code{\link{condorCluster}} or 
#' \code{\link{condorModularityMax}}
#' @return condor.object list has \code{condor.object$qscores} added to it.
#' this includes two data.frames, \code{blue.qscore} and \code{red.qscore}
#' which have the qscore for each red and blue node.
#' 
#' @examples
#' r = c(1,1,1,2,2,2,3,3,3,4,4);
#' b = c(1,2,3,1,2,4,2,3,4,3,4);
#' reds <- c("Alice","Sue","Janine","Mary")
#' blues <- c("Bob","John","Ed","Hank")
#' elist <- data.frame(red=reds[r],blue=blues[b])
#' condor.object <- createCondorObject(elist)
#' condor.object <- condorCluster(condor.object)
#' condor.object <- condorQscore(condor.object)   
#' 
#' @import igraph
#' @import Matrix
#' @export
#' 
condorQscore = function(condor.object){
  
  if(is.null(condor.object$red.memb) | is.null(condor.object$blue.memb)){
    stop("Community Memberships missing. Run condorCluster or condorModularityMax first!")
  }
  bo <- condor.object
  bo$blue.memb <- bo$blue.memb[order(bo$blue.memb[,"blue.names"]),]
  bo$red.memb <- bo$red.memb[order(bo$red.memb[,"red.names"]),]
  bo$Qcoms <- bo$Qcoms[order(bo$Qcoms[,"community"]),]
  condor.object <- bo
  
  #Convert the edgelist to a sparseMatrix object
  esub <- condor.object$edges
  reds = as.integer(factor(esub[,1]))
  blues = as.integer(factor(esub[,2]))
  
  #extract weights, if any
  if(ncol(esub) > 2){
    weights <- esub[, 3]
  }
  else {
    weights <- rep(1, nrow(esub))
  }
  
  #Sparse matrix with the upper right block of the true Adjacency matrix. notices the dimension is reds x blues
  A = sparseMatrix(i=reds,j=blues,x=weights,dims=c(length(unique(reds)),length(unique(blues))),index1=TRUE);
  if(nrow(A) < ncol(A)){A <- t(A)}
  #if(max(blues) > max(reds)){blues <- reds;}
  p = nrow(A)
  q = ncol(A)
  N = p+q
  ki = rowSums(A)
  dj = colSums(A)
  m = sum(ki) # m = sum(dj) too
  R1 = condor.object$red.memb
  T1 = condor.object$blue.memb
  r1 = cbind(as.numeric(factor(R1[,1])),R1[,2])
  t1 = cbind(as.numeric(factor(T1[,1])),T1[,2]) 
  Rtrans = sparseMatrix(i=r1[,2],j=r1[,1],x=1,dims=c(max(r1[,2]),length(unique(r1[,1]))),index1=TRUE);    
  T2 = sparseMatrix(i=t1[,1],j=t1[,2],x=1,dims=c(max(t1[,1]),max(t1[,2])),index1=TRUE);
  Qcoms <- condor.object$Qcoms
  Qjk = vector(length=q)
  for(j in seq_len(max(t1[,1]))){
    if(j %% 1000 == 0){print(paste(j,t1[j,]))}
    Bj = A[,j] - (ki*dj[j])/m;
    Qjk[j] = ((Rtrans[t1[j,2],] %*% Bj)/(2*m))*(1/Qcoms[t1[j,2],1])
  }  
  Qik = vector(length=p)
  for(i in seq_len(max(r1[,1]))){
    if(i %% 1000 == 0){print(i)}
    Bi = A[i,] - (ki[i]*dj)/m;
    Qik[i] = ((Bi %*% T2[,r1[i,2]])/(2*m))*(1/Qcoms[r1[i,2],1])  
  }    
  condor.object$qscores = list(blue.qscore=data.frame(T1,Qjk),red.qscore=data.frame(R1,Qik))
  return(condor.object)
}


#' Create list amenable to analysis using \code{condor} package.
#' 
#' Converts an edge list into a \code{list} which is then an input for 
#' other functions in the \code{condor} package.
#' @param edgelist a data.frame with 'red' nodes in the first column and
#' 'blue' nodes in the second column, representing links from the node in
#' the first column to the node in the second column. There must be more
#' unique 'red' nodes than 'blue' nodes. Optionally, a third column may be
#' provided to create a weighted network.
#' @param return.gcc if TRUE, returns the giant connected component
#' @return G is an igraph graph object with a 'color' attribute
#' based on the colnames of edgelist. This can be accessed via
#' V(g)$color, which returns a vector indicating red/blue. Use V(g)$name
#' with V(g)$color to identify red/blue node names
#' @return edges corresponding to graph G. If return.gcc=TRUE, includes only
#' those edges in the giant connected component.
#' @return Qcoms output from \code{\link{condorCluster}} or 
#' \code{\link{condorModularityMax}}
#' @return modularity \code{NULL} output from \code{\link{condorCluster}} 
#' or \code{\link{condorModularityMax}}
#' @return red.memb \code{NULL} output from \code{\link{condorCluster}} 
#' or \code{\link{condorModularityMax}}
#' @return blue.memb \code{NULL} output from \code{\link{condorCluster}} 
#' or \code{\link{condorModularityMax}}
#' @return qscores \code{NULL} output from \code{\link{condorQscore}} 
#' @examples 
#' r = c(1,1,1,2,2,2,3,3,3,4,4);
#' b = c(1,2,3,1,2,4,2,3,4,3,4);
#' reds <- c("Alice","Sue","Janine","Mary")
#' blues <- c("Bob","John","Ed","Hank")
#' elist <- data.frame(red=reds[r],blue=blues[b])
#' condor.object <- createCondorObject(elist)
#' 
#' @import igraph
#' @export 
#' 
createCondorObject <- function(edgelist,return.gcc=TRUE){
  
  # make sure first to columns are of class character
  edgelist[, 1] <- as.character(edgelist[, 1])
  edgelist[, 2] <- as.character(edgelist[, 2])
  if(any(is.na(edgelist))) {
    stop("NA's detected. Remove these from edgelist")
  }
  if(any(edgelist[, 1] == "") | any(edgelist[, 2] == "")) {
    stop("Empty strings detected. Remove these from edgelist")
  }
  if(sum(edgelist[, 1] %in% edgelist[, 2]) > 0) {
    stop("edgelist contains one or more nodes that appear in both red and blue columns.
         Check to make sure network is truly bipartite and nodes of each type appear in the
         same column of 'edgelist'.")
  }
  
  g <- graph_from_data_frame(edgelist,directed=FALSE)
  blue.indx <- V(g)$name %in% unique(edgelist[, 2])
  V(g)$color <- "red"
  V(g)$color[blue.indx] <- "blue"
  
  if(ncol(edgelist) > 2) {
    message("Weights detected. Building condor object with weighted edges.")
    E(g)$weight <- as.numeric(edgelist[, 3])
    # omit 0 edges
    edgelist <- edgelist[edgelist[, 3] != 0, ]
  }
  
  if(!return.gcc){ g.out <- g}
  if(return.gcc){ gcc <- max.component(g); g.out <- gcc }
  blue.names <- V(g.out)$name[V(g.out)$name %in% unique(edgelist[, 2])]
  red.names <- V(g.out)$name[V(g.out)$name %in% unique(edgelist[, 1])]
  edges <- edgelist[edgelist[, 2] %in% blue.names,]
  
  return(list(G=g.out,edges=edges,Qcoms=NULL,modularity=NULL,red.memb=NULL,
              blue.memb=NULL,qscores=NULL))
  }


max.component = function(g){
  # return largest connected component of the iGraph graph object g
  g.clust = components(g);
  maxclust.id = which(g.clust$csize == max(g.clust$csize))[1];
  h = induced_subgraph(g, which(g.clust$membership == maxclust.id)); # 1-indexed here
  return(h);
}



#' Pollinator-plant interactions
#'
#' A dataset containing the number of interactions 34 plants and 13 pollinators. The variables are as follows:
#'
#' \itemize{
#'   \item pollinator. Species name of insect pollinator
#'   \item plant. Species name of plant
#'   \item interactions. Number of visitors caught on each plant species
#' }
#'
#' @docType data
#' @references \url{https://www.nceas.ucsb.edu/interactionweb/html/small_1976.html}
#' @keywords datasets
#' @name small1976
#' @usage data(small1976)
#' @format A data frame with 442 rows and 3 variables
NULL

globalVariables(c('SNP','blue.memb',"G", "bindx", "blue.com.size","gene",  "red.com.size","red.memb", "rindx"))
