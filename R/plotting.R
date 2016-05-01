#' Comparison of Z scores between two PANDA runs
#'
#' Given two PANDA objects with the same network structure, plot the Z-score comparison.
#' The two PANDA objects should only differ in the gene expression used for the network 
#' constructions or other parameters.
#'
#' @param x PANDA object - output of the \code{panda} function.
#' @param y PANDA object - second PANDA object.
#' @param hex TRUE/FALSE - If TRUE, bin data points to avoid over plotting.
#' @param bins Number of bins to use for plotting.
#' @param addLine  TRUE/FALSE - to add y=x line.
#' @param rank TRUE/FALSE - If TRUE, plot rank of edge weights rather than weight values.
#' @return ggplot comparing the Z-scores between the two networks.
#' @import hexbin
#' @importFrom reshape melt.array
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 ylab
#' @importFrom ggplot2 aes_string
#' @importFrom ggplot2 geom_abline
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 stat_binhex
#' @export
#' @examples
#' 
#' data(pandaResult)
#' data(pandaToyData)
#' pandaRes <- pandaRes2 <- pandaResult
#' plotZ(pandaRes, pandaRes2)
#'
#' \donttest{
#' panda.res1 <- with(pandaToyData, panda(motif, expression, ppi, hamming=1))
#' panda.res2 <- with(pandaToyData, panda(motif, expression + 
#'                    rnorm(prod(dim(expression)),sd=5), ppi, hamming=1))
#' plotZ(panda.res1, panda.res2,addLine=FALSE)
#' }
#'
plotZ <- function(x,y,hex=TRUE,bins=200,addLine=TRUE,rank=FALSE){
  x = melt.array(slot(x,"regNet"))
  y = melt.array(slot(y,"regNet"))
  d <- merge(x,y,by=colnames(x)[1:2])
  if (rank) {
    d$value.x <- rank(-d$value.x)
    d$value.y <- rank(-d$value.y)
  }
  p <- ggplot(d, aes_string(x="value.x", y="value.y")) +
    xlab(sprintf("Network X %s", ifelse(rank, "rank", "Z-score"))) +
    ylab(sprintf("Network Y %s", ifelse(rank, "rank", "Z-score")))
  if (addLine) {
    p <- p + geom_abline(slope=1, col="red")
  }
  if (hex) {
    p <- p + stat_binhex(bins=bins)
  } else {
    p <- p + geom_point(shape=1)
  }
  p
}
#' Check motif
#' 
#' This function adds random false positive edges to the regulatory prior
#' and will check if they become pruned.
#'
#' @param x Model regulatory network.
#' @param motif Motif used to construct the model regulatory network.
#' @param expr Expression matrix used to construct model network.
#' @param ppi PPI used to construct model regulatory network.
#' @param mode a character string - either "augment" to add random edges or "remove" to remove random edges.
#' @param prop numeric specifying number of edges to augment or remove from regulatory prior, as a proportion of the number
#'         of edges in the regulatory prior.
#' @param seed Random seed.
#' @param ... Options for the panda function.
#' @return ggplot heatmap list of indices of net corresponding to each TF
#' @importFrom reshape melt.array
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 ylab
#' @importFrom ggplot2 aes_string
#' @importFrom ggplot2 scale_color_discrete
#' @export 
#' @examples
#' 
#' data(pandaToyData)
#' data(pandaResult)
#' regnet = slot(pandaResult,"regNet")
#' with(pandaToyData, testMotif(regnet, motif, mode="augment", expression, ppi, hamming=1))
#'
testMotif <- function(x,motif,expr,ppi,mode=c("augment","remove"),prop=0.05,seed=1,...) {
  mode <- match.arg(mode)
  net <- melt.array(x)
  colnames(motif) <- colnames(net) <- c("TF", "Gene", "Score")
  set.seed(seed)
  false.idx <- c()
  if (mode=="augment") {
    tf.idx = list()
    motif.table <- table(motif[,1])
    for (tf in names(motif.table)) {
      tf.idx[[tf]] <- which(net$TF==tf)
    }
    # sample net with same distribution of TFs as the motif data
    for (i in 1:length(motif.table)) {
      tf <- names(motif.table)[i]
      false.idx <- c(false.idx, sample(tf.idx[[tf]], motif.table[i]))
    }
    true.idx <- which(paste0(net$TF, net$Gene)%in%paste0(motif$TF, motif$Gene))
    false.idx <- setdiff(false.idx, true.idx)
    false.idx <- sample(false.idx, nrow(motif) * prop)
    motif.mod <- cbind(net[false.idx, c("TF", "Gene")], Score=1)
    motif.mod <- rbind(motif, motif.mod)
  }
  if (mode=="remove") {
    rm.idx <- sample(nrow(motif), nrow(motif) * prop)
    motif.mod <- motif[-rm.idx,]
    # remove randomly selected TFs from PPI
    ppi <- ppi[which(ppi[,1]%in%motif.mod$TF & ppi[,2]%in%motif.mod$TF),]
    false.idx <- which(paste0(net$TF, net$Gene)%in%paste0(motif[rm.idx,]$TF, motif[rm.idx,]$Gene))
  }
  panda.mod <- panda(motif.mod, expr, ppi, ...)
  reg.mod <- slot(panda.mod,"regNet")
  net.mod <- melt.array(reg.mod)
  colnames(net.mod) <- c("TF", "Gene", "Score")
  
  net$x <- 1:nrow(net) %in% false.idx
  d <- merge(net, net.mod, by=c("TF", "Gene"), all.x=TRUE, sort=FALSE)
  d <- d[order(d$x),]
  p <- ggplot(d, aes_string(x="Score.x", y="Score.y", color="x==TRUE")) + geom_point(shape=1) +
    xlab("Model network weights") + ylab("Modified network weights") +
    scale_color_discrete(name="", labels=c("Model edges", "False edges"))
  p
}
#' Plot Z by TF out-degree quantiles
#'
#' Generates a Z-score scatterplot for edges according to the TF outdegree in prior.
#' The two PANDA objects should only differ in the gene expression used for the network 
#' constructions or other parameters.
#'
#' @param x PANDA object - output of the \code{panda} function.
#' @param y PANDA object - second PANDA object.
#' @param motif Motif used to construct the networks.
#' @param hasPrior TRUE/FALSE, If TRUE plots the edges that are given a weight > 0 in the motif,
#'              else plot those given a weight of 0
#' @param cuts either a numeric vector of two or more unique cut points or a
#'        single number (greater than or equal to 2) giving the number
#'        of intervals into which 'x' is to be cut.
#' @param cols Number of columns in layout.
#' @export
#' @return ggplot heatmapfor each TF, get outdegree in regulatory prior
#' @importFrom reshape melt.array
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 ylab
#' @importFrom ggplot2 aes_string
#' @importFrom ggplot2 xlim
#' @importFrom ggplot2 ylim
#' @importFrom ggplot2 ggtitle
#' @importFrom ggplot2 scale_fill_gradient
#' @importFrom ggplot2 stat_binhex
#' @importFrom ggplot2 "%+%"
#' @examples
#'
#' data(pandaResult)
#' data(pandaToyData)
#' plotZbyTF(pandaResult,pandaResult, pandaToyData$motif, hasPrior=TRUE)
#' plotZbyTF(pandaResult,pandaResult, pandaToyData$motif, cuts=2)
#' 
plotZbyTF <- function(x, y, motif,hasPrior=TRUE,cuts=1,cols=2){
  x = slot(x,"regNet")
  y = slot(y,"regNet")
  motif.table <- aggregate(motif[,3], list(TF=motif[,1]), sum)
  motif.counts <- data.frame(TF=motif.table$TF, outdegree=motif.table$x)
  motif.counts <- motif.counts[order(motif.counts$outdegree),]
  if(cuts == 1){
    quant = factor(rep(1,nrow(x)))
  } else {
    quant <- cut(motif.counts$outdegree,breaks=cuts,include.lowest=TRUE)
  }
  tfs.q <- list()
  for (i in 1:cuts){
    tfs.q[[i]] <- motif.counts[quant==levels(quant)[i], ]$TF
  }
  x = melt.array(x)
  y = melt.array(y)
  d <- merge(x,y,by=colnames(x)[1:2])
  colnames(d) <- c("TF", "Gene", "Z.x", "Z.y")
  d$Motif <- paste0(d$TF, d$Gene) %in% paste0(motif[,1], motif[,2])
  pList <- list()
  
  if (!hasPrior) {
    d <- d[d$Motif==0, ]
  } else {
    d <- d[d$Motif>0, ]
  }
  lim <- c(min(d$Z.x, d$Z.y), max(d$Z.x, d$Z.y))
  for (i in 1:cuts){
    p <- ggplot(d[d$TF%in%tfs.q[[i]], ], aes_string(x="Z.x", y="Z.y")) +
      geom_abline(intercept=0, slope=1, colour="red") + ggtitle(sprintf("Q %s",i)) +
      xlim(lim) + ylim(lim)
    if (!hasPrior) {
      p <- p + stat_binhex(bins=200) +
        scale_fill_gradient(name = "count", trans = "log", breaks = 10^(0:10))
    } else {
      p <- p + geom_point(shape=1)
    }
    pList[[i]] <- p
  }
  multiplot(plotlist=pList, cols=cols)
}

#' Multiple plots
#'
#' Multiple plot function as described in: 
#' http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/.
#' If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
#' then plot 1 will go in the upper left, 2 will go in the upper right, and
#' 3 will go all the way across the bottom.
#'
#' @param ... ggplot objects can be passed in or to plotlist (as a list of ggplot objects).
#' @param plotlist NULL - if the plotlist is null, the function will figure out the panel dimensions.
#' @param cols Number of columns in layout.
#' @param layout A matrix specifying the layout. If present, 'cols' is ignored.
#' @export
#' @importFrom grid grid.newpage
#' @importFrom grid grid.layout
#' @importFrom grid viewport
#' @importFrom grid pushViewport
#'
multiplot <- function(...,plotlist=NULL,cols=1,layout=NULL) {
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  numPlots = length(plots)
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  if (numPlots==1) {
    print(plots[[1]])
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
#' Community detection plot
#' 
#' This function performs community detection on an undirected PANDA network.
#' The function optionally returns the graph and community.
#'
#' @param x Toy PANDA output represented as a TF, Gene, and Score.
#' @param scaleEdge Visualization parameter for the edges.
#' @param verbose TRUE/FALSE - Report community structure.
#' @param ... Options for the plot function.
#' @return Optionally return a list with the graph and community.
#' @importFrom igraph graph.data.frame
#' @importFrom igraph fastgreedy.community
#' @importFrom igraph membership
#' @importFrom igraph E<-
#' @importFrom igraph E
#' @importFrom igraph V
#' @importFrom igraph V<-
#' @importFrom igraph layout.fruchterman.reingold
#' @importFrom igraph groups
#' @export 
#' @examples
#' 
#' # start with some toy PANDA output
#' mat <- cbind(rep(1:5, each=10), rep(seq(11,20),5), sample(100, 50)/100)
#' x =plotCommunityDetection(mat)
#' str(x)
#'
#' #example of very different edges
#' set.seed(1)
#' subst <- sample(50,10)
#' mat[subst, 3] <- subst
#' plotCommunityDetection(mat,scaleEdge=0.5)
#'
plotCommunityDetection<-function(x,scaleEdge = 5,verbose=TRUE,...){
  # PANDA networks are directed, but for now let's treat them as undirected
  g <- graph.data.frame(x, directed=F) 
  # or 2^x[,3] (we might want to transform the PANDA z-scores)
  E(g)$weight <- x[,3]
  com <- fastgreedy.community(g,merges=TRUE,modularity=TRUE,membership=TRUE,weights=E(g)$weight)

  if(verbose){
    # make igraph object and perform community structure analysis
    # print max modularity
    cat("Modularity (best split):", max(com$modularity),"\n")
    # print communities to screen
    for (m in 1:length(groups(com))){
      cat("Community", m, ":", names(membership(com)[which(membership(com)==m)]), "\n")
    }
  }

  # color code nodes belonging to the same community
  V(g)$color <- com$membership
  # plot network
  plot(g, vertex.label.cex=1, vertex.size=15, vertex.label.color = "white", vertex.label.font=2,
    edge.width=E(g)$weight*scaleEdge, layout=layout.fruchterman.reingold, vertex.color=V(g)$color,...)
  title(paste("Modularity", round(max(com$modularity),2)))
  invisible(list(graph=g,community=com))
}
