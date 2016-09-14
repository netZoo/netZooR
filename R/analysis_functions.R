#' Bi-partite network analysis tools
#'
#' This function analyzes a bi-partite network, such as a Transcription factor 
#' to gene network derived from the PANDA algorithm.
#'
#' @param net1 starting network, a genes by transcription factors data.frame with scores 
#' for the existence of edges between
#' @param net2 final network, a genes by transcription factors data.frame with scores 
#' for the existence of edges between
#' @param by.tfs logical indicating a transcription factor based transformation.    If 
#' false, gives gene by gene transformation matrix
#' @param remove.diagonal logical for returning a result containing 0s across the diagonal
#' @param method character specifying which algorithm to use, default='kabsch'
#' @return matrix object corresponding to transition matrix
#' @keywords keywords
#' @import MASS
#' @importFrom testthat test_that
#' @export
#' @examples
#' data(yeast)
#' cc.net <- monsterNI(yeast$motif,yeast$exp.cc)
#' sr.net <- monsterNI(yeast$motif,yeast$exp.sr)
#' transformation.matrix(cc.net, sr.net)
transformation.matrix <- function(network.1, network.2, by.tfs=TRUE, standardize=FALSE, 
                                remove.diagonal=TRUE, method="ols"){
    require(MASS)
    if(is.list(network.1)&&is.list(network.2)){
        if(by.tfs){
            net1 <- t(network.1$reg.net)
            net2 <- t(network.2$reg.net)
        } else {
            net1 <- network.1$reg.net
            net2 <- network.2$reg.net
        }
    } else if(is.matrix(network.1)&&is.matrix(network.2)){
        if(by.tfs){
            net1 <- t(network.1)
            net2 <- t(network.2)
        } else {
            net1 <- network.1
            net2 <- network.2
        }
    } else {
        stop("Networks must be lists or matrices")
    }
    #gene.trans.matrix <- svd(net2)$v %*% diag(1/svd(net2)$d) %*% t(svd(net2)$u) %*% net1
    if (method == "kabsch"){
        tf.trans.matrix <- kabsch(net1,net2)
    }
    if (method == "old"){
        svd.net2 <- svd(net2)
        tf.trans.matrix <- svd.net2$v %*% diag(1/svd.net2$d) %*% t(svd.net2$u) %*% net1
    }
    if (method == "ols"){
        # 9/14/15
        # re-rewrote 'same column priority' feature
        net2.star <- sapply(1:ncol(net1), function(i,x,y){
            lm(y[,i]~x[,i])$resid
        }, net1, net2)
        tf.trans.matrix <- ginv(t(net1)%*%net1)%*%t(net1)%*%net2.star
        colnames(tf.trans.matrix) <- colnames(net1)
        rownames(tf.trans.matrix) <- colnames(net1)
        print("Using OLS method")

    }
    if (method == "L1"){
            net2.star <- sapply(1:ncol(net1), function(i,x,y){
                    lm(y[,i]~x[,i])$resid
            }, net1, net2)
            tf.trans.matrix <- sapply(1:ncol(net1), function(i){
                    z <- optL1(net2.star[,i], net1, fold=5, minlambda1=1, 
                            maxlambda1=2, model="linear", standardize=TRUE)
                    coefficients(z$fullfit, "penalized")
            })
            colnames(tf.trans.matrix) <- rownames(tf.trans.matrix)
            print("Using L1 method")

    }
    if (standardize){
        tf.trans.matrix <- apply(tf.trans.matrix, 1, function(x){
            #     x.zero <- (x-mean(x))
            x/sum(abs(x))
        })
    }

    if (remove.diagonal){
        diag(tf.trans.matrix) <- 0
    }
    # Add column labels
    colnames(tf.trans.matrix) <- rownames(tf.trans.matrix)
    tf.trans.matrix
}
plottm <- function(melt.tm,title="Transition Matrix"){
    ggplot(melt.tm, aes(y=Var1,x=Var2)) +
        ggtitle(title) +
        geom_tile(aes(fill=value)) +
        scale_fill_gradient(low="blue", high="yellow") +
        xlab("") +
        ylab("") +
        theme(axis.text.x = element_text(angle = 90,hjust=1,vjust=0.5))
}
kabsch <- function(P,Q){

    P <- apply(P,2,function(x){
        x - mean(x)
    })
    Q <- apply(Q,2,function(x){
        x - mean(x)
    })
    covmat <- cov(P,Q)
    P.bar <- colMeans(P)
    Q.bar <- colMeans(Q)
    num.TFs <- ncol(P)        #n
    num.genes <- nrow(P)    #m

    #     covmat <- (t(P)%*%Q - P.bar%*%t(Q.bar)*(num.genes))

    svd.res <- svd(covmat-num.TFs*Q.bar%*%t(P.bar))

    # Note the scalar multiplier in the middle.
    # NOT A MISTAKE!
    c.k <- colSums(P %*% svd.res$v * Q %*% svd.res$u) - 
        num.genes*(P.bar%*%svd.res$v)*(Q.bar%*%svd.res$u)

    E <- diag(c(sign(c.k)))

    W <- svd.res$v %*% E %*% t(svd.res$u)
    rownames(W) <- colnames(P)
    colnames(W) <- colnames(P)
    W
}

#' Sum of squared off-diagonal mass
#'
#' This function calculates the off-diagonal sum of squared mass for a transition matrix
#'
#' @param tm a transition matrix for two bipartite networks
#' @keywords keywords
#' @export
#' @return vector containined the sum of squared off diagonal mass, dTFI
#' @examples
#' mat <- matrix(rnorm(100),ncol=10)
#' ssodm(mat)
#' 
ssodm <-    function(tm){
    diag(tm)<-0
    apply(tm,1,function(x){t(x)%*%x})
}

#' Transformation matrix plot
#'
#' This function plots a hierachically clustered heatmap and 
#' corresponding dendrogram of a transaction matrix
#'
#' @param x monster Object
#' @param method distance metric for hierarchical clustering.    
#' Default is "Pearson correlation"
#' @keywords keywords
#' @export
#' @import ggplot2
#' @import grid
#' @return ggplot2 object for transition matrix heatmap
#' @examples
#' data(yeast)
#' cc.net <- monsterNI(yeast$motif,yeast$exp.cc)
#' sr.net <- monsterNI(yeast$motif,yeast$exp.sr)
#' transformation.matrix(cc.net, sr.net)
hcl.heatmap.plot <- function(x, method="pearson"){
    assert_that(class(x)=="monster")
    x <- x@tm
    if(method=="pearson"){
        dist.func <- function(y) as.dist(cor(y))
    } else {
        dist.func <- dist
    }
    # x <- as.matrix(scale(mtcars))
    x <- scale(x)
    dd.col <- as.dendrogram(hclust(dist.func(x)))
    col.ord <- order.dendrogram(dd.col)

    dd.row <- as.dendrogram(hclust(dist.func(t(x))))
    row.ord <- order.dendrogram(dd.row)

    xx <- x[col.ord, row.ord]
    xx_names <- attr(xx, "dimnames")
    df <- as.data.frame(xx)
    colnames(df) <- xx_names[[2]]
    df$Var1 <- xx_names[[1]]
    df$Var1 <- with(df, factor(Var1, levels=Var1, ordered=TRUE))
    mdf <- reshape2::melt(df)


    ddata_x <- dendro_data(dd.row)
    ddata_y <- dendro_data(dd.col)

    ### Set up a blank theme
    theme_none <- theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title.x = element_text(colour=NA),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.line = element_blank()
        #axis.ticks.length = element_blank()
    )
    ### Set up a blank theme
    theme_heatmap <- theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title.x = element_text(colour=NA),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.line = element_blank()
        #axis.ticks.length = element_blank()
    )
    ### Create plot components ###
    # Heatmap
    p1 <- ggplot(mdf, aes(x=variable, y=Var1)) +
        geom_tile(aes(fill=value)) + 
        scale_fill_gradient2() + 
        theme(axis.text.x = element_text(angle = 90, hjust = 1))

    # Dendrogram 1
    p2 <- ggplot(segment(ddata_x)) +
        geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) +
        theme_none + theme(axis.title.x=element_blank())

    # Dendrogram 2
    p3 <- ggplot(segment(ddata_y)) +
        geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) +
        coord_flip() + theme_none

    ### Draw graphic ###

    grid.newpage()
    print(p1, vp=viewport(0.80, 0.8, x=0.400, y=0.40))
    print(p2, vp=viewport(0.73, 0.2, x=0.395, y=0.90))
    print(p3, vp=viewport(0.20, 0.8, x=0.910, y=0.43))
}

#' Principal Components plot of transformation matrix
#'
#' This function plots the first two principal components for a 
#' transaction matrix
#'
#' @param tm a transition matrix for a bipartite network
#' @param title The title of the plot
#' @param clusters A vector indicating the colors to be plotted for each node
#' @param alpha A vector indicating the level of transparency to be plotted 
#' for each node
#' @return ggplot2 object for transition matrix PCA
#' @keywords keywords
#' @import ggdendro
#' @export
#' @examples
#' data(yeast)
#' monsterRes <- monster(yeast$exp.ko,c(rep(1,42),rep(0,49),rep(NA,15)),
#'     yeast$motif, nullPerms=10, numMaxCores=4)
#' # Color the nodes according to cluster membership
#' clusters <- kmeans(slot(monsterRes, 'tm'),3)$cluster 
#' transitionPCAPlot(monsterRes, 
#' title="PCA Plot of Transition - Cell Cycle vs Stress Response", 
#' clusters=clusters)
transitionPCAPlot <-    function(monsterObj, 
                                title="PCA Plot of Transition", 
                                clusters=1, alpha=1){
    require(ggplot2)
    tm.pca <- princomp(monsterObj@tm)
    odsm <- apply(monsterObj@tm,2,function(x){t(x)%*%x})
    odsm.scaled <- 2*(odsm-mean(odsm))/sd(odsm)+4
    scores.pca <- as.data.frame(tm.pca$scores)
    scores.pca <- cbind(scores.pca,'node.names'=rownames(scores.pca))
    ggplot(data = scores.pca, aes(x = Comp.1, y = Comp.2, label = node.names)) +
        geom_hline(yintercept = 0, colour = "gray65") +
        geom_vline(xintercept = 0, colour = "gray65") +
        geom_text(size = odsm.scaled, alpha=alpha, color=clusters) +
        ggtitle(title)
}

#' This function uses igraph to plot the transition matrix as a network
#'
#' @param monsterObj Monster Object
#' @param numEdges The number of edges to display
#' @param numTopTFs The number of TFs to display, ranked by largest dTFI
#' @return igraph object for transition matrix
#' @keywords keywords
#' @import igraph
#' @export
#' @examples
#' data(yeast)
#' monsterRes <- monster(yeast$exp.ko,
#'     c(rep(1,42),rep(0,49),rep(NA,15)),
#'     yeast$motif, nullPerms=10, numMaxCores=4)
#' transitionNetworkPlot(monsterRes)
#' 
transitionNetworkPlot <- function(monsterObj, numEdges=100, numTopTFs=10){
    require(reshape2)
    require(igraph)
    ## Calculate p-values for off-diagonals
    transitionSigmas <- function(tm.observed, tm.null){
        tm.null.mean <- apply(simplify2array(tm.null), 1:2, mean)
        tm.null.sd <- apply(simplify2array(tm.null), 1:2, sd)
        sigmas <- (tm.observed - tm.null.mean)/tm.null.sd
    }
    
    tm.sigmas <- transitionSigmas(monsterObj@tm, monsterObj@nullTM)
    diag(tm.sigmas) <- 0
    tm.sigmas.melt <- melt(tm.sigmas)
    
    adjMat <- monsterObj@tm
    diag(adjMat) <- 0
    adjMat.melt <- melt(adjMat)
    
    adj.combined <- merge(tm.sigmas.melt, adjMat.melt, by=c("Var1","Var2"))
    
    # adj.combined[,1] <- mappings[match(adj.combined[,1], mappings[,1]),2]
    # adj.combined[,2] <- mappings[match(adj.combined[,2], mappings[,1]),2]
    
    dTFI_pVals_All <- 1-2*abs(.5-calculate.tm.p.values(monsterObj, 
                                                method="z-score"))
    topTFsIncluded <- names(sort(dTFI_pVals_All)[1:numTopTFs])
    topTFIndices <- 2>(is.na(match(adj.combined[,1],topTFsIncluded)) + 
        is.na(match(adj.combined[,2],topTFsIncluded)))
    adj.combined <- adj.combined[topTFIndices,]
    adj.combined <- adj.combined[
        abs(adj.combined[,4])>=sort(abs(adj.combined[,4]),decreasing=TRUE)[numEdges],]
    tfNet <- graph.data.frame(adj.combined, directed=TRUE)
    vSize <- -log(dTFI_pVals_All)
    vSize[vSize<0] <- 0
    vSize[vSize>3] <- 3
    
    V(tfNet)$size <- vSize[V(tfNet)$name]*5
    V(tfNet)$color <- "yellow"
    E(tfNet)$width <- (abs(E(tfNet)$value.x))*15/max(abs(E(tfNet)$value.x))
    E(tfNet)$color<-ifelse(E(tfNet)$value.x>0, "blue", "red")
    
    plot.igraph(tfNet, edge.arrow.size=2, vertex.label.cex= 1.5, vertex.label.color= "black",main="")
}

#' This function plots the Off diagonal mass of an 
#' observed Transition Matrix compared to a set of null TMs
#'
#' @param monsterObj Monster Object
#' @param logical indicating whether to reorder transcription
#' factors according to their statistical significance and to 
#' rescale the values observed to be standardized by the null
#' distribution 
#' @param plot.title String specifying the plot title
#' @param highlight.tfs vector specifying a set of transcription 
#' factors to highlight in the plot
#' @return ggplot2 object for transition matrix comparing observed 
#' distribution to that estimated under the null 
#' @keywords keywords
#' @export
#' @examples
#' 1+1
dTFIPlot <- function(monsterObj, rescale=FALSE, plot.title=NA, highlight.tfs=NA){
    require(ggplot2)
    if(is.na(plot.title)){
        plot.title <- "Differential TF Involvement"
    }
    num.iterations <- length(monsterObj@nullTM)
    # Calculate the off-diagonal squared mass for each transition matrix
    null.SSODM <- lapply(monsterObj@nullTM,function(x){
        apply(x,2,function(y){t(y)%*%y})
    })
    null.ssodm.matrix <- matrix(unlist(null.SSODM),ncol=num.iterations)
    null.ssodm.matrix <- t(apply(null.ssodm.matrix,1,sort))

    ssodm <- apply(monsterObj@tm,2,function(x){t(x)%*%x})

    # Get p-value (rank of observed within null ssodm)
    #     p.values <- sapply(1:length(ssodm),function(i){
    #         1-findInterval(ssodm[i], null.ssodm.matrix[i,])/num.iterations
    #     })
    p.values <- 1-pnorm(sapply(1:length(ssodm),function(i){
        (ssodm[i]-mean(null.ssodm.matrix[i,]))/sd(null.ssodm.matrix[i,])
    }))
    t.values <- sapply(1:length(ssodm),function(i){
            (ssodm[i]-mean(null.ssodm.matrix[i,]))/sd(null.ssodm.matrix[i,])
    })

    # Process the data for ggplot2
    combined.mat <- cbind(null.ssodm.matrix, ssodm)
    colnames(combined.mat) <- c(rep('Null',num.iterations),"Observed")


    if (rescale){
        combined.mat <- t(apply(combined.mat,1,function(x){
            (x-mean(x[-(num.iterations+1)]))/sd(x[-(num.iterations+1)])
        }))
        x.axis.order <- rownames(monsterObj@nullTM[[1]])[order(-t.values)]
        x.axis.size    <- 10 # pmin(15,7-log(p.values[order(p.values)]))
    } else {
        x.axis.order <- rownames(monsterObj@nullTM[[1]])
        x.axis.size    <- pmin(15,7-log(p.values))
    }

    null.SSODM.melt <- reshape2::melt(combined.mat)[,-1][,c(2,1)]
    null.SSODM.melt$TF<-rep(rownames(monsterObj@nullTM[[1]]),num.iterations+1)

    ## Plot the data
    ggplot(null.SSODM.melt, aes(x=TF, y=value))+
        geom_point(aes(color=factor(Var2), alpha = .5*as.numeric(factor(Var2))), size=2) +
        scale_color_manual(values = c("blue", "red")) +
        scale_alpha(guide = "none") +
        scale_x_discrete(limits = x.axis.order ) +
        theme_classic() +
        theme(legend.title=element_blank(),
            axis.text.x = element_text(colour = 1+x.axis.order%in%highlight.tfs, 
            angle = 90, hjust = 1, 
            size=x.axis.size,face="bold")) +
        ylab("dTFI") +
        ggtitle(plot.title)

}

#' Calculate p-values for a tranformation matrix
#'
#' This function calculates the significance of an observed
#' transition matrix given a set of null transition matrices
#'
#' @param monsterObj Monster Object
#' @param method one of 'z-score' or 'non-parametric'
#' @return vector of p-values for each transcription factor
#' @keywords keywords
#' @export
#' @examples
#' 1+1
calculate.tm.p.values <- function(monsterObj, method="z-score"){
    num.iterations <- length(monsterObj@nullTM)
    # Calculate the off-diagonal squared mass for each transition matrix
    null.SSODM <- lapply(monsterObj@nullTM,function(x){
        apply(x,1,function(y){t(y)%*%y})
    })
    null.ssodm.matrix <- matrix(unlist(null.SSODM),ncol=num.iterations)
    null.ssodm.matrix <- t(apply(null.ssodm.matrix,1,sort))

    ssodm <- apply(monsterObj@tm,1,function(x){t(x)%*%x})

    # Get p-value (rank of observed within null ssodm)
    if(method=="non-parametric"){
        p.values <- sapply(1:length(ssodm),function(i){
            1-findInterval(ssodm[i], null.ssodm.matrix[i,])/num.iterations
        })
    } else if (method=="z-score"){
        p.values <- pnorm(sapply(1:length(ssodm),function(i){
            (ssodm[i]-mean(null.ssodm.matrix[i,]))/sd(null.ssodm.matrix[i,])
        }))
    } else {
        print('Undefined method')
    }
    p.values
}


