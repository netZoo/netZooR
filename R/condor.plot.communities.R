#' Plot adjacency matrix with links grouped and colored by community
#' 
#' This function will generate the network link 'heatmap' with colored dots
#' representing within-community links and black dots between-community 
#' links
#' @param condor.object output of either \code{\link{condor.cluster}} or 
#' \code{\link{condor.modularity.max}}
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
#' condor.object <- create.condor.object(elist)
#' condor.object <- condor.cluster(condor.object)
#' condor.plot.communities(condor.object,
#' color_list=c("darkgreen","darkorange"),point.size=2,
#' xlab="Women",ylab="Men")
#' @import data.table
#' @export
#'  
condor.plot.communities = function(condor.object,color_list,point.size=0.01,
                                   xlab="SNP",ylab="Gene"){
  
    dt0 <- data.table(condor.object$edges)
    setnames(dt0,c("red","blue"),c("SNP","gene"))
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
    
    #setkeyv(eqtl1,c("SNP","blue.memb","gene","red.memb"))
    if(nlevels(eqtl_block$SNP) != length(unique(eqtl_block$SNP))){
        print("warning: empty levels in SNP column. This may cause silent issues with plotting.")}
    #select all links that connect nodes in the same community
    setkey(eqtl_block,"blue.memb","red.memb")
    #make new index for each node that will correspond to it's row/col number
    red_tmp <- data.table(rindx=1:nlevels(eqtl_block$SNP),SNP=unique(eqtl_block$SNP))
    red_indx <- merge(red_tmp,unique(eqtl_block,by="SNP")[,c("SNP","red.memb"),with=FALSE],by="SNP")
    red_indx[,red.com.size:=length(unique(SNP)),by=red.memb]
    red_indx[red.com.size > 1,rindx:=sample(x=rindx),by=red.memb][,red.memb:=NULL,]
    setkey(red_indx,"SNP")
    blue_tmp <- data.table(bindx=1:nlevels(eqtl_block$gene),gene=unique(eqtl_block$gene))
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
    axis(1,at=lpts,labels=1:length(color_list),lwd.ticks=-0.1,cex.axis=1.25,padj=0.25,font=2)
    #dev.off()
}
