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
#' @examples
#' r = c(1,1,1,2,2,2,3,3,3,4,4);
#' b = c(1,2,3,1,2,4,2,3,4,3,4);
#' reds <- c("Alice","Sue","Janine","Mary")
#' blues <- c("Bob","John","Ed","Hank")
#' elist <- data.frame(red=reds[r],blue=blues[b])
#' condor.object <- create.condor.object(elist)
#' condor.object <- condor.cluster(condor.object)
#' condor.object <- condor.qscore(condor.object) 
#' q_in <- condor.object$qscores$red.qscore
#' out <- condor.core.enrich(c("Alice","Mary"),q=q_in,perm=TRUE,plot.hist=TRUE)
#' @export
#' 
condor.core.enrich = function(test_nodes,q,perm=FALSE,plot.hist=FALSE,nsamp=1000){
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
    for(i in 1:nsamp){
        #randomly choose n values from all, assign the rest to B_rand
        aind = sample(1:length(all),n)
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
    for(i in 1:nsamp){
        #randomly choose n values from all, assign the rest to B_rand
        aind = sample(1:length(all),n)
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
        rect(min(h1$breaks),0,max(c(h1$breaks,w_true))+0.1,1.03*max(h1$counts),lwd=2)
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
