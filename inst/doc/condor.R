### R code from vignette source 'condor.Rnw'

###################################################
### code chunk number 1: condor.Rnw:31-33
###################################################
library(condor)
library(igraph)


###################################################
### code chunk number 2: condor.Rnw:36-41
###################################################
r = c(1,1,1,2,2,2,3,3,3,4,4);
b = c(1,2,3,1,2,4,2,3,4,3,4);
reds <- c("Alice","Sue","Janine","Mary")
blues <- c("Bob","John","Ed","Hank")
elist <- data.frame(red=reds[r], blue=blues[b])


###################################################
### code chunk number 3: condor.Rnw:44-45
###################################################
condor.object <- create.condor.object(elist)


###################################################
### code chunk number 4: condor.Rnw:48-49
###################################################
names(condor.object)


###################################################
### code chunk number 5: condor.Rnw:52-55
###################################################
condor.object <- condor.cluster(condor.object)
print(condor.object$red.memb)
print(condor.object$blue.memb)


###################################################
### code chunk number 6: condor.Rnw:59-62
###################################################
gtoy = graph.edgelist(as.matrix(elist),directed=FALSE)
set.graph.attribute(gtoy, "layout", layout.kamada.kawai(gtoy))
V(gtoy)[c(reds,blues)]$color <- c(rep("red",4),rep("blue",4))


###################################################
### code chunk number 7: condor.Rnw:65-66
###################################################
plot(gtoy,vertex.label.dist=2)


###################################################
### code chunk number 8: condor.Rnw:70-71
###################################################
condor.object <- condor.qscore(condor.object)


###################################################
### code chunk number 9: condor.Rnw:75-78
###################################################
q_women <- condor.object$qscores$red.qscore
core_stats <- condor.core.enrich(test_nodes=c("Alice","Mary"),
                                 q=q_women,perm=TRUE,plot.hist=TRUE)


###################################################
### code chunk number 10: condor.Rnw:82-85
###################################################
data(small1976)
condor.object <- create.condor.object(small1976)
condor.object <- condor.cluster(condor.object, project=F)


###################################################
### code chunk number 11: condor.Rnw:88-89
###################################################
condor.plot.heatmap(condor.object)


