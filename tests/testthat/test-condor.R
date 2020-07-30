context("test CONDOR result")

test_that("CONDOR functions work", {
  
  load("./testDataset.RData")
  
  r = c(1,1,1,2,2,2,3,3,3,4,4);
  b = c(1,2,3,1,2,4,2,3,4,3,4);
  reds <- c("Alice","Sue","Janine","Mary")
  blues <- c("Bob","John","Ed","Hank")
  elist <- data.frame(red=reds[r], blue=blues[b])
  
  condor.object <- create.condor.object(elist)
  # check attribute names
  expect_equal(names(condor.object),c("G","edges","Qcoms","modularity","red.memb","blue.memb","qscores"))
  
  condor.object <- condor.cluster(condor.object,project = F)
  
  # check project=T when clustering
  
  expect_error(condor.cluster(condor.object,project = T),NA)
  
  # check modularity
  expect_equal(condor.object$modularity,as.numeric(c("0.231404958677686","0.231404958677686")),tolerance=1e-7)
  
  # check community membership
  condor.red.memb <- condor.object$red.mem
  condor.blue.memb <- condor.object$blue.memb
  
  expect_equal(condor.object$red.memb, condor.red.memb)
  expect_equal(condor.object$blue.memb, condor.blue.memb)
  
  # check modularity contribution of each node's
  condor.object <- condor.qscore(condor.object)
  q_women <- condor.object$qscores$red.qscore
  expect_equal(condor.object$qscores$red.qscore, q_women)
  
  # check condor.core.enrich
  out <- suppressWarnings(condor.core.enrich(c("Alice","Mary"),q=q_women,perm=T,plot.hist=T))
  expect_equal(out$analytical.pvals[1,1],0.6065307, tolerance=1e-7)
  graphics.off()
  
  # create a new condor.object aka condor.object2 for following tests:
  condor.object2<- create.condor.object(elist)
  
  # randomly assign blues to their own community
  T0 <- data.frame(nodes=blues,coms=1:4)
  
  # check matrix modularity
  condor.object3 <- condor.matrix.modularity(condor.object2,T0=T0)
  expect_equal(condor.object3$modularity,as.numeric(c("0.198347107438017","0.231404958677686","0.231404958677686"),tolerance=1e-7))
  
  # check condor.modularity.max
  condor.object4 <- condor.modularity.max(condor.object2,T0=T0)
  expect_equal(condor.object4$modularity,as.numeric(c("0.140495867768595","0.140495867768595"),tolerance=1e-7))
  

  
  # check community plot
  condor.object5 <- condor.cluster(condor.object2,project = F)
  expect_error(condor.plot.communities(condor.object5,color_list=c("darkgreen","darkorange"),
                                       point.size=2, xlab="Women",ylab="Men"),NA)

  graphics.off()
  
  
  # check heatmap
  
   data("small1976")
   condor.object6 <- create.condor.object(small1976)
   condor.object6 <- condor.cluster(condor.object6, project=FALSE)
   par(mar=c(1,1,1,1))
   expect_error(condor.plot.heatmap(condor.object6),NA)
   graphics.off()
  
  
  
 
})
