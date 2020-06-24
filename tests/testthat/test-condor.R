context("test CONDOR result")

test_that("CONDOR functions work", {
  
  load("./testDataset.RData")
  
  r = c(1,1,1,2,2,2,3,3,3,4,4);
  b = c(1,2,3,1,2,4,2,3,4,3,4);
  reds <- c("Alice","Sue","Janine","Mary")
  blues <- c("Bob","John","Ed","Hank")
  elist <- data.frame(red=reds[r], blue=blues[b])
  
  condor.object <- create.condor.object(elist)

  expect_equal(names(condor.object),c("G","edges","Qcoms","modularity","red.memb","blue.memb","qscores"))
  
  condor.object <- condor.cluster(condor.object,project = F)
  
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
  
})
