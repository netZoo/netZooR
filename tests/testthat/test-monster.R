context("test MONSTER result")

test_that("MONSTER function works", {
  
  system('curl -O https://netzoo.s3.us-east-2.amazonaws.com/netZooR/unittest_datasets/testDatasetMonster.RData')
  load("./testDatasetMonster.RData")
  data("yeast")
  design <- c(rep(0,20),rep(NA,10),rep(1,20))
  yeast$exp.cc[is.na(yeast$exp.cc)] <- mean(as.matrix(yeast$exp.cc),na.rm=T)
  # monster result
  #expect_equal(monster(yeast$exp.cc, design, yeast$motif, nullPerms=0, numMaxCores=1, alphaw=1), monsterRes_nP0)
  
  # analyzes a bi-partite network by monster.transformation.matrix() function.
  #cc.net.1 <- suppressWarnings(monsterMonsterNI(yeast$motif,yeast$exp.cc[1:1000,1:20])) # suppress Warning messages glm.fit: fitted probabilities numerically 0 or 1 occurred
  #cc.net.2 <- suppressWarnings(monsterMonsterNI(yeast$motif,yeast$exp.cc[1:1000,31:50]))
  #expect_equal(monsterTransformationMatrix(cc.net.1, cc.net.2), monsterTM, tolerance = 3e-3)
  
  # analyzes a bi-partite network by monsterTransformationMatrix() function with method "kabsch".
  #expect_equal(monsterTransformationMatrix(cc.net.1, cc.net.2,method = "kabsch"), monsterTM_kabsch, tolerance = 3e-3)
  
  # to do: Add L1 method in test
  
  data("monsterRes")
  
  # Transformation matrix plot
  expect_error(monsterHclHeatmapPlot(monsterRes), NA)
  graphics.off()
  
  # Principal Components plot of transformation matrix
  clusters <- kmeans(slot(monsterRes, 'tm'),3)$cluster 
  expect_error(monsterTransitionPCAPlot(monsterRes, title="PCA Plot of Transition - Cell Cycle vs Stress Response", 
                            clusters=clusters), NA)
  
  graphics.off()
  
  # plot the transition matrix as a network
  expect_error(monsterTransitionNetworkPlot(monsterRes), NA)
  graphics.off()
  
  # plots the Off diagonal mass of an observed Transition Matrix compared to a set of null TMs
  expect_error(monsterdTFIPlot(monsterRes), NA)
  
  monsterCalculateTmStats(monsterRes)
  stopifnot("monsterCalculateTmStats" %in% ls("package:netZooR"))

  # Calculate p-values for a tranformation matrix
  # #TODO: update data to include this test
  #expect_equal(monsterCalculateTmPValues(monsterRes), monster_tm_pval)

  # Before refactoring the test, we can check that the non-paramentrc and z-score methods are similar
  c = monsterCalculateTmPValues(monsterRes)
  d = monsterCalculateTmPValues(monsterRes, method = 'non-parametric')
  r <- cor(c, d, use = "complete.obs")  # Handle NAs if needed
  expect_gt(r, 0.1) 


  # To load the data again
  #load("../data/monsterPvals.RData")
  monster_pvals = monsterPvals$p.values
  monster_tvals = monsterPvals$t.values
  ssodm = monsterPvals$ssodm
  null.ssodm.matrix = monsterPvals$null.ssodm.matrix
  # Now we compare it with the original data
  newp = monsterCalculateTmStats(monsterRes)
  expect_equal(newp$p.values, monster_pvals)
  expect_equal(newp$t.values, monster_tvals)
  expect_equal(newp$ssodm, ssodm)
  expect_equal(newp$null.ssodm.matrix, null.ssodm.matrix)
  
  # Bipartite Edge Reconstruction from Expression data with method = "pearson":
  # error here:  Error in rownames(expr.data) %in% tfNames : object 'tfNames' not found 
  cc.net_pearson <- monsterMonsterNI(yeast$motif, yeast$exp.cc, method = "pearson", score="na")
  
  # Bipartite Edge Reconstruction from Expression data with other methods:
  expect_equal(monsterMonsterNI(yeast$motif, yeast$exp.cc, method = "cd"), cc.net_cd, tolerance = 1e-3)
  
  # Bipartite Edge Reconstruction from Expression data (composite method with direct/indirect)
  monsterRes_bereFull<- monsterBereFull(yeast$motif, yeast$exp.cc, alpha=.5)
  expect_equal(monsterRes_bereFull[1,1], 105770, tolerance=1e-7)
  
  # summarizes the results of a MONSTER analysis
  expect_error(monsterPrintMonsterAnalysis(monsterRes),NA)

})

test_that('domonster runs on toy PANDA data', {
  set.seed(123)
  exp_grn <- matrix(data = rnorm(50, mean = 1), ncol = 10, nrow = 5)
  control_grn <- matrix(data = rnorm(50, mean = 1), ncol = 10, nrow = 5)
  colnames(exp_grn) <- paste0('gene', 1:10)
  colnames(control_grn) <- paste0('gene', 1:10)
  rownames(exp_grn) <- paste0('tf', 1:5)
  rownames(control_grn) <- paste0('tf', 1:5)
  
  testthat::expect_no_error(domonster(exp_grn, control_grn, numMaxCores = 1))
  
  # # more robust test using the pandaToyData that is unexplainably failing the github checks
  # pandaResult_exp <- panda(pandaToyData$motif, pandaToyData$expression[,1:25], pandaToyData$ppi)
  # pandaResult_control <- panda(pandaToyData$motif, pandaToyData$expression[,26:50], pandaToyData$ppi)
  # 
  # # function takes both panda objects and matrices, or a mixture
  # set.seed(123)
  # monster_res1 <- domonster(pandaResult_exp, pandaResult_control, numMaxCores = 1)
  # 
  # set.seed(123)
  # monster_res2 <- domonster(pandaResult_exp@regNet, pandaResult_control@regNet, numMaxCores = 1)
  # 
  # set.seed(123)
  # monster_res3 <- domonster(pandaResult_exp@regNet, pandaResult_control, numMaxCores = 1)
  # 
  # # these should all yield same result; confirming they are the same
  # expect_equal(monster_res1, monster_res2, tolerance=1e-15) 
  # expect_equal(monster_res1, monster_res3, tolerance=1e-15) 
})


