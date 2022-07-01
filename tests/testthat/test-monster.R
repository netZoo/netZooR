context("test MONSTER result")

test_that("MONSTER function works", {
  
  system('curl -O https://netzoo.s3.us-east-2.amazonaws.com/netZooR/unittest_datasets/testDatasetMonster.RData')
  load("./testDatasetMonster.RData")
  data("yeast")
  design <- c(rep(0,20),rep(NA,10),rep(1,20))
  yeast$exp.cc[is.na(yeast$exp.cc)] <- mean(as.matrix(yeast$exp.cc),na.rm=T)
  # monster result
  expect_equal(monster(yeast$exp.cc, design, yeast$motif, nullPerms=0, numMaxCores=1, alphaw=1), monsterRes_nP0)
  
  # analyzes a bi-partite network by monster.transformation.matrix() function.
  cc.net.1 <- suppressWarnings(monsterMonsterNI(yeast$motif,yeast$exp.cc[1:1000,1:20])) # suppress Warning messages glm.fit: fitted probabilities numerically 0 or 1 occurred
  cc.net.2 <- suppressWarnings(monsterMonsterNI(yeast$motif,yeast$exp.cc[1:1000,31:50]))
  expect_equal(monsterTransformationMatrix(cc.net.1, cc.net.2), monsterTM, tolerance = 2e-3)
  
  # analyzes a bi-partite network by monsterTransformationMatrix() function with method "kabsch".
  expect_equal(monsterTransformationMatrix(cc.net.1, cc.net.2,method = "kabsch"), monsterTM_kabsch, tolerance = 2e-3)
  
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
  
  # Calculate p-values for a tranformation matrix
  expect_equal(monsterCalculateTmPValues(monsterRes), monster_tm_pval)
  
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



