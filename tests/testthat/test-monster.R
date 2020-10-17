context("test MONSTER result")

test_that("MONSTER function works", {
  
  system('curl -O https://netzoo.s3.us-east-2.amazonaws.com/netZooR/unittest_datasets/testDatasetMonster.RData')
  load("./testDatasetMonster.RData")
  data("yeast")
  design <- c(rep(0,20),rep(NA,10),rep(1,20))
  yeast$exp.cc[is.na(yeast$exp.cc)] <- mean(as.matrix(yeast$exp.cc),na.rm=T)
  # monster result
  # caused server build failed for unknown reason.
  expect_equal(monster(yeast$exp.cc, design, yeast$motif, nullPerms=0, numMaxCores=1), monsterRes_nP0)
  
  # analyzes a bi-partite network by monster.transformation.matrix() function.
 
  #cc.net.1 <- suppressWarnings(monster.monsterNI(yeast$motif,yeast$exp.cc[1:1000,1:20])) # suppress Warning messages glm.fit: fitted probabilities numerically 0 or 1 occurred
  #cc.net.2 <- suppressWarnings(monster.monsterNI(yeast$motif,yeast$exp.cc[1:1000,31:50]))
  # error:  Error in svd(X) : infinite or missing values in 'x' 
  #expect_equal(monster.transformation.matrix(cc.net.1, cc.net.2), monsterTM)
  
  # analyzes a bi-partite network by monster.transformation.matrix() function with method "kabsch".
  # error in server
  #expect_equal(monster.transformation.matrix(cc.net.1, cc.net.2,method = "kabsch"), monsterTM_kabsch)
  
  # analyzes a bi-partite network by monster.transformation.matrix() function with method "L1".
  # to do: error  Error: $ operator not defined for this S4 class 
  #monsterTM_L1 <- monster.transformation.matrix(cc.net.1, cc.net.2,method = "L1")
  
  data("monsterRes")
  
  # Transformation matrix plot
  expect_error(monster.hcl.heatmap.plot(monsterRes), NA)
  graphics.off()
  
  # Principal Components plot of transformation matrix
  clusters <- kmeans(slot(monsterRes, 'tm'),3)$cluster 
  expect_error(monster.transitionPCAPlot(monsterRes, title="PCA Plot of Transition - Cell Cycle vs Stress Response", 
                            clusters=clusters), NA)
  
  graphics.off()
  
  # plot the transition matrix as a network
  expect_error(monster.transitionNetworkPlot(monsterRes), NA)
  graphics.off()
  
  # plots the Off diagonal mass of an observed Transition Matrix compared to a set of null TMs
  expect_error(monster.dTFIPlot(monsterRes), NA)
  
  # Calculate p-values for a tranformation matrix
  expect_equal(monster.calculate.tm.p.values(monsterRes), monster_tm_pval)
  
  # Bipartite Edge Reconstruction from Expression data with method = "pearson":
  # error here:  Error in rownames(expr.data) %in% tfNames : object 'tfNames' not found 
  cc.net_pearson <- monster.monsterNI(yeast$motif, yeast$exp.cc, method = "pearson", score="na")
  
  # Bipartite Edge Reconstruction from Expression data with other methods:
  expect_equal(monster.monsterNI(yeast$motif, yeast$exp.cc, method = "cd"), cc.net_cd, tolerance = 1e-3)
  
  # Bipartite Edge Reconstruction from Expression data (composite method with direct/indirect)
  monsterRes_bereFull<- monster.bereFull(yeast$motif, yeast$exp.cc, alpha=.5)
  expect_equal(monsterRes_bereFull[1,1], 105770, tolerance=1e-7)
  
  # summarizes the results of a MONSTER analysis
  expect_error(monster.print.monsterAnalysis(monsterRes),NA)

})



