context("test MONSTER result")

test_that("MONSTER function works", {
  
  netzoo_download("unittest_datasets/testDatasetMonster.RData")
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
  monsterRes@logging <- FALSE
  
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
  
  stats <- monsterCalculateTmStats(monsterRes)
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
  cc.net_pearson <- monsterMonsterNI(yeast$motif, yeast$exp.cc, method = "pearson", score="na",
                                     logging = FALSE)
  
  # Bipartite Edge Reconstruction from Expression data with other methods:
  expect_equal(monsterMonsterNI(yeast$motif, yeast$exp.cc, method = "cd", logging = FALSE), cc.net_cd, tolerance = 1e-3)
  
  # Bipartite Edge Reconstruction from Expression data (composite method with direct/indirect)
  monsterRes_bereFull<- monsterBereFull(yeast$motif, yeast$exp.cc, alpha=.5)
  monsterRes_bereFull2<- monsterMonsterNI(yeast$motif, yeast$exp.cc, method = "BERE", alphaw = 0.5)
  expect_equal(monsterRes_bereFull[1,1], 105770, tolerance=1e-7)
  expect_equal(monsterRes_bereFull2[1,1], 105770, tolerance=1e-7)
  
  # summarizes the results of a MONSTER analysis
  expect_error(monsterPrintMonsterAnalysis(monsterRes),NA)
  
  # Test that MONSTER runs without error (buildNet mode).
  testthat::expect_no_error(monster(expr = yeast$exp.cc, design = design, motif = yeast$motif, 
                                    nullPerms=3, ni_method = "BERE", mode = "buildNet",
                                    method = "ols", remove.diagonal = TRUE, logging = FALSE))
  testthat::expect_no_error(monster(expr = yeast$exp.cc, design = design, motif = yeast$motif, 
                                    nullPerms=3, ni_method = "BERE", mode = "buildNet",
                                    method = "ols", remove.diagonal = FALSE, logging = FALSE))
  testthat::expect_no_error(monster(expr = yeast$exp.cc, design = design, motif = yeast$motif, 
                                    nullPerms=3, ni_method = "BERE", mode = "buildNet",
                                    method = "kabsch", remove.diagonal = TRUE, logging = FALSE))
  # Running time for this method is too long.
  #testthat::expect_no_error(monster(expr = yeast$exp.cc, design = design, motif = yeast$motif, 
  #                                  nullPerms=3, ni_method = "BERE", mode = "buildNet",
  #                                  method = "L1", remove.diagonal = TRUE, logging = FALSE))
  testthat::expect_no_error(monster(expr = yeast$exp.cc, design = design, motif = yeast$motif, 
                                    nullPerms=3, ni_method = "BERE", mode = "buildNet",
                                    method = "orig", remove.diagonal = TRUE, logging = FALSE))
  testthat::expect_no_error(monster(expr = yeast$exp.cc, design = design, motif = yeast$motif, 
                                    nullPerms=3, ni_method = "pearson", mode = "buildNet",
                                    method = "ols", remove.diagonal = TRUE, logging = FALSE))
  testthat::expect_error(monster(expr = yeast$exp.cc, design = design, motif = yeast$motif, 
                                    nullPerms=3, ni_method = "cd", mode = "buildNet",
                                    method = "ols", remove.diagonal = TRUE, logging = FALSE),
                         paste("Supported values for ni_method are BERE, pearson, and lda.",
                               "cd is invalid."))
  testthat::expect_no_error(monster(expr = yeast$exp.cc, design = design, motif = yeast$motif, 
                                    nullPerms=3, ni_method = "lda", mode = "buildNet",
                                    method = "ols", remove.diagonal = TRUE, logging = FALSE))
  
  # Test that MONSTER runs without error (regNet mode.)
  set.seed(123)
  exp_grn <- matrix(data = rnorm(50, mean = 1), ncol = 10, nrow = 5)
  control_grn <- matrix(data = rnorm(50, mean = 1), ncol = 10, nrow = 5)
  colnames(exp_grn) <- paste0('gene', 1:10)
  colnames(control_grn) <- paste0('gene', 1:10)
  rownames(exp_grn) <- paste0('tf', 1:5)
  rownames(control_grn) <- paste0('tf', 1:5)
  combdf <- as.data.frame(cbind(control_grn, exp_grn))
  combdes <- c(rep(0, ncol(control_grn)), rep(1, ncol(exp_grn)))
  testthat::expect_no_error(monster(expr = combdf,
                                           design = combdes,
                                           motif = NA,
                                           mode = 'regNet',
                                           nullPerms = 3,
                                           numMaxCores = 1, logging = FALSE))
  
  # Test that MONSTER with null generated from control data gives the same
  # result as MONSTER with null loaded from a file.
  set.seed(1)
  networks <- as.data.frame(matrix(runif(50), ncol = 10, nrow = 5))
  colnames(networks) <- paste0('gene', 1:10)
  rownames(networks) <- paste0('tf', 1:5)
  networksIn <- networks
  networksIn$tf <- rownames(networks)
  str(networks)
  nullNetworks <- reshape2::melt(
    as.matrix(networksIn[,combdes == 0])
  )
  colnames(nullNetworks) <- c("tf", "gene", "score")
  nullNetworks$tf <- as.character(nullNetworks$tf)
  nullNetworks$gene <- as.character(nullNetworks$gene)
  for(i in 3:12){
    nullNetworks[,i] <- runif(50, min = 0, max = 0.1)
  }
  monsterControlNull <- monster(expr = combdf,
                                    design = combdes,
                                    motif = NA,
                                    mode = 'regNet', nullNetworks = nullNetworks,
                                    nullPerms = 10, nullModelType = "nullNetwork",
                                    numMaxCores = 1, logging = FALSE)
  # monsterControlNullMultiCore <- monster(expr = combdf,
  #                               design = combdes,
  #                               motif = NA,
  #                               mode = 'regNet', nullNetworks = nullNetworks,
  #                               nullPerms = 10, nullModelType = "nullNetwork",
  #                               numMaxCores = 3, logging = FALSE)
  nullNetworksOffset <- nullNetworks
  nullNetworksOffset[,3:12] <- nullNetworksOffset[,3:12] + 5
  nullNetworks <- nullNetworksOffset
  rhdf5::h5createFile("nullNetworks.h5")
  rhdf5::h5createGroup("nullNetworks.h5", "matrices")
  
  # Make a list of networks.
  for(i in seq(3, 12)){
    wide <- unclass(xtabs(score ~ tf + gene, data = data.frame(tf = nullNetworksOffset[,1],
                                                      gene = nullNetworksOffset[,2],
                                                      score = nullNetworksOffset[,i])))
    if(i == 3){
      rhdf5::h5write(rownames(wide), "nullNetworks.h5", paste0("matrices/tfs"))
      rhdf5::h5write(colnames(wide), "nullNetworks.h5", paste0("matrices/genes"))
    }
    rhdf5::h5write(wide, "nullNetworks.h5", paste0("matrices/", as.character(i-2)))
  }
  
  # Run MONSTER.
  monsterFileNull <- monster(expr = combdf,
                                design = combdes,
                                motif = NA,
                                mode = 'regNet', nullNetworks = "nullNetworks.h5",
                                nullPerms = 10, nullModelType = "nullNetwork",
                                numMaxCores = 1, logging = FALSE)
  # monsterFileNullMultiCore <- monster(expr = combdf,
  #                                        design = combdes,
  #                                        motif = NA,
  #                                        mode = 'regNet', nullNetworks = "nullNetworks.h5",
  #                                        nullPerms = 10, nullModelType = "nullNetwork",
  #                                        numMaxCores = 3, logging = FALSE)
  expect_equal(c(monsterFileNull@tm), c(monsterControlNull@tm))
  for(i in 1:10){
    expect_equal(c(monsterFileNull@nullTM[[i]]), c(monsterControlNull@nullTM[[i]]))
  }
  expect_equal(monsterFileNull@numGenes, monsterControlNull@numGenes)
  expect_equal(c(monsterFileNull@numSamples), c(monsterControlNull@numSamples))
  expect_equal(monsterFileNull@logging, monsterControlNull@logging)
  
   monsterControlNull <- monster(expr = combdf,
                                    design = combdes,
                                    motif = NA,
                                    mode = 'regNet', nullNetworks = nullNetworks,
                                    nullPerms = 10, nullModelType = "nullNetwork",
                                    numMaxCores = 1, logging = FALSE)
  # monsterControlNullMultiCore <- monster(expr = combdf,
  #                               design = combdes,
  #                               motif = NA,
  #                               mode = 'regNet', nullNetworks = nullNetworks,
  #                               nullPerms = 10, nullModelType = "nullNetwork",
  #                               numMaxCores = 3, logging = FALSE)
  nullNetworksOffset <- nullNetworks
  nullNetworksOffset[,3:12] <- nullNetworksOffset[,3:12] + 5
  nullNetworks <- nullNetworksOffset
  rhdf5::h5createFile("nullNetworks.h5")
  rhdf5::h5createGroup("nullNetworks.h5", "matrices")
  
  # Make a list of networks.
  for(i in seq(3, 12)){
    wide <- unclass(xtabs(score ~ tf + gene, data = data.frame(tf = nullNetworksOffset[,1],
                                                      gene = nullNetworksOffset[,2],
                                                      score = nullNetworksOffset[,i])))
    if(i == 3){
      rhdf5::h5write(rownames(wide), "nullNetworks.h5", paste0("matrices/tfs"))
      rhdf5::h5write(colnames(wide), "nullNetworks.h5", paste0("matrices/genes"))
    }
    rhdf5::h5write(wide, "nullNetworks.h5", paste0("matrices/", as.character(i-2)))
  }
  
  # Run MONSTER.
  monsterFileNull <- monster(expr = combdf,
                                design = combdes,
                                motif = NA,
                                mode = 'regNet', nullNetworks = "nullNetworks.h5",
                                nullPerms = 10, nullModelType = "nullNetwork",
                                numMaxCores = 1, logging = FALSE)
  # monsterFileNullMultiCore <- monster(expr = combdf,
  #                                        design = combdes,
  #                                        motif = NA,
  #                                        mode = 'regNet', nullNetworks = "nullNetworks.h5",
  #                                        nullPerms = 10, nullModelType = "nullNetwork",
  #                                        numMaxCores = 3, logging = FALSE)
  expect_equal(c(monsterFileNull@tm), c(monsterControlNull@tm))
  for(i in 1:10){
    expect_equal(c(monsterFileNull@nullTM[[i]]), c(monsterControlNull@nullTM[[i]]))
  }
  expect_equal(monsterFileNull@numGenes, monsterControlNull@numGenes)
  expect_equal(c(monsterFileNull@numSamples), c(monsterControlNull@numSamples))
  expect_equal(monsterFileNull@logging, monsterControlNull@logging)
  
  # Get the gene transition matrix instead.
  monsterControlNullGene <- monster(expr = combdf,
                                design = combdes,
                                motif = NA,
                                mode = 'regNet', nullNetworks = nullNetworks,
                                nullPerms = 10, nullModelType = "nullNetwork",
                                numMaxCores = 1, logging = FALSE, by.tfs = FALSE,
                                method = "orig")
  
  # Run MONSTER.
  monsterFileNullGene <- monster(expr = combdf,
                             design = combdes,
                             motif = NA,
                             mode = 'regNet', nullNetworks = "nullNetworks.h5",
                             nullPerms = 10, nullModelType = "nullNetwork",
                             numMaxCores = 1, logging = FALSE, by.tfs = FALSE,
                             method = "orig")
  expect_equal(c(monsterFileNullGene@tm), c(monsterControlNullGene@tm))
  for(i in 1:10){
    expect_equal(c(monsterFileNullGene@nullTM[[i]]), c(monsterControlNullGene@nullTM[[i]]))
  }
  expect_equal(monsterFileNullGene@numGenes, monsterControlNullGene@numGenes)
  expect_equal(c(monsterFileNullGene@numSamples), c(monsterControlNullGene@numSamples))
  expect_equal(monsterFileNullGene@logging, monsterControlNullGene@logging)
  
  # Remove data file.
  unlink("nullNetworks.h5")
  unlink("./testDatasetMonster.RData")
})

test_that("monsterCalculateTmStatsPerCell function works", {
  
  # Initialize data.
  set.seed(1)
  monsterTMNeg <- matrix(rnorm(mean = 0, sd = 1, n = 100), nrow = 10)
  rownames(monsterTMNeg) <- paste0("TF", 1:10)
  colnames(monsterTMNeg) <- paste0("TF", 1:10)
  monsterTMPos <- monsterTMNeg
  monsterTMPos[5,6] <- 200
  monsterTMPos[1,2] <- 300
  monsterNullTM <- lapply(1:1000, function(i){
    retval <- matrix(rnorm(mean = 0, sd = 1, n = 100), nrow = 10)
    rownames(retval) <- paste0("TF", 1:10)
    colnames(retval) <- paste0("TF", 1:10)
    return(retval)
  })
  monsterNegObj <- new("monsterAnalysis", tm = monsterTMNeg, nullTM = monsterNullTM,
                       numGenes = 1000, numSamples = 1000, logging = FALSE)
  monsterPosObj <- new("monsterAnalysis", tm = monsterTMPos, nullTM = monsterNullTM,
                       numGenes = 1000, numSamples = 1000, logging = FALSE)

  # Error if the data are not of the correct type.
  expect_error(monsterCalculateTmStatsPerCell(monsterObj = NULL, method = "z-score"),
               "Input must be an object of class 'monsterAnalysis'.")
  expect_error(monsterCalculateTmStatsPerCell(monsterObj = monsterNegObj, method = "google"),
               "Valid methods include 'z-score' and 'non-parametric'. Invalid method: google")
  
  # Ensure that nothing is significant in the negative case.
  monsterNegStats <- monsterCalculateTmStatsPerCell(monsterObj = monsterNegObj, method = "z-score")
  expect_all_true(c(monsterNegStats$p.adj[which(!is.na(monsterNegStats$p.adj))]) > 0.05)
  monsterNegStats <- monsterCalculateTmStatsPerCell(monsterObj = monsterNegObj, method = "non-parametric")
  expect_all_true(c(monsterNegStats$p.adj[which(!is.na(monsterNegStats$p.adj))]) > 0.05)
  
  # Ensure that the right values are significant in the positive case.
  monsterPosStats <- monsterCalculateTmStatsPerCell(monsterObj = monsterPosObj, method = "z-score")
  expect_lt(monsterPosStats$p.adj[5,6], 0.05)
  expect_lt(monsterPosStats$p.adj[1,2], 0.05)
  monsterPosStats$p.adj[5,6] <- 1
  monsterPosStats$p.adj[1,2] <- 1
  expect_all_true(c(monsterPosStats$p.adj[which(!is.na(monsterPosStats$p.adj))]) > 0.05)
  monsterPosStats2 <- monsterCalculateTmStatsPerCell(monsterObj = monsterPosObj, method = "non-parametric")
  expect_lt(monsterPosStats2$p.adj[5,6], 0.05)
  expect_lt(monsterPosStats2$p.adj[1,2], 0.05)
  monsterPosStatsMod2 <- monsterPosStats2
  monsterPosStatsMod2$p.adj[5,6] <- 1
  monsterPosStatsMod2$p.adj[1,2] <- 1
  expect_all_true(c(monsterPosStatsMod2$p.adj[which(!is.na(monsterPosStatsMod2$p.adj))]) > 0.05)
  
  # Set diagonals to 0 and check for equality.
  monsterPosObj2 <- monsterPosObj
  diag(monsterPosObj2@tm) <- 0
  for(i in 1:length(monsterPosObj2@nullTM)){
    diag(monsterPosObj2@nullTM[[i]]) <- 0
  }
  monsterPosStatsNoDiag <- monsterCalculateTmStatsPerCell(monsterObj = monsterPosObj2, method = "non-parametric")
  expect_equal(c(monsterPosStats2$p.adj), c(monsterPosStatsNoDiag$p.adj))
})

test_that("GeneratePermutationNull function works", {
  
  # Put together MONSTER inputs.
  set.seed(1)
  networks <- as.data.frame(matrix(runif(24), nrow = 3))
  rownames(networks) <- paste0("tf", 1:3)
  colnames(networks) <- rep(paste0("gene", 1:4), 2)
  design <- c(rep(0, 4), rep(1, 4))
  motif <- data.frame(tf = c("tf1", rep("tf2", 4), "tf1", "tf2"), 
                      gene = c("gene2", "gene1", "gene2", "gene3", "gene4", "tf2", "tf1"), 
                      score = c(1, 1, rep(0, 3), 1, 1))
  
  # Test that the MONSTER function does not allow GeneratePermutationNull
  # to be called without the correct inputs.
  expect_error(monster(expr = networks, design = design, motif = NA, nullModelType = "aj;fdaipofhnweapiofn"),
               "Only 'permutation' and 'nullNetwork' are permitted values for nullModelType.")
  
  # Test that the permutation null has the same format and that
  # the values are not equal for regNet case.
  permNull <- GeneratePermutationNull(expr = networks, iterations = 10, mode = "regNet", logging = FALSE)
  expect_equal(colnames(permNull[[5]]), colnames(networks))
  expect_equal(rownames(permNull[[5]]), rownames(networks))
  expect_true(!identical(permNull[[5]], networks))
  
  # Test the same for buildNet case. Note that this requires a different input format.
  expression <- as.data.frame(matrix(runif(60), nrow = 6))
  rownames(expression) <- c("tf1", "tf2", paste0("gene", 1:4))
  colnames(expression) <- paste0("samp", 1:10)
  design2 <- c(rep(0, 5), rep(1, 5))
  permNull <- GeneratePermutationNull(expr = expression, iterations = 10, mode = "buildNet", logging = FALSE)
  expect_equal(colnames(permNull[[5]]), colnames(expression))
  expect_equal(rownames(permNull[[5]]), rownames(expression))
  expect_true(!identical(permNull[[5]], expression))
  
  # Test that MONSTER runs without error.
  expect_no_error(monster(expr = networks, design = design, motif = motif, nullModelType = "permutation",
                                          mode = "regNet", nullPerms = 10, logging = FALSE))
  # We cannot test multi-core processes in the build environment.
  #expect_no_error(suppressWarnings(monster(expr = networks, design = design, motif = motif, nullModelType = "permutation",
  #                        mode = "regNet", numMaxCores = 3, nullPerms = 10)))
  expect_no_error(monster(expr = expression, design = design2, motif = motif, nullModelType = "permutation",
                          mode = "buildNet", nullPerms = 10, logging = FALSE))
  # We cannot test multi-core processes in the build environment.
  #expect_no_error(suppressWarnings(monster(expr = expression, design = design2, motif = motif, nullModelType = "permutation",
  #                                         mode = "buildNet", numMaxCores = 3, nullPerms = 10)))
})

test_that("GenerateNullFromControl function works", {
  
  # Set up inputs.
  set.seed(1)
  networks <- as.data.frame(matrix(runif(24), nrow = 3))
  rownames(networks) <- paste0("tf", 1:3)
  colnames(networks) <- rep(paste0("gene", 1:4), 2)
  design <- c(rep(0, 4), rep(1, 4))
  networksIn <- networks
  networksIn$tf <- rownames(networks)
  nullNetworks <- reshape2::melt(
    networksIn[,design == 0],
    id.vars = "tf",
    variable.name = "gene",
    value.name = "score"
  )
  nullNetworks$gene <- as.character(nullNetworks$gene)
  for(i in 3:12){
    nullNetworks[,i] <- runif(12, min = 0, max = 0.1)
  }
  
  # Test that the MONSTER function does not allow GenerateNullFromControl
  # to be called without the correct inputs.
  nullNetworksWrong1 <- nullNetworksWrong2 <- nullNetworksWrong3 <- nullNetworksWrong4 <- nullNetworks
  nullNetworksWrong1[,1] <- 0
  nullNetworksWrong2[,2] <- 0
  nullNetworksWrong3[,3] <- rep("", nrow(nullNetworksWrong3))
  nullNetworksWrong4[1,1] <- "afpujewipaojefwp"
  
  expect_error(monster(expr = networks, design = design, motif = NA, nullModelType = "nullNetwork", mode = "buildNet",
                       nullNetworks = nullNetworks),
               "Cannot run 'nullNetwork' model type in buildNet mode")
  expect_error(monster(expr = networks, design = design, motif = NA, nullModelType = "nullNetwork", mode = "regNet",
                       nullNetworks = NA),
               paste("Must provide null networks when null model type is 'nullNetwork'. These can",
                     "be generated using the function GenerateNullPANDADistribution() if you are",
                     "using PANDA networks."), fixed = TRUE)
  expect_error(monster(expr = networks, design = design, motif = NA, nullModelType = "nullNetwork", mode = "regNet",
                       nullNetworks = nullNetworksWrong1),
               paste("The null network provided must have source and target nodes in the first two columns",
                     "and numeric scores in all remaining columns."))
  expect_error(monster(expr = networks, design = design, motif = NA, nullModelType = "nullNetwork", mode = "regNet",
                       nullNetworks = nullNetworksWrong2),
               paste("The null network provided must have source and target nodes in the first two columns",
                     "and numeric scores in all remaining columns."))
  expect_error(monster(expr = networks, design = design, motif = NA, nullModelType = "nullNetwork", mode = "regNet",
                       nullNetworks = nullNetworksWrong3),
               paste("The null network provided must have source and target nodes in the first two columns",
                     "and numeric scores in all remaining columns."))
  expect_error(monster(expr = networks, design = design, motif = NA, nullModelType = "nullNetwork", mode = "regNet",
                       nullNetworks = nullNetworksWrong4),
               "The node lists differ between the input networks and the provided null networks.")
  
  # Test that the control is still the same and the format is correct.
  null <- GenerateNullFromControl(concatNet = networks, iterations = 10, design = design, nullNetworks = nullNetworks,
                                  logging = FALSE)
  expect_equal(colnames(null[[5]]), colnames(networks))
  expect_equal(rownames(null[[5]]), rownames(networks))
  expect_true(!identical(null[[5]], networks))
  expect_true(identical(null[[5]][,design == 0], as.matrix(networks[,design == 0])))
  
  # Test that MONSTER runs without error.
  expect_no_error(monster(expr = networks, design = design, motif = NA, nullModelType = "nullNetwork",
                          mode = "regNet", nullPerms = 10, nullNetworks = nullNetworks, logging = FALSE))
  # We cannot test multi-core processes in the build environment.
  #expect_no_error(suppressWarnings(monster(expr = networks, design = design, motif = NA, nullModelType = "nullNetwork",
  #                                         mode = "regNet", numMaxCores = 3, nullPerms = 10, nullNetworks = nullNetworks)))
  
  # Test that a warning is issued when the number of iterations does not match the number of null networks and that the
  # number of generated null networks matches the input null.
  numNets <- ncol(nullNetworks) - 2
  expect_warning(GenerateNullFromControl(concatNet = networks, iterations = 100, design = design, nullNetworks = nullNetworks,
                                         logging = FALSE),
                 paste("Warning: You have specified", 100, "iterations but",
                       "provided", numNets, "null networks. MONSTER",
                       "will generate", numNets, "iterations instead."))
  null <- suppressWarnings(GenerateNullFromControl(concatNet = networks, iterations = 100, design = design, nullNetworks = nullNetworks,
                                                   logging = FALSE))
  expect_equal(length(null), numNets)
  
})

test_that('domonster runs on toy PANDA data', {
  set.seed(123)
  exp_grn <- matrix(data = rnorm(50, mean = 1), ncol = 10, nrow = 5)
  control_grn <- matrix(data = rnorm(50, mean = 1), ncol = 10, nrow = 5)
  colnames(exp_grn) <- paste0('gene', 1:10)
  colnames(control_grn) <- paste0('gene', 1:10)
  rownames(exp_grn) <- paste0('tf', 1:5)
  rownames(control_grn) <- paste0('tf', 1:5)
  
  testthat::expect_no_error(domonster(exp_grn, control_grn, numMaxCores = 1, logging = FALSE))
  
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

test_that("monsterCheckDataType() converts data.frame to matrix", {
  df <- data.frame(a = 1:3, b = 4:6)
  result <- monsterCheckDataType(df)
  expect_true(is.matrix(result))
  expect_equal(dim(result), c(3, 2))
})

test_that("monsterCheckDataType() passes through matrix unchanged", {
  mat <- matrix(1:6, nrow = 2)
  result <- monsterCheckDataType(mat)
  expect_true(is.matrix(result))
  expect_equal(result, mat)
})

test_that("monsterCheckDataType() converts ExpressionSet", {
  expr_mat <- matrix(rnorm(20), nrow = 4,
                     dimnames = list(paste0("g", 1:4), paste0("s", 1:5)))
  eset <- Biobase::ExpressionSet(assayData = expr_mat)
  result <- monsterCheckDataType(eset)
  expect_true(is.matrix(result))
  expect_equal(dim(result), c(4, 5))
})

test_that("monsterCheckDataType() errors on invalid input", {
  expect_error(monsterCheckDataType("foo"), "must be a data.frame")
  expect_error(monsterCheckDataType(42), "must be a data.frame")
  expect_error(monsterCheckDataType(list(1, 2, 3)), "must be a data.frame")
})

test_that("kabsch() returns square transformation matrix", {
  set.seed(42)
  n <- 20
  p <- 3
  P <- matrix(rnorm(n * p), nrow = n, ncol = p)
  colnames(P) <- paste0("TF", 1:p)
  
  # Create Q as a rotation of P with some noise
  Q <- P %*% matrix(c(0, -1, 0, 1, 0, 0, 0, 0, 1), 3, 3) + matrix(rnorm(n * p, sd = 0.1), nrow = n, ncol = p)
  colnames(Q) <- paste0("TF", 1:p)
  
  W <- kabsch(P, Q)
  expect_true(is.matrix(W))
  expect_equal(nrow(W), p)
  expect_equal(ncol(W), p)
  expect_true(all(is.finite(W)))
})

test_that("monsterTransformationMatrix() works with ols method", {
  data("yeast")
  yeast$exp.cc[is.na(yeast$exp.cc)] <- mean(as.matrix(yeast$exp.cc), na.rm = TRUE)
  cc.net.1 <- suppressWarnings(monsterMonsterNI(yeast$motif, yeast$exp.cc[1:1000, 1:20]))
  cc.net.2 <- suppressWarnings(monsterMonsterNI(yeast$motif, yeast$exp.cc[1:1000, 31:50]))
  
  tm <- monsterTransformationMatrix(cc.net.1, cc.net.2, method = "ols")
  expect_true(is.matrix(tm))
  expect_equal(nrow(tm), ncol(tm))
})

test_that("monsterTransformationMatrix() works with kabsch method", {
  data("yeast")
  yeast$exp.cc[is.na(yeast$exp.cc)] <- mean(as.matrix(yeast$exp.cc), na.rm = TRUE)
  cc.net.1 <- suppressWarnings(monsterMonsterNI(yeast$motif, yeast$exp.cc[1:1000, 1:20]))
  cc.net.2 <- suppressWarnings(monsterMonsterNI(yeast$motif, yeast$exp.cc[1:1000, 31:50]))
  
  tm <- monsterTransformationMatrix(cc.net.1, cc.net.2, method = "kabsch")
  expect_true(is.matrix(tm))
})


