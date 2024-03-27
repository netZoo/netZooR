# unit-tests for FERRET
context("test FERRET functions")

test_that("[FERRET] LoadResults() function yields expected results", {
  # Create a dummy directory and dummy networks.
  if(!file.exists("tmp")){
    dir.create("tmp")
  }
  network1 <- data.frame(source = c("tfA", "tfB", "tfC"),
                         target = c("geneA", "geneB", "geneC"),
                         score = c(1,1,1))
  write.csv(network1, "tmp/network1.csv")
  network2 <- data.frame(source = c("tfA", "tfB", "tfC"),
                        target = c("geneA", "geneB", "geneC"),
                        score = c(1,-1,1))
  write.csv(network2, "tmp/network2.csv")
  network3 <- data.frame(source = c("tfA", "tfB", "tfC"),
                         target = c("geneA", "geneB", "geneC"),
                         score = c(0,1,1))
  write.csv(network3, "tmp/network3.csv")
  
  # Create a temporary data where everything is saved as RDS.
  dir.create("tmpRDS")
  saveRDS(network1, "tmpRDS/network1.RDS")
  
  # Create an empty directory.
  dir.create("tmpEmpty")
  
  # Load the data.
  dataPoor <- LoadResults("tmp")
  dataInhibitory <- LoadResults("tmp", "inhibitory")

  # Check that the results are of the correct type.
  expect_true(is(dataPoor, "FERRET_Results"))
  expect_true(is(dataInhibitory, "FERRET_Results"))
  
  # Check that values are correct.
  expect_identical(dataPoor@results[[1]]$source, network1$source)
  expect_identical(dataPoor@results[[2]]$target, network2$target)
  expect_equal(dataPoor@results[[3]]$score, network3$score)
  expect_equal(dataPoor@directory, "tmp")
  expect_equal(dataPoor@interpretationOfNegative, "poor")
  expect_identical(dataInhibitory@results[[1]]$source, network1$source)
  expect_identical(dataInhibitory@results[[2]]$target, network2$target)
  expect_equal(dataInhibitory@results[[3]]$score, network3$score)
  expect_equal(dataInhibitory@directory, "tmp")
  expect_equal(dataInhibitory@interpretationOfNegative, "inhibitory")
  
  # Create a network where the scores are not numeric.
  networkChar <- network1
  networkChar$score <- c("Hello", "world", "!")
  dir.create("tmpChar")
  write.csv(networkChar, "tmpChar/network1.csv")
  
  # Create a network without scores.
  networkNoScore <- network1[,1:2]
  dir.create("tmpNoScore")
  write.csv(networkNoScore, "tmpNoScore/network1.csv")
  
  # Check that the correct errors are thrown.
  suppressWarnings({
    expect_error(LoadResults("($&@_(*&$_@"), "Invalid directory: ($&@_(*&$_@", fixed = TRUE)
    expect_error(LoadResults("tmpChar"), "File has an invalid format: network1.csv")
    expect_error(LoadResults("tmpNoScore"), "File has an invalid format: network1.csv")
    expect_error(LoadResults("tmpRDS"), "File has an invalid format: network1.RDS")
    expect_error(LoadResults("tmpEmpty"), "The directory tmpEmpty is empty!")
    expect_error(LoadResults("tmp", "OogaBooga"), paste("Valid values for interpretationOfNegative are 'poor' or 'inhibitory.'",
                                                        "You entered: OogaBooga"))
  })
  
  # Delete directories created by this function.
  unlink("tmp", recursive = TRUE)
  unlink("tmpRDS", recursive = TRUE)
  unlink("tmpEmpty", recursive = TRUE)
  unlink("tmpChar", recursive = TRUE)
  unlink("tmpNoScore", recursive = TRUE)
})
test_that("[FERRET] BuildComparisonObject() function yields expected results", {
  # Create a FERRET_Results object with dummy networks.
  network1 <- data.frame(source = c("tfA", "tfB", "tfC"),
                         target = c("geneA", "geneB", "geneC"),
                         score = c(1,1,1))
  network2 <- data.frame(source = c("tfA", "tfB", "tfC"),
                         target = c("geneA", "geneB", "geneC"),
                         score = c(1,-1,1))
  network3 <- data.frame(source = c("tfA", "tfB", "tfC"),
                         target = c("geneA", "geneB", "geneC"),
                         score = c(0,1,1))
  resultList <- list(network1, network2, network3)
  names(resultList) <- c("network1", "network2", "network3")
  result <- methods::new("FERRET_Results",results = resultList, directory = "tmp")
  
  # Build a comparison object to test the ideal case. Check that it contains
  # the objects as expected.
  comparison <- BuildComparisonObject("network1", "network2", "network3", result)
  expect_true(is(comparison, "FERRET_Comparisons"))
  expect_equal(length(comparison@ingroup), 1)
  expect_equal(length(comparison@outgroup), 1)
  expect_true(is(comparison@ingroup[[1]], "FERRET_Comparison"))
  expect_true(is(comparison@outgroup[[1]], "FERRET_Comparison"))
  expect_equal(comparison@ingroup[[1]]@source, "network1")
  expect_equal(comparison@ingroup[[1]]@target, "network2")
  expect_equal(comparison@outgroup[[1]]@source, "network1")
  expect_equal(comparison@outgroup[[1]]@target, "network3")
  
  # Build a comparison object to test the cases in which the inputs are not
  # of the expected types.
  incorrectTypeMsg <- paste("SourceNetwork must be a name (string). ingroupToCompare and outgroupToCompare must be",
                            "a list of names (strings). results must be an object of type FERRET_Results.")
  expect_error(BuildComparisonObject(NULL, "network2", "network3", result), incorrectTypeMsg, fixed = TRUE)
  expect_error(BuildComparisonObject("network1", NULL, "network3", result), incorrectTypeMsg, fixed = TRUE)
  expect_error(BuildComparisonObject("network1", "network2", NULL, result), incorrectTypeMsg, fixed = TRUE)
  expect_error(BuildComparisonObject("network1", "network2", "network3", NULL), incorrectTypeMsg, fixed = TRUE)
  
  # Build a comparison object to test the cases in which a network is being
  # compared to itself or a network is in both the in-group and the out-group.
  sourceOverlapMsg <- "The source network cannot be in the in-group or the out-group."
  overlapMsg <- "The in-group and the out-group cannot overlap."
  expect_error(BuildComparisonObject("network1", c("network1", "network2"), 
                                     "network3", result), sourceOverlapMsg, fixed = TRUE)
  expect_error(BuildComparisonObject("network1", "network2", 
                                     c("network1", "network3"), result), sourceOverlapMsg, fixed = TRUE)
  expect_error(BuildComparisonObject("network1", c("network3", "network2"), 
                                     "network3", result), overlapMsg, fixed = TRUE)
})
test_that("[FERRET] AUCTrapezoid() function yields expected results",{
  # Test for a straight line.
  x <- c(0, 0.5, 1)
  y <- c(0, 0.5, 1)
  expect_equal(AUCTrapezoid(x,y), 0.5)
  
  # Test for a straight line with y-intercept = 0.5.
  x <- c(0, 0.5, 1)
  y <- c(0.5, 0.75, 1)
  expect_equal(AUCTrapezoid(x,y), 0.75)
  
  # Test for a straight line with x-intercept = 0.5.
  x <- c(0.5, 0.75, 1)
  y <- c(0, 0.5, 1)
  expect_equal(AUCTrapezoid(x,y), 0.25)
  
  # Test for a curve with x-intercept = 0.5.
  x <- c(0.5, 0.75, 1)
  y <- c(0, 0.75, 1)
  expect_equal(AUCTrapezoid(x,y), 0.3125)
  
  # Test for a straight line where max(X) < max(Y).
  x <- c(0, 0.25, 0.5)
  y <- c(0, 0.5, 1)
  expect_equal(AUCTrapezoid(x,y), 0.75)
  
  # Test for a dip in Y.
  x <- c(0, 0.5, 1)
  y <- c(1, 0.5, 1)
  expect_equal(AUCTrapezoid(x,y), 0.75)
  
  # Test for different minimum cutoffs.
  x <- c(0.5, 0.75, 1)
  y <- c(0.5, 0.75, 1)
  expect_equal(AUCTrapezoid(x,y), 0.5)
  
  # Test for all zeros.
  x <- c(0,0,0)
  y <- c(0,0,0)
  expect_equal(AUCTrapezoid(x,y), NA)
  
  # Test for all ones.
  x <- c(1,1,1)
  y <- c(1,1,1)
  expect_equal(AUCTrapezoid(x,y), NA)
  
  # Test for no range in Y.
  x <- c(0,0.5,1)
  y <- c(1,1,1)
  expect_equal(AUCTrapezoid(x,y), 1)
  
  # Test for no range in X, and X is the max.
  x <- c(1,1,1)
  y <- c(0,0.5,1)
  expect_equal(AUCTrapezoid(x,y), 0)
  
  # Test for no range in X, and X is in the middle.
  x <- c(0.5, 0.5, 0.5)
  y <- c(0, 0.5, 1)
  expect_equal(AUCTrapezoid(x,y), 0.5)
  
  # Test when the X axis is out-of-order.
  x <- c(1, 0.5, 0)
  y <- c(1, 0.5, 0)
  expect_equal(AUCTrapezoid(x,y), 0.5)
  
  # Test for error with negative values.
  x <- c(-1,2,3)
  y <- c(-1,2,3)
  expect_error(AUCTrapezoid(x,y), "ERROR: You have input negative similarities.")
})
test_that("[FERRET] PlotROC() function handles errors appropriately",{
  
  # Check for the appropriate input.
  expect_error(PlotROC("sim", 0.5, 
                       xlab = "Outgroup", ylab = "Ingroup", main = "Plot"),
               "The averageSims parameter must be of type FERRET_Similarities.")
  
  # Test for error with negative values.
  x <- c(-1,2,3)
  y <- c(-1,2,3)
  sim <- methods::new("FERRET_Similarities", ingroup = y, outgroup = x)
  expect_error(PlotROC(sim, 0.5, xlab = "Outgroup", ylab = "Ingroup", main = "Plot"),
               "ERROR: You have input negative similarities.")
})
test_that("[FERRET] ObtainNetworkCutoffs() function yields expected results",{
  # Error should be thrown when input is not in the correct format.
  expect_error(ObtainNetworkCutoffs(NULL, 1), "Results must be of type FERRET_Results.")
  
  # Error should be thrown whenever there are both negatives and positives.
  network1 <- data.frame(source = c("tfA", "tfB", "tfC"),
                         target = c("geneA", "geneB", "geneC"),
                         score = c(1,1,1))
  network2 <- data.frame(source = c("tfA", "tfB", "tfC"),
                         target = c("geneA", "geneB", "geneC"),
                         score = c(1,-1,1))
  results <- methods::new("FERRET_Results", results = list(network1, network2))
  expect_error(ObtainNetworkCutoffs(results, 1), "Negative scores must be adjusted before calculating cutoffs.")
  
  # Check that number of cutoffs is formatted appropriately.
  network2 <- data.frame(source = c("tfA", "tfB", "tfC"),
                         target = c("geneA", "geneB", "geneC"),
                         score = c(1,1,1))
  results <- methods::new("FERRET_Results", results = list(network1, network2))
  expect_error(ObtainNetworkCutoffs(results, 0), "The number of cutoffs must be an integer and must be above 1. You entered: 0")
  expect_error(ObtainNetworkCutoffs(results, -1), "The number of cutoffs must be an integer and must be above 1. You entered: -1")
  expect_error(ObtainNetworkCutoffs(results, 0.5), "The number of cutoffs must be an integer and must be above 1. You entered: 0.5")
  expect_error(ObtainNetworkCutoffs(results, "WOOHOO"), "The number of cutoffs must be an integer and must be above 1. You entered: WOOHOO")
  
  # If all values are the same, only return one cutoff.
  expect_equal(ObtainNetworkCutoffs(results, 10), 1)
  
  # Test for a single network.
  network1 <- data.frame(source = c("tfA", "tfB", "tfC"),
                         target = c("geneA", "geneB", "geneC"),
                         score = c(1,0,2))
  results <- methods::new("FERRET_Results", results = list(network1))
  expect_equal(ObtainNetworkCutoffs(results, 10), seq(0,2,0.2))
  
  # Test for multiple networks but where the min and max are in the same network.
  network1 <- data.frame(source = c("tfA", "tfB", "tfC"),
                         target = c("geneA", "geneB", "geneC"),
                         score = c(1,0,2))
  network2 <- data.frame(source = c("tfA", "tfB", "tfC"),
                         target = c("geneA", "geneB", "geneC"),
                         score = c(0.5,0.6,0.7))
  results <- methods::new("FERRET_Results", results = list(network1, network2))
  expect_equal(ObtainNetworkCutoffs(results, 10), seq(0,2,0.2))
  
  # Test for multiple networks where the min is in one and the max in the other.
  network1 <- data.frame(source = c("tfA", "tfB", "tfC"),
                         target = c("geneA", "geneB", "geneC"),
                         score = c(1,0,0.7))
  network2 <- data.frame(source = c("tfA", "tfB", "tfC"),
                         target = c("geneA", "geneB", "geneC"),
                         score = c(0.5,0.6,2))
  results <- methods::new("FERRET_Results", results = list(network1, network2))
  expect_equal(ObtainNetworkCutoffs(results, 10), seq(0,2,0.2))
  
  # Check for one cutoff.
  results <- methods::new("FERRET_Results", results = list(network1, network2))
  expect_equal(ObtainNetworkCutoffs(results, 1), seq(0,2,2))
})
test_that("[FERRET] ScaleNetworksByPercentile() function yields expected results",{
  
  # Set up networks.
  network1 <- data.frame(source = c("tfA", "tfB", "tfC"),
                         target = c("geneA", "geneB", "geneC"),
                         score = c(1,0,2))
  network2 <- data.frame(source = c("tfA", "tfB", "tfC"),
                         target = c("geneA", "geneB", "geneC"),
                         score = c(1,1,1))
  network3 <- data.frame(source = c("tfA", "tfB", "tfC"),
                         target = c("geneA", "geneB", "geneC"),
                         score = c(0,1,1))
  
  # Check that the network input is formatted appropriately.
  expect_error(ScaleNetworksByPercentile(network1, 3), 
               "networks parameter must be a list of networks with scores in the 3rd column!")
  expect_error(ScaleNetworksByPercentile("WOOHOO", 3), 
               "networks parameter must be a list of networks with scores in the 3rd column!")
  expect_error(ScaleNetworksByPercentile(data.frame(source = c("tfA", "tfB", "tfC"),
                                                    target = c("geneA", "geneB", "geneC"),
                                                    score = c("one", "ttwo", "three")), 3), 
               "networks parameter must be a list of networks with scores in the 3rd column!")
  
  # Check that number of cutoffs is formatted appropriately.
  expect_error(ScaleNetworksByPercentile(list(network1, network2), 0), 
               "The number of cutoffs must be an integer and must be above 1. You entered: 0")
  expect_error(ScaleNetworksByPercentile(list(network1, network2), -1), 
               "The number of cutoffs must be an integer and must be above 1. You entered: -1")
  expect_error(ScaleNetworksByPercentile(list(network1, network2), 0.5), 
               "The number of cutoffs must be an integer and must be above 1. You entered: 0.5")
  expect_error(ScaleNetworksByPercentile(list(network1, network2), "WOOHOO"), 
               "The number of cutoffs must be an integer and must be above 1. You entered: WOOHOO")
  expect_error(ScaleNetworksByPercentile(list(network1, network2, network3), 3), 
               "The number of unique positive values in each network must be greater than the number of cutoffs.")
  expect_error(ScaleNetworksByPercentile(list(network1, network3), 2), 
               "The number of unique positive values in each network must be greater than the number of cutoffs.")
  
  
  # If there is only one cutoff, the cutoff should be at 50%.
  scaledNetworks <- ScaleNetworksByPercentile(list(network1, network3), 1)
  expect_equal(scaledNetworks[[1]][,3], c(0.5, 0, 0.5))
  expect_equal(scaledNetworks[[2]][,3], c(0, 0.5, 0.5))
  
  # If there are two cutoffs, they should be at approximately 0.33 and 0.67.
  scaledNetworks <- ScaleNetworksByPercentile(list(network1), 2)
  expect_equal(scaledNetworks[[1]][,3], c(0.33, 0, 0.67), tolerance = 0.005)
})
test_that("[FERRET] GetEdgesAboveCutoff() function yields expected results",{
  # Test that the correct values are returned.
  network1 <- data.frame(source = c("tfA", "tfB", "tfC"),
                         target = c("geneA", "geneB", "geneC"),
                         score = c(1,0,0.7))
  expect_equal(GetEdgesAboveCutoff(network1, 0), network1)
  expect_equal(GetEdgesAboveCutoff(network1, -1), network1)
  expect_equal(GetEdgesAboveCutoff(network1, 0.5), network1[-2,])
  expect_equal(GetEdgesAboveCutoff(network1, 0.8), network1[1,])
  expect_equal(GetEdgesAboveCutoff(network1, 1), network1[1,])
  
  # Test that an empty data frame is returned where appropriate.
  emptyDf <- as.data.frame(matrix(nrow = 0,ncol = 3))
  colnames(emptyDf) <- colnames(network1)
  emptyDf$source <- as.character(emptyDf$source)
  emptyDf$target <- as.character(emptyDf$target)
  emptyDf$score <- as.numeric(emptyDf$score)
  expect_equal(GetEdgesAboveCutoff(network1, 2), emptyDf)
})
test_that("[FERRET] GetNegative() function yields expected results",{
  # Check that the function returns all negatives for a single network.
  network1 <- data.frame(source = c("tfA", "tfB", "tfC"),
                         target = c("geneA", "geneB", "geneC"),
                         score = c(1,-1,1))
  rownames(network1) <- network1$source
  negOnly <- data.frame(source = c("tfB"),
             target = c("geneB"),
             score = c(1))
  rownames(negOnly) <- c("tfB")
  res <- list(net = network1)
  FERRET_ResultsObj <- methods::new("FERRET_Results",results = res, directory = "tmp")
  expect_equal(GetNegative(FERRET_ResultsObj)@results[[1]], negOnly)
  
  # Check that, if a network has no negatives, nothing is returned.
  network1 <- data.frame(source = c("tfA", "tfB", "tfC"),
                         target = c("geneA", "geneB", "geneC"),
                         score = c(1,1,1))
  rownames(network1) <- network1$source
  emptyDf <- as.data.frame(matrix(nrow = 0,ncol = 3))
  colnames(emptyDf) <- colnames(network1)
  emptyDf$source <- as.character(emptyDf$source)
  emptyDf$target <- as.character(emptyDf$target)
  emptyDf$score <- as.numeric(emptyDf$score)
  rownames(emptyDf) <- as.character(rownames(emptyDf))
  res <- list(net = network1)
  FERRET_ResultsObj <- methods::new("FERRET_Results",results = res, directory = "tmp")
  expect_equal(GetNegative(FERRET_ResultsObj)@results[[1]], emptyDf)
  
  # Check that, if a network has only negatives, the whole network is returned.
  network1_allNeg <- data.frame(source = c("tfA", "tfB", "tfC"),
                         target = c("geneA", "geneB", "geneC"),
                         score = c(-1,-1,-1))
  rownames(network1_allNeg) <- network1_allNeg$source
  res <- list(net = network1_allNeg)
  FERRET_ResultsObj <- methods::new("FERRET_Results",results = res, directory = "tmp")
  expect_equal(GetNegative(FERRET_ResultsObj)@results[[1]], network1)
})
test_that("[FERRET] GetPositive() function yields expected results",{
  # Check that the function returns all negatives for a single network.
  network1 <- data.frame(source = c("tfA", "tfB", "tfC"),
                         target = c("geneA", "geneB", "geneC"),
                         score = c(1,-1,1))
  rownames(network1) <- network1$source
  posOnly <- data.frame(source = c("tfA", "tfC"),
                        target = c("geneA", "geneC"),
                        score = c(1))
  rownames(posOnly) <- c("tfA", "tfC")
  res <- list(net = network1)
  FERRET_ResultsObj <- methods::new("FERRET_Results",results = res, directory = "tmp")
  expect_equal(GetPositive(FERRET_ResultsObj)@results[[1]], posOnly)
  
  # Check that, if a network has no negatives, nothing is returned.
  network1 <- data.frame(source = c("tfA", "tfB", "tfC"),
                         target = c("geneA", "geneB", "geneC"),
                         score = c(-1,-1,-1))
  rownames(network1) <- network1$source
  emptyDf <- as.data.frame(matrix(nrow = 0,ncol = 3))
  colnames(emptyDf) <- colnames(network1)
  emptyDf$source <- as.character(emptyDf$source)
  emptyDf$target <- as.character(emptyDf$target)
  emptyDf$score <- as.numeric(emptyDf$score)
  rownames(emptyDf) <- as.character(rownames(emptyDf))
  res <- list(net = network1)
  FERRET_ResultsObj <- methods::new("FERRET_Results",results = res, directory = "tmp")
  expect_equal(GetPositive(FERRET_ResultsObj)@results[[1]], emptyDf)
  
  # Check that, if a network has only negatives, the whole network is returned.
  network1_allPos <- data.frame(source = c("tfA", "tfB", "tfC"),
                                target = c("geneA", "geneB", "geneC"),
                                score = c(1,1,1))
  rownames(network1_allPos) <- network1_allPos$source
  res <- list(net = network1_allPos)
  FERRET_ResultsObj <- methods::new("FERRET_Results",results = res, directory = "tmp")
  expect_equal(GetPositive(FERRET_ResultsObj)@results[[1]], network1_allPos)
})
test_that("[FERRET] JaccardSim() function yields expected results",{
  
  # If one network is empty, overlap should be 0.
  network1 <- data.frame(source = c("tfA", "tfB", "tfC"),
                         target = c("geneA", "geneB", "geneC"),
                         score = c(0.2,0.4,0.6))
  rownames(network1) <- paste(network1$source, network1$target, sep = "_")
  emptyDf <- as.data.frame(matrix(nrow = 0,ncol = 3))
  colnames(emptyDf) <- colnames(network1)
  emptyDf$source <- as.character(emptyDf$source)
  emptyDf$target <- as.character(emptyDf$target)
  emptyDf$score <- as.numeric(emptyDf$score)
  rownames(emptyDf) <- as.character(rownames(emptyDf))
  expect_equal(JaccardSim(network1, emptyDf), 0)
  
  # If both networks are empty, overlap should be NA.
  expect_equal(JaccardSim(emptyDf, emptyDf), NA)
  
  # If one network's weights are half of the others but they are otherwise
  # the same, the overlap should be 0.5.
  network2 <- data.frame(source = c("tfA", "tfB", "tfC"),
                         target = c("geneA", "geneB", "geneC"),
                         score = c(0.1,0.2,0.3))
  rownames(network2) <- paste(network2$source, network2$target, sep = "_")
  expect_equal(JaccardSim(network1, network2), 0.5)
  
  # If one network contains exactly half of the other network's edges,
  # the score should be 0.5.
  network1 <- data.frame(source = c("tfA", "tfB"),
                         target = c("geneA", "geneB"),
                         score = c(0.4,0.4))
  rownames(network1) <- paste(network1$source, network1$target, sep = "_")
  network2 <- data.frame(source = c("tfA"),
                         target = c("geneA"),
                         score = c(0.4))
  rownames(network2) <- paste(network2$source, network2$target, sep = "_")
  expect_equal(JaccardSim(network1, network2), 0.5)
  
  # If exactly half of each network overlaps with the other, the score should be 0.33.
  network1 <- data.frame(source = c("tfA", "tfB"),
                         target = c("geneA", "geneB"),
                         score = c(0.4,0.4))
  rownames(network1) <- paste(network1$source, network1$target, sep = "_")
  network2 <- data.frame(source = c("tfA", "tfC"),
                         target = c("geneA", "geneC"),
                         score = c(0.4, 0.4))
  rownames(network2) <- paste(network2$source, network2$target, sep = "_")
  expect_equal(JaccardSim(network1, network2), 1/3)
})
test_that("[FERRET] DegreeSim() function yields expected results",{
  
  # Validate calculations when one network is empty.
  network1 <- data.frame(source = c("tf1", "tf2", "tf3", "tf1", "tf2", "tf3"),
                         target = c("geneA", "geneA", "geneA", "geneB", "geneB", "geneB"),
                         score = c(1,1,1,1,1,1))
  emptyDf <- as.data.frame(matrix(nrow = 0,ncol = 3))
  colnames(emptyDf) <- colnames(network1)
  emptyDf$source <- as.character(emptyDf$source)
  emptyDf$target <- as.character(emptyDf$target)
  emptyDf$score <- as.numeric(emptyDf$score)
  rownames(emptyDf) <- as.character(rownames(emptyDf))
  expect_equal(DegreeSim(network1, emptyDf), 0)
  
  # If both networks are empty, overlap should be NA.
  expect_equal(DegreeSim(emptyDf, emptyDf), NA)
  
  # Validate calculations on a simple example. Here, the degree distributions are
  # as follows:
  # out-degree:
  #     net1  net2  net3  net4
  # tf1 1/3   2/7   2/7   2/7
  # tf2 1/3   2/7   2/7   2/7
  # tf3 1/3   2/7   2/7   3/7
  # tf4 0     1/7   1/7   0
  # in-degree:
  #       net1  net2  net3  net4
  # geneA 1/2   3/7   3/7   3/7
  # geneB 1/2   3/7   4/7   3/7
  # geneC 0     1/7   0     1/7 
  network1 <- data.frame(source = c("tf1", "tf2", "tf3", "tf1", "tf2", "tf3"),
                         target = c("geneA", "geneA", "geneA", "geneB", "geneB", "geneB"),
                         score = c(1,1,1,1,1,1))
  network2 <- data.frame(source = c("tf1", "tf2", "tf3", "tf1", "tf2", "tf3", "tf4"),
                         target = c("geneA", "geneA", "geneA", "geneB", "geneB", "geneB", "geneC"),
                         score = c(1,1,1,1,1,1,1))
  network3 <- data.frame(source = c("tf1", "tf2", "tf3", "tf1", "tf2", "tf3", "tf4"),
                         target = c("geneA", "geneA", "geneA", "geneB", "geneB", "geneB", "geneB"),
                         score = c(1,1,1,1,1,1,1))
  network4 <- data.frame(source = c("tf1", "tf2", "tf3", "tf1", "tf2", "tf3", "tf3"),
                         target = c("geneA", "geneA", "geneA", "geneB", "geneB", "geneB", "geneC"),
                         score = c(1,1,1,1,1,1,1))
  outDiff12 <- (((1/3) - (2/7)) * 3 + (1/7)) / 2
  outDiff13 <- outDiff12
  outDiff14 <- (((1/3) - (2/7)) * 2 + (3/7) - (1/3)) / 2
  outDiff23 <- 0
  outDiff24 <- (2/7) / 2
  outDiff34 <- (2/7) / 2
  inDiff12 <- (((1/2) - (3/7)) * 2 + (1/7)) / 2
  inDiff13 <- ((1/2) - (3/7) + (4/7) - (1/2)) / 2
  inDiff14 <- inDiff12
  inDiff23 <- (2/7) / 2
  inDiff24 <- 0
  inDiff34 <- (2/7) / 2
  expect_equal(DegreeSim(network1, network1), 1)
  expect_equal(DegreeSim(network1, network2), (2 * (1 - outDiff12) * (1 - inDiff12)) /  ((1 - outDiff12) + (1 - inDiff12)))
  expect_equal(DegreeSim(network2, network1), (2 * (1 - outDiff12) * (1 - inDiff12)) /  ((1 - outDiff12) + (1 - inDiff12)))
  expect_equal(DegreeSim(network1, network3), (2 * (1 - outDiff13) * (1 - inDiff13)) /  ((1 - outDiff13) + (1 - inDiff13)))
  expect_equal(DegreeSim(network1, network4), (2 * (1 - outDiff14) * (1 - inDiff14)) /  ((1 - outDiff14) + (1 - inDiff14)))
  expect_equal(DegreeSim(network2, network3), (2 * (1 - outDiff23) * (1 - inDiff23)) /  ((1 - outDiff23) + (1 - inDiff23)))
  expect_equal(DegreeSim(network2, network4), (2 * (1 - outDiff24) * (1 - inDiff24)) /  ((1 - outDiff24) + (1 - inDiff24)))
  expect_equal(DegreeSim(network3, network4), (2 * (1 - outDiff34) * (1 - inDiff34)) /  ((1 - outDiff34) + (1 - inDiff34)))
  
  # Check that when one network is identical to another but weights are on a different
  # scale, value doesn't change.
  network1Half <- network1
  network1Half$score <- network1$score / 2
  expect_equal(DegreeSim(network1, network2), DegreeSim(network1Half, network2))
  
  # Check that when one weight is different, value changes. We should now have:
  # out-degree:
  #     net1  net2Changed
  # tf1 1/3   4/13 
  # tf2 1/3   4/13
  # tf3 1/3   4/13
  # tf4 0     1/13 
  # in-degree:
  #       net1  net2 
  # geneA 1/2   6/13 
  # geneB 1/2   6/13 
  # geneC 0     1/13 
  network2Changed <- network2
  network2Changed[7,"score"] <- 0.5
  outDiff12Changed <- (((1/3) - (4/13)) * 3 + (1/13)) / 2
  inDiff12Changed <- (((1/2) - (6/13)) * 2 + (1/13)) / 2
  expect_equal(DegreeSim(network1, network2Changed), (2 * (1 - outDiff12Changed) * (1 - inDiff12Changed)) /  ((1 - outDiff12Changed) + (1 - inDiff12Changed)))
  
  # Check that adding a self-edge causes expected behavior. We should now have:
  # out-degree:
  #     net1  net1Changed
  # tf1 1/3   3/7 
  # tf2 1/3   2/7
  # tf3 1/3   2/7
  # in-degree:
  #       net1  net1Changed 
  # tf1   0     1/7
  # geneA 1/2   3/7 
  # geneB 1/2   3/7 
  network1SelfEdge <- rbind(network1, data.frame(source = "tf1", target = "tf1", score = 1))
  outDiff11SelfEdge <- (((1/3) - (2/7)) * 2 + ((3/7) - (1/3))) / 2
  inDiff11SelfEdge <- ((1/7) + ((1/2) - (3/7)) * 2) / 2
  expect_equal(DegreeSim(network1, network1SelfEdge), (2 * (1 - outDiff11SelfEdge) * (1 - inDiff11SelfEdge)) /  ((1 - outDiff11SelfEdge) + (1 - inDiff11SelfEdge)))
  
  # Check that adding a "backwards" edge causes expected behavior. We should now have:
  # out-degree:
  #       net1  net1Backwards
  # tf1   1/3   2/7 
  # tf2   1/3   2/7
  # tf3   1/3   2/7
  # geneB 0     1/7
  # in-degree:
  #       net1  net1Backwards 
  # tf1   0     1/7
  # geneA 1/2   3/7 
  # geneB 1/2   3/7 
  network1Backwards <- rbind(network1, data.frame(source = "geneB", target = "tf1", score = 1))
  outDiff11Backwards <- (((1/3) - (2/7)) * 3 + (1/7)) / 2
  inDiff11Backwards <- ((1/7) + ((1/2) - (3/7)) * 2) / 2
  expect_equal(DegreeSim(network1, network1Backwards), (2 * (1 - outDiff11Backwards) * (1 - inDiff11Backwards)) /  ((1 - outDiff11Backwards) + (1 - inDiff11Backwards)))
})
test_that("[FERRET] ModularitySim() function yields expected results",{
  # Network to test
  net1 <- data.frame(source = c(rep("tf1", 4), rep("tf2", 4), rep("tf3", 3), rep("tf4", 3), "tf3"),
                     target = c(rep(c("gene1", "gene2", "gene3", "gene4"), 2),
                                rep(c("gene5", "gene6", "gene7"), 2), "gene3"),
                     score = rep(1, 15))
  net2 <- data.frame(source = c(rep("tf1", 4), rep("tf2", 4), rep("tf3", 3), rep("tf4", 3), "tf4"),
                     target = c(rep(c("gene1", "gene2", "gene3", "gene4"), 2),
                                rep(c("gene5", "gene6", "gene7"), 2), "gene3"),
                     score = rep(1, 15))
  net3 <- data.frame(source = c("tf1", "tf2", rep("tf3", 3), rep("tf4", 5)),
                     target = c("gene1", "gene2", rep(c("gene3", "gene4", "gene5"), 2), "gene6", "gene7"),
                     score = rep(1, 10))
  
  # Check that number of dimensions is correct for different values of k.
  expect_gt(ModularitySim(net1, net2), ModularitySim(net1, net3))
  
  # Check that ModularitySim to an empty network is 0.
  emptyDf <- as.data.frame(matrix(nrow = 0,ncol = 3))
  colnames(emptyDf) <- colnames(net1)
  emptyDf$source <- as.character(emptyDf$source)
  emptyDf$target <- as.character(emptyDf$target)
  emptyDf$score <- as.numeric(emptyDf$score)
  rownames(emptyDf) <- as.character(rownames(emptyDf))
  expect_equal(ModularitySim(net1, emptyDf), 0)
  
  # Check that similarity between an empty network and itself is NA.
  expect_equal(ModularitySim(emptyDf, emptyDf), NA)
  
  # Check that ModularitySim between a network and itself is high with a backward edge.
  net1Backwards <- rbind(net1, data.frame(source = "gene5", target = "tf1", score = 1))
  expect_gt(ModularitySim(net1, net1Backwards), ModularitySim(net1, net3))
  
  # Check that ModularitySim between a network and itself is high with a self-loop.
  net1SelfLoop <- rbind(net1, data.frame(source = "tf1", target = "tf1", score = 1))
  expect_gt(ModularitySim(net1, net1SelfLoop), ModularitySim(net1, net3))
  
  # Check that ModularitySim between a network and itself is 1.
  expect_equal(ModularitySim(net1, net1), 1)
})
test_that("[FERRET] ComputeProjection() function yields expected results",{
  # Check networks with similar clustering vs. dissimilar clustering
  # and make sure that their similarities are as expected.
  net1 <- data.frame(source = c(rep("tf1", 4), rep("tf2", 4), rep("tf3", 3), rep("tf4", 3), "tf3"),
                     target = c(rep(c("gene1", "gene2", "gene3", "gene4"), 2),
                                rep(c("gene5", "gene6", "gene7"), 2), "gene3"),
                     score = rep(1, 15))
    
  # Check invalid values of k.
  expect_error(ComputeProjection(net1, 0, c(net1$source, net1$target)),
               "Number of eigenvectors k must be an integer greater than 0!")
  expect_error(ComputeProjection(net1, -1, c(net1$source, net1$target)),
               "Number of eigenvectors k must be an integer greater than 0!")
  expect_error(ComputeProjection(net1, 0.6, c(net1$source, net1$target)),
               "Number of eigenvectors k must be an integer greater than 0!")
  expect_error(ComputeProjection(net1, "hello", c(net1$source, net1$target)),
               "Number of eigenvectors k must be an integer greater than 0!")
  expect_error(ComputeProjection(net1, 300, c(net1$source, net1$target)),
               "Number of eigenvectors must be less than the dimensionality of the adjacency matrix!")
  expect_error(ComputeProjection(net1, 11, c(net1$source, net1$target)),
               "Number of eigenvectors must be less than the dimensionality of the adjacency matrix!")
  
  # Check that projection has expected number of dimensions for different values of k.
  expect_equal(dim(ComputeProjection(net1, 1, c(net1$source, net1$target))),
               c(length(unique(c(net1$source, net1$target))), 1))
  expect_equal(dim(ComputeProjection(net1, 3, c(net1$source, net1$target))),
               c(length(unique(c(net1$source, net1$target))), 3))
  expect_equal(dim(ComputeProjection(net1, 5, c(net1$source, net1$target))),
               c(length(unique(c(net1$source, net1$target))), 5))
  
  # Make sure it still works when we add additional nodes.
  suppressWarnings({
    extendedNodes <- c(c(net1$source, net1$target, "something",
                         "somethingElse", "yetAnotherThing"))
    expect_equal(dim(ComputeProjection(net1, 1, extendedNodes)),
                 c(length(unique(extendedNodes)), 1))
    expect_equal(dim(ComputeProjection(net1, 5, extendedNodes)),
                 c(length(unique(extendedNodes)), 5))
    expect_equal(dim(ComputeProjection(net1, 11, extendedNodes)),
                 c(length(unique(extendedNodes)), 11))
    expect_equal(dim(ComputeProjection(net1, 12, extendedNodes)),
                 c(length(unique(extendedNodes)), 12))
    expect_error(ComputeProjection(net1, 14, extendedNodes),
                 "Number of eigenvectors must be less than the dimensionality of the adjacency matrix!")
  })
  
  
  # Check performance with an empty matrix.
  suppressWarnings({
    emptyDf <- as.data.frame(matrix(nrow = 0,ncol = 3))
    colnames(emptyDf) <- colnames(net1)
    emptyDf$source <- as.character(emptyDf$source)
    emptyDf$target <- as.character(emptyDf$target)
    emptyDf$score <- as.numeric(emptyDf$score)
    rownames(emptyDf) <- as.character(rownames(emptyDf))
    expect_error(ComputeProjection(emptyDf, 2, c()),
                 "Number of eigenvectors must be less than the dimensionality of the adjacency matrix!")
    extraNodes <- unique(extendedNodes)
    projectionZero <- matrix(data = rep(0, length(extraNodes) * 2),
                             nrow = length(extraNodes))
    colnames(projectionZero) <- c("PC1", "PC2")
    expect_equal(ComputeProjection(emptyDf, 2, extraNodes), projectionZero)
  })
  
})
test_that("[FERRET] SubspaceSim() function yields expected results",{
  # Check the same networks as we checked for ModularitySim().
  net1 <- data.frame(source = c(rep("tf1", 4), rep("tf2", 4), rep("tf3", 3), rep("tf4", 3), "tf3"),
                     target = c(rep(c("gene1", "gene2", "gene3", "gene4"), 2),
                                rep(c("gene5", "gene6", "gene7"), 2), "gene3"),
                     score = rep(1, 15))
  net2 <- data.frame(source = c(rep("tf1", 4), rep("tf2", 4), rep("tf3", 3), rep("tf4", 3), "tf4"),
                     target = c(rep(c("gene1", "gene2", "gene3", "gene4"), 2),
                                rep(c("gene5", "gene6", "gene7"), 2), "gene3"),
                     score = rep(1, 15))
  net3 <- data.frame(source = c("tf1", "tf2", rep("tf3", 3), rep("tf4", 5)),
                     target = c("gene1", "gene2", rep(c("gene3", "gene4", "gene5"), 2), "gene6", "gene7"),
                     score = rep(1, 10))
  allNodes <- unique(c(net1[,"source"], net1[,"target"], net2[,"source"], net2[,"source"],
                net3[,"source"], net3[,"target"]))
  projection1 <- ComputeProjection(net1, 5, allNodes)
  projection2 <- ComputeProjection(net2, 5, allNodes)
  projection3 <- ComputeProjection(net3, 5, allNodes)
  expect_gt(SubspaceSim(projection1, projection2, 5), 
            SubspaceSim(projection1, projection3, 5))
  
  # Check that SubspaceSim() to the zero matrix is 0.
  emptyDf <- as.data.frame(matrix(nrow = 0,ncol = 3))
  colnames(emptyDf) <- colnames(net1)
  emptyDf$source <- as.character(emptyDf$source)
  emptyDf$target <- as.character(emptyDf$target)
  emptyDf$score <- as.numeric(emptyDf$score)
  rownames(emptyDf) <- as.character(rownames(emptyDf))
  projectionZero <- ComputeProjection(emptyDf, 2, allNodes)
  expect_equal(SubspaceSim(projection1, projectionZero, 5), 0)
  expect_equal(SubspaceSim(projection2, projectionZero,  5), 0)
  expect_equal(SubspaceSim(projection3, projectionZero, 5), 0)
  
  # Check that SubspaceSim between a network and itself is 1.
  expect_equal(SubspaceSim(projection1, projection1, 5), 1)
  expect_equal(SubspaceSim(projection3, projection3, 5), 1)
  
  # Check that SubspaceSim between an empty projection and itself is NA.
  expect_equal(SubspaceSim(projectionZero, projectionZero, 5), NA)
  
  # Check that SubspaceSim between a network and itself is high with a self-loop.
  net1SelfLoop <- rbind(net1, data.frame(source = "tf1", target = "tf1", score = 1))
  projection1SelfLoop <- ComputeProjection(net1SelfLoop, 5, allNodes)
  expect_gt(SubspaceSim(projection1, projection1SelfLoop, 5), 
                     SubspaceSim(projection1, projection3, 5))
})
test_that("[FERRET] ComputeSimilarityForGroup() function yields expected results",{
  # Check that the average similarity for the group is output as expected for each
  # type of similarity. Net1 should be similar to net2 and net4 but dissimilar to
  # net3 and net5.
  net1 <- data.frame(source = c(rep("tf1", 4), rep("tf2", 4), rep("tf3", 3), rep("tf4", 3), "tf3"),
                     target = c(rep(c("gene1", "gene2", "gene3", "gene4"), 2),
                                rep(c("gene5", "gene6", "gene7"), 2), "gene3"),
                     score = seq(0.1, 1.5, 0.1))
  rownames(net1) <- paste(net1$source, net1$target, sep = "_")
  net2 <- data.frame(source = c(rep("tf1", 4), rep("tf2", 4), rep("tf3", 3), rep("tf4", 3), "tf4"),
                     target = c(rep(c("gene1", "gene2", "gene3", "gene4"), 2),
                                rep(c("gene5", "gene6", "gene7"), 2), "gene3"),
                     score = seq(0.1, 1.5, 0.1))
  rownames(net2) <- paste(net2$source, net2$target, sep = "_")
  net3 <- data.frame(source = c("tf1", "tf2", rep("tf3", 3), rep("tf4", 5)),
                     target = c("gene1", "gene2", rep(c("gene3", "gene4", "gene5"), 2), "gene6", "gene7"),
                     score = seq(0.1, 1.0, 0.1))
  rownames(net3) <- paste(net3$source, net3$target, sep = "_")
  net4 <- data.frame(source = c(rep("tf1", 4), rep("tf2", 4), rep("tf3", 3), rep("tf4", 3), "tf4"),
                     target = c(rep(c("gene1", "gene2", "gene3", "gene4"), 2),
                                rep(c("gene5", "gene6", "gene7"), 2), "gene4"),
                     score = seq(0.1, 1.5, 0.1))
  rownames(net4) <- paste(net4$source, net4$target, sep = "_")
  net5 <- data.frame(source = c("tf1", "tf2", rep("tf3", 3), rep("tf4", 5)),
                     target = c("gene2", "gene1", rep(c("gene3", "gene4", "gene5"), 2), "gene6", "gene7"),
                     score = seq(0.1, 1.0, 0.1))
  rownames(net5) <- paste(net5$source, net5$target, sep = "_")
  k = 5
  allNodes <- c(net1$source, net1$target, net2$source, net2$target, net3$target,
                net4$source, net4$target, net5$source, net5$target)
  identityMat <- data.frame(source = allNodes, target = allNodes, score = rep(1, length(allNodes)))
  projection1 <- ComputeProjection(net1, k, allNodes)
  projection2 <- ComputeProjection(net2, k, allNodes)
  projection3 <- ComputeProjection(net3, k, allNodes)
  projection4 <- ComputeProjection(net4, k, allNodes)
  projection5 <- ComputeProjection(net5, k, allNodes)
  projectionI <- ComputeProjection(identityMat, k, allNodes)
  resultList <- list(net1,net2,net3,net4,net5)
  names(resultList) <- c("net1", "net2", "net3", "net4", "net5")
  projectionList <- list(projection1, projection2, projection3, projection4, projection5)
  names(projectionList) <- c("net1", "net2", "net3", "net4", "net5")
  FERRET_ResultsObj <- methods::new("FERRET_Results",results = resultList, directory = "tmp")
  FERRET_Comparisons <- methods::new("FERRET_Comparisons", ingroup = c(new("FERRET_Comparison", source = "net1", target = "net2"),
                                                                       new("FERRET_Comparison", source = "net1", target = "net4")),
                                     outgroup = c(new("FERRET_Comparison", source = "net1", target = "net3"),
                                                  new("FERRET_Comparison", source = "net1", target = "net5")))
  ingroupSim <- ComputeSimilarityForGroup(results = FERRET_ResultsObj, comparisons = FERRET_Comparisons@ingroup,
                            metric = "jaccard", k = k)
  outgroupSim <- ComputeSimilarityForGroup(results = FERRET_ResultsObj, comparisons = FERRET_Comparisons@outgroup,
                                          metric = "jaccard", k = k)
  expect_gt(ingroupSim, outgroupSim)
  ingroupSim <- ComputeSimilarityForGroup(results = FERRET_ResultsObj, comparisons = FERRET_Comparisons@ingroup,
                                          metric = "degree", k = k)
  outgroupSim <- ComputeSimilarityForGroup(results = FERRET_ResultsObj, comparisons = FERRET_Comparisons@outgroup,
                                           metric = "degree", k = k)
  expect_gt(ingroupSim, outgroupSim)
  ingroupSim <- ComputeSimilarityForGroup(results = FERRET_ResultsObj, comparisons = FERRET_Comparisons@ingroup,
                                          metric = "modularity", k = k)
  outgroupSim <- ComputeSimilarityForGroup(results = FERRET_ResultsObj, comparisons = FERRET_Comparisons@outgroup,
                                           metric = "modularity", k = k)
  expect_gt(ingroupSim, outgroupSim)
  ingroupSim <- ComputeSimilarityForGroup(results = FERRET_ResultsObj, comparisons = FERRET_Comparisons@ingroup,
                                          projections = projectionList, projectionIdentity = projectionI,
                                          metric = "subspace", k = k)
  outgroupSim <- ComputeSimilarityForGroup(results = FERRET_ResultsObj, comparisons = FERRET_Comparisons@outgroup,
                                           projections = projectionList, projectionIdentity = projectionI,
                                           metric = "subspace", k = k)
  expect_gt(ingroupSim, outgroupSim)
  
  # Check when there is an uneven number of comparisons.
  FERRET_Comparisons <- methods::new("FERRET_Comparisons", ingroup = c(new("FERRET_Comparison", source = "net3", target = "net5")),
                                     outgroup = c(new("FERRET_Comparison", source = "net3", target = "net1"),
                                                  new("FERRET_Comparison", source = "net3", target = "net2"),
                                                  new("FERRET_Comparison", source = "net3", target = "net4")))
  ingroupSim <- ComputeSimilarityForGroup(results = FERRET_ResultsObj, comparisons = FERRET_Comparisons@ingroup,
                                          metric = "jaccard", k = k)
  outgroupSim <- ComputeSimilarityForGroup(results = FERRET_ResultsObj, comparisons = FERRET_Comparisons@outgroup,
                                           metric = "jaccard", k = k)
  expect_gt(ingroupSim, outgroupSim)
  ingroupSim <- ComputeSimilarityForGroup(results = FERRET_ResultsObj, comparisons = FERRET_Comparisons@ingroup,
                                          metric = "degree", k = k)
  outgroupSim <- ComputeSimilarityForGroup(results = FERRET_ResultsObj, comparisons = FERRET_Comparisons@outgroup,
                                           metric = "degree", k = k)
  expect_gt(ingroupSim, outgroupSim)
  ingroupSim <- ComputeSimilarityForGroup(results = FERRET_ResultsObj, comparisons = FERRET_Comparisons@ingroup,
                                          metric = "modularity", k = k)
  outgroupSim <- ComputeSimilarityForGroup(results = FERRET_ResultsObj, comparisons = FERRET_Comparisons@outgroup,
                                           metric = "modularity", k = k)
  expect_gt(ingroupSim, outgroupSim)
  ingroupSim <- ComputeSimilarityForGroup(results = FERRET_ResultsObj, comparisons = FERRET_Comparisons@ingroup,
                                          projections = projectionList, projectionIdentity = projectionI,
                                          metric = "subspace", k = k)
  outgroupSim <- ComputeSimilarityForGroup(results = FERRET_ResultsObj, comparisons = FERRET_Comparisons@outgroup,
                                           projections = projectionList, projectionIdentity = projectionI,
                                           metric = "subspace", k = k)
  expect_gt(ingroupSim, outgroupSim)
  
  # Test that an error is thrown if the comparisons and result object lists do not
  # correspond.
  FERRET_Comparisons_invalid <- methods::new("FERRET_Comparisons", ingroup = c(new("FERRET_Comparison", source = "net1", target = "net2"),
                                                                       new("FERRET_Comparison", source = "net1", target = "yoohoo")),
                                     outgroup = c(new("FERRET_Comparison", source = "net1", target = "yoohoo"),
                                                  new("FERRET_Comparison", source = "net1", target = "net5")))
  expect_error(ComputeSimilarityForGroup(results = FERRET_ResultsObj, comparisons = FERRET_Comparisons_invalid@ingroup,
                                         metric = "jaccard", k = k), "Invalid network 'yoohoo' in comparisons object!")
  expect_error(ComputeSimilarityForGroup(results = FERRET_ResultsObj, comparisons = FERRET_Comparisons_invalid@outgroup,
                                         metric = "modularity", k = k), "Invalid network 'yoohoo' in comparisons object!")
  expect_error(ComputeSimilarityForGroup(results = FERRET_ResultsObj, comparisons = FERRET_Comparisons_invalid@ingroup,
                                         metric = "degree", k = k), "Invalid network 'yoohoo' in comparisons object!")
  expect_error(ComputeSimilarityForGroup(results = FERRET_ResultsObj, comparisons = FERRET_Comparisons_invalid@outgroup,
                                         metric = "subspace", k = k), "Invalid network 'yoohoo' in comparisons object!")
  
  # Test that an error is thrown if an invalid similarity metric is passed.
  expect_error(ComputeSimilarityForGroup(results = FERRET_ResultsObj, comparisons = FERRET_Comparisons@outgroup,
                                         metric = "myAmazingMetric", k = k), "Invalid metric to use for similarity: myAmazingMetric")
})
test_that("[FERRET] ComputeSimilarities() function yields expected results",{
  # Test that numeric values are returned for all cutoffs and that the results
  # are of type FERRET_Similarities.
  net1 <- data.frame(source = c(rep("tf1", 4), rep("tf2", 4), rep("tf3", 3), rep("tf4", 3), "tf3"),
                     target = c(rep(c("gene1", "gene2", "gene3", "gene4"), 2),
                                rep(c("gene5", "gene6", "gene7"), 2), "gene3"),
                     score = seq(0.1, 1.5, 0.1))
  rownames(net1) <- paste(net1$source, net1$target, sep = "_")
  net2 <- data.frame(source = c(rep("tf1", 4), rep("tf2", 4), rep("tf3", 3), rep("tf4", 3), "tf4"),
                     target = c(rep(c("gene1", "gene2", "gene3", "gene4"), 2),
                                rep(c("gene5", "gene6", "gene7"), 2), "gene3"),
                     score = seq(0.1, 1.5, 0.1))
  rownames(net2) <- paste(net2$source, net2$target, sep = "_")
  net3 <- data.frame(source = c("tf1", "tf2", rep("tf3", 3), rep("tf4", 5)),
                     target = c("gene1", "gene2", rep(c("gene3", "gene4", "gene5"), 2), "gene6", "gene7"),
                     score = seq(0.1, 1.0, 0.1))
  rownames(net3) <- paste(net3$source, net3$target, sep = "_")
  net4 <- data.frame(source = c(rep("tf1", 4), rep("tf2", 4), rep("tf3", 3), rep("tf4", 3), "tf4"),
                     target = c(rep(c("gene1", "gene2", "gene3", "gene4"), 2),
                                rep(c("gene5", "gene6", "gene7"), 2), "gene4"),
                     score = seq(0.1, 1.5, 0.1))
  rownames(net4) <- paste(net4$source, net4$target, sep = "_")
  net5 <- data.frame(source = c("tf1", "tf2", rep("tf3", 3), rep("tf4", 5)),
                     target = c("gene2", "gene1", rep(c("gene3", "gene4", "gene5"), 2), "gene6", "gene7"),
                     score = seq(0.1, 1.0, 0.1))
  rownames(net5) <- paste(net5$source, net5$target, sep = "_")
  k = 5
  allNodes <- c(net1$source, net1$target, net2$source, net2$target, net3$target,
                net4$source, net4$target, net5$source, net5$target)
  identityMat <- data.frame(source = allNodes, target = allNodes, score = rep(1, length(allNodes)))
  resultList <- list(net1,net2,net3,net4,net5)
  names(resultList) <- c("net1", "net2", "net3", "net4", "net5")
  FERRET_ResultsObj <- methods::new("FERRET_Results",results = resultList, directory = "tmp")
  FERRET_Comparisons <- methods::new("FERRET_Comparisons", ingroup = c(new("FERRET_Comparison", source = "net1", target = "net2"),
                                                                       new("FERRET_Comparison", source = "net1", target = "net4")),
                                     outgroup = c(new("FERRET_Comparison", source = "net1", target = "net3"),
                                                  new("FERRET_Comparison", source = "net1", target = "net5")))
  suppressWarnings({
    similarities <- ComputeSimilarities(results = FERRET_ResultsObj, comparisons = FERRET_Comparisons,
                                        cutoffs = seq(0.1, 1.5, 0.1), metric = "jaccard", k = k)
    expect_true(is(similarities, "FERRET_Similarities"))
    expect_equal(names(similarities@ingroup), as.character(seq(0.1, 1.5, 0.1)))
    expect_equal(names(similarities@outgroup), as.character(seq(0.1, 1.5, 0.1)))
    similarities <- ComputeSimilarities(results = FERRET_ResultsObj, comparisons = FERRET_Comparisons,
                                        cutoffs = seq(0.1, 1.5, 0.1), metric = "degree", k = k)
    expect_true(is(similarities, "FERRET_Similarities"))
    expect_equal(names(similarities@ingroup), as.character(seq(0.1, 1.5, 0.1)))
    expect_equal(names(similarities@outgroup), as.character(seq(0.1, 1.5, 0.1)))
    similarities <- ComputeSimilarities(results = FERRET_ResultsObj, comparisons = FERRET_Comparisons,
                                        cutoffs = seq(0.1, 1.5, 0.1), metric = "modularity", k = k)
    expect_true(is(similarities, "FERRET_Similarities"))
    expect_equal(names(similarities@ingroup), as.character(seq(0.1, 1.5, 0.1)))
    expect_equal(names(similarities@outgroup), as.character(seq(0.1, 1.5, 0.1)))
    similarities <- ComputeSimilarities(results = FERRET_ResultsObj, comparisons = FERRET_Comparisons,
                                        cutoffs = seq(0.1, 1.5, 0.1), metric = "subspace", k = k)
    expect_true(is(similarities, "FERRET_Similarities"))
    expect_equal(names(similarities@ingroup), as.character(seq(0.1, 1.5, 0.1)))
    expect_equal(names(similarities@outgroup), as.character(seq(0.1, 1.5, 0.1)))
    
    # Check that we get NA values if we extend past the range of both network weights.
    similarities <- ComputeSimilarities(results = FERRET_ResultsObj, comparisons = FERRET_Comparisons,
                                        cutoffs = seq(0.1, 2.0, 0.1), metric = "jaccard", k = k)
    expect_true(is(similarities, "FERRET_Similarities"))
    expect_equal(unname(similarities@ingroup[as.character(seq(1.6, 2.0, 0.1))]), as.numeric(rep(NA, 5)))
    expect_equal(unname(similarities@ingroup[as.character(seq(1.6, 2.0, 0.1))]), as.numeric(rep(NA, 5)))
    expect_equal(length(which(is.na(similarities@ingroup[as.character(seq(0.1, 1.5, 0.1))]))), 0)
    expect_equal(length(which(is.na(similarities@outgroup[as.character(seq(0.1, 1.5, 0.1))]))), 0)
    similarities <- ComputeSimilarities(results = FERRET_ResultsObj, comparisons = FERRET_Comparisons,
                                        cutoffs = seq(0.1, 2.0, 0.1), metric = "degree", k = k)
    expect_true(is(similarities, "FERRET_Similarities"))
    expect_equal(unname(similarities@ingroup[as.character(seq(1.6, 2.0, 0.1))]), as.numeric(rep(NA, 5)))
    expect_equal(unname(similarities@ingroup[as.character(seq(1.6, 2.0, 0.1))]), as.numeric(rep(NA, 5)))
    expect_equal(length(which(is.na(similarities@ingroup[as.character(seq(0.1, 1.5, 0.1))]))), 0)
    expect_equal(length(which(is.na(similarities@outgroup[as.character(seq(0.1, 1.5, 0.1))]))), 0)
    similarities <- ComputeSimilarities(results = FERRET_ResultsObj, comparisons = FERRET_Comparisons,
                                        cutoffs = seq(0.1, 2.0, 0.1), metric = "modularity", k = k)
    expect_true(is(similarities, "FERRET_Similarities"))
    expect_equal(unname(similarities@ingroup[as.character(seq(1.6, 2.0, 0.1))]), as.numeric(rep(NA, 5)))
    expect_equal(unname(similarities@ingroup[as.character(seq(1.6, 2.0, 0.1))]), as.numeric(rep(NA, 5)))
    expect_equal(length(which(is.na(similarities@ingroup[as.character(seq(0.1, 1.5, 0.1))]))), 0)
    expect_equal(length(which(is.na(similarities@outgroup[as.character(seq(0.1, 1.5, 0.1))]))), 0)
    similarities <- ComputeSimilarities(results = FERRET_ResultsObj, comparisons = FERRET_Comparisons,
                                      cutoffs = seq(0.1, 2.0, 0.1), metric = "subspace", k = k)
    expect_true(is(similarities, "FERRET_Similarities"))
    expect_equal(unname(similarities@ingroup[as.character(seq(1.6, 2.0, 0.1))]), as.numeric(rep(NA, 5)))
    expect_equal(unname(similarities@ingroup[as.character(seq(1.6, 2.0, 0.1))]), as.numeric(rep(NA, 5)))
    expect_equal(length(which(is.na(similarities@ingroup[as.character(seq(0.1, 1.5, 0.1))]))), 0)
    expect_equal(length(which(is.na(similarities@outgroup[as.character(seq(0.1, 1.5, 0.1))]))), 0)
    
    # Check that we get the same repeating values if we extend below the range of both network weights.
    lowRange <- seq(-1, 0, 0.1)
    similarities <- ComputeSimilarities(results = FERRET_ResultsObj, comparisons = FERRET_Comparisons,
                                        cutoffs = lowRange, metric = "jaccard", k = k)
    expect_true(is(similarities, "FERRET_Similarities"))
    expect_equal(range(similarities@ingroup)[1], range(similarities@ingroup)[2])
    expect_equal(range(similarities@outgroup)[1], range(similarities@outgroup)[2])
    similarities <- ComputeSimilarities(results = FERRET_ResultsObj, comparisons = FERRET_Comparisons,
                                        cutoffs = lowRange, metric = "degree", k = k)
    expect_true(is(similarities, "FERRET_Similarities"))
    expect_equal(range(similarities@ingroup)[1], range(similarities@ingroup)[2])
    expect_equal(range(similarities@outgroup)[1], range(similarities@outgroup)[2])
    similarities <- ComputeSimilarities(results = FERRET_ResultsObj, comparisons = FERRET_Comparisons,
                                        cutoffs = lowRange, metric = "modularity", k = k)
    expect_true(is(similarities, "FERRET_Similarities"))
    expect_equal(range(similarities@ingroup)[1], range(similarities@ingroup)[2])
    expect_equal(range(similarities@outgroup)[1], range(similarities@outgroup)[2])
    similarities <- ComputeSimilarities(results = FERRET_ResultsObj, comparisons = FERRET_Comparisons,
                                        cutoffs = lowRange, metric = "subspace", k = k)
    expect_true(is(similarities, "FERRET_Similarities"))
    expect_equal(range(similarities@ingroup)[1], range(similarities@ingroup)[2])
    expect_equal(range(similarities@outgroup)[1], range(similarities@outgroup)[2])
    
    # Check that we still get values for a single cutoff.
    singleCutoff <- 0.5
    similarities <- ComputeSimilarities(results = FERRET_ResultsObj, comparisons = FERRET_Comparisons,
                                        cutoffs = singleCutoff, metric = "jaccard", k = k)
    expect_true(is(similarities, "FERRET_Similarities"))
    expect_equal(length(similarities@ingroup), 1)
    expect_equal(length(similarities@outgroup), 1)
    similarities <- ComputeSimilarities(results = FERRET_ResultsObj, comparisons = FERRET_Comparisons,
                                        cutoffs = singleCutoff, metric = "degree", k = k)
    expect_true(is(similarities, "FERRET_Similarities"))
    expect_equal(length(similarities@ingroup), 1)
    expect_equal(length(similarities@outgroup), 1)
    similarities <- ComputeSimilarities(results = FERRET_ResultsObj, comparisons = FERRET_Comparisons,
                                        cutoffs = singleCutoff, metric = "modularity", k = k)
    expect_true(is(similarities, "FERRET_Similarities"))
    expect_equal(length(similarities@ingroup), 1)
    expect_equal(length(similarities@outgroup), 1)
    similarities <- ComputeSimilarities(results = FERRET_ResultsObj, comparisons = FERRET_Comparisons,
                                        cutoffs = singleCutoff, metric = "subspace", k = k)
    expect_true(is(similarities, "FERRET_Similarities"))
    expect_equal(length(similarities@ingroup), 1)
    expect_equal(length(similarities@outgroup), 1)
    
    # Check that we get an error if we don't specify any cutoffs.
    expect_error(ComputeSimilarities(results = FERRET_ResultsObj, comparisons = FERRET_Comparisons,
                                     cutoffs = NULL, metric = "jaccard", k = k), 
                 "You must specify at least one cutoff of numeric type!")
    expect_error(ComputeSimilarities(results = FERRET_ResultsObj, comparisons = FERRET_Comparisons,
                                     cutoffs = c(), metric = "jaccard", k = k), 
                 "You must specify at least one cutoff of numeric type!")
    
    # Check that we get an error for cutoffs not of numeric type.
    expect_error(ComputeSimilarities(results = FERRET_ResultsObj, comparisons = FERRET_Comparisons,
                                     cutoffs = "cut", metric = "jaccard", k = k), 
                 "You must specify at least one cutoff of numeric type!")
  })
})
test_that("[FERRET] ComputeRobustnessForOneEdgeType() function yields expected results",{
  net1 <- data.frame(source = c(rep("tf1", 4), rep("tf2", 4), rep("tf3", 3), rep("tf4", 3), "tf3"),
                     target = c(rep(c("gene1", "gene2", "gene3", "gene4"), 2),
                                rep(c("gene5", "gene6", "gene7"), 2), "gene3"),
                     score = seq(0.1, 1.5, 0.1))
  rownames(net1) <- paste(net1$source, net1$target, sep = "_")
  net2 <- data.frame(source = c(rep("tf1", 4), rep("tf2", 4), rep("tf3", 3), rep("tf4", 3), "tf4"),
                     target = c(rep(c("gene1", "gene2", "gene3", "gene4"), 2),
                                rep(c("gene5", "gene6", "gene7"), 2), "gene3"),
                     score = seq(0.1, 1.5, 0.1))
  rownames(net2) <- paste(net2$source, net2$target, sep = "_")
  net3 <- data.frame(source = c("tf1", "tf2", rep("tf3", 3), rep("tf4", 5)),
                     target = c("gene1", "gene2", rep(c("gene3", "gene4", "gene5"), 2), "gene6", "gene7"),
                     score = seq(0.1, 1.0, 0.1))
  rownames(net3) <- paste(net3$source, net3$target, sep = "_")
  net4 <- data.frame(source = c(rep("tf1", 4), rep("tf2", 4), rep("tf3", 3), rep("tf4", 3), "tf4"),
                     target = c(rep(c("gene1", "gene2", "gene3", "gene4"), 2),
                                rep(c("gene5", "gene6", "gene7"), 2), "gene4"),
                     score = seq(0.1, 1.5, 0.1))
  rownames(net4) <- paste(net4$source, net4$target, sep = "_")
  net5 <- data.frame(source = c("tf1", "tf2", rep("tf3", 3), rep("tf4", 5)),
                     target = c("gene2", "gene1", rep(c("gene3", "gene4", "gene5"), 2), "gene6", "gene7"),
                     score = seq(0.1, 1.0, 0.1))
  k=5
  rownames(net5) <- paste(net5$source, net5$target, sep = "_")
  resultList <- list(net1,net2,net3,net4,net5)
  names(resultList) <- c("net1", "net2", "net3", "net4", "net5")
  FERRET_ResultsObj <- methods::new("FERRET_Results",results = resultList, directory = "tmp")
  FERRET_Comparisons <- methods::new("FERRET_Comparisons", ingroup = c(new("FERRET_Comparison", source = "net1", target = "net2"),
                                                                            new("FERRET_Comparison", source = "net1", target = "net4")),
                                          outgroup = c(new("FERRET_Comparison", source = "net1", target = "net3"),
                                                       new("FERRET_Comparison", source = "net1", target = "net5")))
  
  # Check that robustness is formulated correctly.
  suppressWarnings({
    expect_equal(names(ComputeRobustnessForOneEdgeType(results = FERRET_ResultsObj, comparisons = FERRET_Comparisons, k = k,
                                                       plotCurve = FALSE)@auc),
                 c("Jaccard", "Degree", "Modularity", "Subspace"))
    expect_equal(names(ComputeRobustnessForOneEdgeType(results = FERRET_ResultsObj, comparisons = FERRET_Comparisons, k = k,
                                                       plotCurve = FALSE,
                                                       metric = c("jaccard", "modularity", "subspace"))@auc),
                 c("Jaccard", "Modularity", "Subspace"))
    expect_equal(names(ComputeRobustnessForOneEdgeType(results = FERRET_ResultsObj, comparisons = FERRET_Comparisons, k = k,
                                                       plotCurve = FALSE,
                                                       metric = "modularity")@auc),
                 c("Modularity"))
  })
  
  # Check errors when incorrect metrics are input.
  expect_error(ComputeRobustnessForOneEdgeType(results = FERRET_ResultsObj, comparisons = FERRET_Comparisons, k = k,
                                               plotCurve = FALSE,
                                               metric = c()),
               "Metric must be one or more of the following: jaccard, modularity, subspace, degree")
  expect_error(ComputeRobustnessForOneEdgeType(results = FERRET_ResultsObj, comparisons = FERRET_Comparisons, k = k,
                                               plotCurve = FALSE,
                                               metric = NULL),
               "Metric must be one or more of the following: jaccard, modularity, subspace, degree")
  expect_error(ComputeRobustnessForOneEdgeType(results = FERRET_ResultsObj, comparisons = FERRET_Comparisons, k = k,
                                               plotCurve = FALSE,
                                               metric = "mySuperAwesomeMetric"),
               "Metric must be one or more of the following: jaccard, modularity, subspace, degree")
  
  # Check that robustness is high when we expect it to be high.
  suppressWarnings({
    FERRET_Comparisons_High <- methods::new("FERRET_Comparisons", ingroup = c(new("FERRET_Comparison", source = "net1", target = "net2"),
                                                                         new("FERRET_Comparison", source = "net1", target = "net4")),
                                       outgroup = c(new("FERRET_Comparison", source = "net1", target = "net3"),
                                                    new("FERRET_Comparison", source = "net1", target = "net5")))
    robustnessHigh <- ComputeRobustnessForOneEdgeType(results = FERRET_ResultsObj, comparisons = FERRET_Comparisons_High, k = k,
                                                  plotCurve = FALSE)@auc
    FERRET_Comparisons_Low <- methods::new("FERRET_Comparisons", outgroup = c(new("FERRET_Comparison", source = "net1", target = "net2"),
                                                                         new("FERRET_Comparison", source = "net1", target = "net4")),
                                       ingroup = c(new("FERRET_Comparison", source = "net1", target = "net3"),
                                                    new("FERRET_Comparison", source = "net1", target = "net5")))
    robustnessLow <- ComputeRobustnessForOneEdgeType(results = FERRET_ResultsObj, comparisons = FERRET_Comparisons_Low, k = k,
                                                      plotCurve = FALSE)@auc
    expect_gt(robustnessHigh["Jaccard"], robustnessLow["Jaccard"])
    expect_gt(robustnessHigh["Degree"], robustnessLow["Degree"])
    expect_gt(robustnessHigh["Modularity"], robustnessLow["Modularity"])
    expect_gt(robustnessHigh["Subspace"], robustnessLow["Subspace"])
  })
})
test_that("[FERRET] ComputeRobustnessAUC() function yields expected results",{
  net1 <- data.frame(source = c(rep("tf1", 4), rep("tf2", 4), rep("tf3", 3), rep("tf4", 3), "tf3"),
                     target = c(rep(c("gene1", "gene2", "gene3", "gene4"), 2),
                                rep(c("gene5", "gene6", "gene7"), 2), "gene3"),
                     score = seq(0.1, 1.5, 0.1))
  rownames(net1) <- paste(net1$source, net1$target, sep = "_")
  net2 <- data.frame(source = c(rep("tf1", 4), rep("tf2", 4), rep("tf3", 3), rep("tf4", 3), "tf4"),
                     target = c(rep(c("gene1", "gene2", "gene3", "gene4"), 2),
                                rep(c("gene5", "gene6", "gene7"), 2), "gene3"),
                     score = seq(0.1, 1.5, 0.1))
  rownames(net2) <- paste(net2$source, net2$target, sep = "_")
  net3 <- data.frame(source = c("tf1", "tf2", rep("tf3", 3), rep("tf4", 5)),
                     target = c("gene1", "gene2", rep(c("gene3", "gene4", "gene5"), 2), "gene6", "gene7"),
                     score = seq(0.1, 1.0, 0.1))
  rownames(net3) <- paste(net3$source, net3$target, sep = "_")
  net4 <- data.frame(source = c(rep("tf1", 4), rep("tf2", 4), rep("tf3", 3), rep("tf4", 3), "tf4"),
                     target = c(rep(c("gene1", "gene2", "gene3", "gene4"), 2),
                                rep(c("gene5", "gene6", "gene7"), 2), "gene4"),
                     score = seq(0.1, 1.5, 0.1))
  rownames(net4) <- paste(net4$source, net4$target, sep = "_")
  net5 <- data.frame(source = c("tf1", "tf2", rep("tf3", 3), rep("tf4", 5)),
                     target = c("gene2", "gene1", rep(c("gene3", "gene4", "gene5"), 2), "gene6", "gene7"),
                     score = seq(0.1, 1.0, 0.1))
  rownames(net5) <- paste(net5$source, net5$target, sep = "_")
  resultList <- list(net1,net2,net3,net4,net5)
  names(resultList) <- c("net1", "net2", "net3", "net4", "net5")
  FERRET_ResultsObj <- methods::new("FERRET_Results",results = resultList, directory = "tmp")
  FERRET_Comparisons <- methods::new("FERRET_Comparisons", ingroup = c(new("FERRET_Comparison", source = "net1", target = "net2"),
                                                                            new("FERRET_Comparison", source = "net1", target = "net4")),
                                          outgroup = c(new("FERRET_Comparison", source = "net1", target = "net3"),
                                                       new("FERRET_Comparison", source = "net1", target = "net5")))
  
  # Check handling of invalid inputs.
  expect_error(ComputeRobustnessAUC(results = FERRET_ResultsObj, comparisons = "hi")@auc,
               paste("Inputs must be of type FERRET_Results and FERRET_Comparisons. Construct these using",
                     "LoadResults() and BuildComparisonObject(), respectively."), fixed = TRUE)
  expect_error(ComputeRobustnessAUC(results = "hi", comparisons = FERRET_Comparisons)@auc,
              paste("Inputs must be of type FERRET_Results and FERRET_Comparisons. Construct these using",
                    "LoadResults() and BuildComparisonObject(), respectively."), fixed = TRUE)
  
  # Check that par is reset after plotting.
  suppressWarnings({
    robustness <- ComputeRobustnessAUC(results = FERRET_ResultsObj, comparisons = FERRET_Comparisons)@auc
    robustnessPercentile <- ComputeRobustnessAUC(results = FERRET_ResultsObj, comparisons = FERRET_Comparisons,
                                                 mode = "Percentile", numberOfCutoffs = 5)@auc
  })
  expect_equal(par()$mfrow, c(1,1))
  
  # Compute robustness when all edges are inhibitory and check that it is the same.
  net1Neg <- net1
  net1Neg[,3] <- -1 * net1[,3]
  net1Neg[,1] <- paste0(net1[,1], "neg")
  net1Neg[,2] <- paste0(net1[,2], "neg")
  net2Neg <- net2
  net2Neg[,3] <- -1 * net2[,3]
  net2Neg[,1] <- paste0(net2[,1], "neg")
  net2Neg[,2] <- paste0(net2[,2], "neg")
  net3Neg <- net3
  net3Neg[,3] <- -1 * net3[,3]
  net3Neg[,1] <- paste0(net3[,1], "neg")
  net3Neg[,2] <- paste0(net3[,2], "neg")
  net4Neg <- net4
  net4Neg[,3] <- -1 * net4[,3]
  net4Neg[,1] <- paste0(net4[,1], "neg")
  net4Neg[,2] <- paste0(net4[,2], "neg")
  net5Neg <- net5
  net5Neg[,3] <- -1 * net5[,3]
  net5Neg[,1] <- paste0(net5[,1], "neg")
  net5Neg[,2] <- paste0(net5[,2], "neg")
  resultList <- list(net1Neg,net2Neg,net3Neg,net4Neg,net5Neg)
  names(resultList) <- c("net1Neg", "net2Neg", "net3Neg", "net4Neg", "net5Neg")
  FERRET_ResultsObj <- methods::new("FERRET_Results",results = resultList, directory = "tmp")
  FERRET_Comparisons <- methods::new("FERRET_Comparisons", ingroup = c(new("FERRET_Comparison", source = "net1Neg", target = "net2Neg"),
                                                                       new("FERRET_Comparison", source = "net1Neg", target = "net4Neg")),
                                     outgroup = c(new("FERRET_Comparison", source = "net1Neg", target = "net3Neg"),
                                                  new("FERRET_Comparison", source = "net1Neg", target = "net5Neg")))
  expect_error(ComputeRobustnessAUC(results = FERRET_ResultsObj, comparisons = FERRET_Comparisons)@auc,
               paste("If negatives are present, valid values for interpretationOfNegative in a FERRET_Results object are",
                     "'inhibitory' and 'poor'"))
  FERRET_ResultsObj@interpretationOfNegative <- "someOtherInterpretation"
  expect_error(ComputeRobustnessAUC(results = FERRET_ResultsObj, comparisons = FERRET_Comparisons)@auc,
               paste("If negatives are present, valid values for interpretationOfNegative in a FERRET_Results object are",
                     "'inhibitory' and 'poor'"))
  FERRET_ResultsObj@interpretationOfNegative <- "inhibitory"
  suppressWarnings({
    robustnessNeg <- ComputeRobustnessAUC(results = FERRET_ResultsObj, comparisons = FERRET_Comparisons, plotCurve = FALSE)@auc
    robustnessNegPercentile <- ComputeRobustnessAUC(results = FERRET_ResultsObj, comparisons = FERRET_Comparisons, plotCurve = FALSE,
                                                    mode = "Percentile", numberOfCutoffs = 5)@auc
  })
  expect_equal(robustness, robustnessNeg)
  expect_equal(robustnessPercentile, robustnessNegPercentile)
  
  # Check that, if we concatenate the negatives and positives, we get the same robustness.
  net1Comp <- rbind(net1, net1Neg)
  net2Comp <- rbind(net2, net2Neg)
  net3Comp <- rbind(net3, net3Neg)
  net4Comp <- rbind(net4, net4Neg)
  net5Comp <- rbind(net5, net5Neg)
  resultList <- list(net1Comp,net2Comp,net3Comp,net4Comp,net5Comp)
  names(resultList) <- c("net1Comp", "net2Comp", "net3Comp", "net4Comp", "net5Comp")
  FERRET_ResultsObj <- methods::new("FERRET_Results",results = resultList, directory = "tmp")
  FERRET_Comparisons <- methods::new("FERRET_Comparisons", ingroup = c(new("FERRET_Comparison", source = "net1Comp", target = "net2Comp"),
                                                                       new("FERRET_Comparison", source = "net1Comp", target = "net4Comp")),
                                     outgroup = c(new("FERRET_Comparison", source = "net1Comp", target = "net3Comp"),
                                                  new("FERRET_Comparison", source = "net1Comp", target = "net5Comp")))
  FERRET_ResultsObj@interpretationOfNegative <- "inhibitory"
  suppressWarnings({
    robustnessComp <- ComputeRobustnessAUC(results = FERRET_ResultsObj, comparisons = FERRET_Comparisons, plotCurve = FALSE)@auc
    robustnessCompPercentile <- ComputeRobustnessAUC(results = FERRET_ResultsObj, comparisons = FERRET_Comparisons, plotCurve = FALSE,
                                                    mode = "Percentile", numberOfCutoffs = 5)@auc
  })
  expect_equal(robustness, robustnessComp)
  expect_equal(robustnessPercentile, robustnessCompPercentile)
  
  # Check that scaling works for negative edges that represent a poor relationship.
  net1 <- data.frame(source = c(rep("tf1", 4), rep("tf2", 4), rep("tf3", 3), rep("tf4", 3), "tf3"),
                     target = c(rep(c("gene1", "gene2", "gene3", "gene4"), 2),
                                rep(c("gene5", "gene6", "gene7"), 2), "gene3"),
                     score = seq(0.0, 1.4, 0.1))
  rownames(net1) <- paste(net1$source, net1$target, sep = "_")
  net2 <- data.frame(source = c(rep("tf1", 4), rep("tf2", 4), rep("tf3", 3), rep("tf4", 3), "tf4"),
                     target = c(rep(c("gene1", "gene2", "gene3", "gene4"), 2),
                                rep(c("gene5", "gene6", "gene7"), 2), "gene3"),
                     score = seq(0.0, 1.4, 0.1))
  rownames(net2) <- paste(net2$source, net2$target, sep = "_")
  net3 <- data.frame(source = c("tf1", "tf2", rep("tf3", 3), rep("tf4", 5)),
                     target = c("gene1", "gene2", rep(c("gene3", "gene4", "gene5"), 2), "gene6", "gene7"),
                     score = seq(0.0, 0.9, 0.1))
  rownames(net3) <- paste(net3$source, net3$target, sep = "_")
  net4 <- data.frame(source = c(rep("tf1", 4), rep("tf2", 4), rep("tf3", 3), rep("tf4", 3), "tf4"),
                     target = c(rep(c("gene1", "gene2", "gene3", "gene4"), 2),
                                rep(c("gene5", "gene6", "gene7"), 2), "gene4"),
                     score = seq(0.0, 1.4, 0.1))
  rownames(net4) <- paste(net4$source, net4$target, sep = "_")
  net5 <- data.frame(source = c("tf1", "tf2", rep("tf3", 3), rep("tf4", 5)),
                     target = c("gene2", "gene1", rep(c("gene3", "gene4", "gene5"), 2), "gene6", "gene7"),
                     score = seq(0.0, 0.9, 0.1))
  rownames(net5) <- paste(net5$source, net5$target, sep = "_")
  resultList <- list(net1,net2,net3,net4,net5)
  names(resultList) <- c("net1", "net2", "net3", "net4", "net5")
  FERRET_ResultsObj <- methods::new("FERRET_Results",results = resultList, directory = "tmp")
  FERRET_Comparisons <- methods::new("FERRET_Comparisons", ingroup = c(new("FERRET_Comparison", source = "net1", target = "net2"),
                                                                       new("FERRET_Comparison", source = "net1", target = "net4")),
                                     outgroup = c(new("FERRET_Comparison", source = "net1", target = "net3"),
                                                  new("FERRET_Comparison", source = "net1", target = "net5")))
  FERRET_ResultsObj@interpretationOfNegative <- "poor"
  suppressWarnings({
    robustness <- ComputeRobustnessAUC(results = FERRET_ResultsObj, comparisons = FERRET_Comparisons, plotCurve = FALSE)@auc
    robustnessPercentile <- ComputeRobustnessAUC(results = FERRET_ResultsObj, comparisons = FERRET_Comparisons, plotCurve = FALSE,
                                                 mode = "Percentile", numberOfCutoffs = 5)@auc
  })
  net1Poor <- net1
  net1Poor[,3] <- net1[,3] - 0.3
  net2Poor <- net2
  net2Poor[,3] <- net2[,3] - 0.3
  net3Poor <- net3
  net3Poor[,3] <- net3[,3] - 0.3
  net4Poor <- net4
  net4Poor[,3] <- net4[,3] - 0.3
  net5Poor <- net5
  net5Poor[,3] <- net5[,3] - 0.3
  resultList <- list(net1Poor,net2Poor,net3Poor,net4Poor,net5Poor)
  names(resultList) <- c("net1Poor", "net2Poor", "net3Poor", "net4Poor", "net5Poor")
  FERRET_ResultsObj <- methods::new("FERRET_Results",results = resultList, directory = "tmp")
  FERRET_Comparisons <- methods::new("FERRET_Comparisons", ingroup = c(new("FERRET_Comparison", source = "net1Poor", target = "net2Poor"),
                                                                       new("FERRET_Comparison", source = "net1Poor", target = "net4Poor")),
                                     outgroup = c(new("FERRET_Comparison", source = "net1Poor", target = "net3Poor"),
                                                  new("FERRET_Comparison", source = "net1Poor", target = "net5Poor")))
  FERRET_ResultsObj@interpretationOfNegative <- "poor"
  suppressWarnings({
    robustnessPoor <- ComputeRobustnessAUC(results = FERRET_ResultsObj, comparisons = FERRET_Comparisons, plotCurve = FALSE)@auc
    robustnessPoorPercentile <- ComputeRobustnessAUC(results = FERRET_ResultsObj, comparisons = FERRET_Comparisons, plotCurve = FALSE,
                                                     mode = "Percentile", numberOfCutoffs = 5)@auc
  })
  expect_equal(robustness, robustnessPoor)
  expect_equal(robustnessPercentile, robustnessPoorPercentile)
})
test_that("[FERRET] WriteRobustnessAUC() and ReadRobustnessAUC() functions yield expected results",{
  # Check that we get the correct number of outgroup and ingroup results.
  ingroupVec <- c(0.1, 0.2, 0.3)
  names(ingroupVec) <- c(0.4, 0.5, 0.6)
  robustnessAUC <- methods::new("FERRET_ROC_AUC", auc = c(Metric1 = 0.6, Metric2 = 0.4, Metric3 = 0.8), 
                                roc = list("Metric1" = methods::new("FERRET_Similarities", ingroup = ingroupVec, outgroup = ingroupVec / 10),
                                           "Metric2" = methods::new("FERRET_Similarities", ingroup = ingroupVec, outgroup = ingroupVec / 10),
                                           "Metric3" = methods::new("FERRET_Similarities", ingroup = ingroupVec, outgroup = ingroupVec / 10)))
  WriteRobustnessAUC(robustnessAUC, "tmp.csv")
  robustnessWritten <- read.csv("tmp.csv", row.names = 1)
  expect_equal(colnames(robustnessWritten), c("Metric1", "Metric2", "Metric3", "cutoff", "inOrOut"))
  expect_equal(rownames(robustnessWritten), c("AUC", "IngroupROC_0.4", "IngroupROC_0.5", "IngroupROC_0.6",
                                             "OutgroupROC_0.4", "OutgroupROC_0.5", "OutgroupROC_0.6"))
  expect_true(is.na(robustnessWritten["AUC", "cutoff"]))
  expect_true(is.na(robustnessWritten["AUC", "inOrOut"]))
  expect_equal(robustnessWritten[2:4,"Metric1"], robustnessWritten[5:7, "Metric1"] * 10)
  expect_equal(robustnessWritten[2:7,"Metric1"], robustnessWritten[2:7,"Metric2"])
  expect_equal(robustnessWritten[2:7,"Metric3"], robustnessWritten[2:7,"Metric2"])
  expect_equal(robustnessWritten[2:7,"cutoff"], rep(c(0.4, 0.5, 0.6), 2))
  expect_true(is.na(robustnessWritten[1,"cutoff"]))
  expect_true(is.na(robustnessWritten[1,"inOrOut"]))
  robustnessAUCRead <- ReadRobustnessAUC("tmp.csv")
  expect_equal(robustnessAUC, robustnessAUCRead)
  unlink("tmp.csv")
  
  # Check the same when there is only one cutoff.
  ingroupVec <- 0.1
  names(ingroupVec) <- 0.4
  robustnessAUC <- methods::new("FERRET_ROC_AUC", auc = c(Metric1 = 0.6, Metric2 = 0.4, Metric3 = 0.8), 
                                roc = list("Metric1" = methods::new("FERRET_Similarities", ingroup = ingroupVec, outgroup = ingroupVec / 10),
                                           "Metric2" = methods::new("FERRET_Similarities", ingroup = ingroupVec, outgroup = ingroupVec / 10),
                                           "Metric3" = methods::new("FERRET_Similarities", ingroup = ingroupVec, outgroup = ingroupVec / 10)))
  WriteRobustnessAUC(robustnessAUC, "tmp.csv")
  robustnessWritten <- read.csv("tmp.csv", row.names = 1)
  expect_equal(colnames(robustnessWritten), c("Metric1", "Metric2", "Metric3", "cutoff", "inOrOut"))
  expect_equal(rownames(robustnessWritten), c("AUC", "IngroupROC_0.4","OutgroupROC_0.4"))
  expect_true(is.na(robustnessWritten["AUC", "cutoff"]))
  expect_true(is.na(robustnessWritten["AUC", "inOrOut"]))
  expect_equal(robustnessWritten[2:3,"Metric1"], c(0.1, 0.01))
  expect_equal(robustnessWritten[2:3,"Metric1"], robustnessWritten[2:3,"Metric2"])
  expect_equal(robustnessWritten[2:3,"Metric3"], robustnessWritten[2:3,"Metric2"])
  expect_equal(robustnessWritten[2:3,"cutoff"], rep(0.4, 2))
  expect_equal(robustnessWritten[2:3,"inOrOut"], c("Ingroup", "Outgroup"))
  robustnessAUCRead <- ReadRobustnessAUC("tmp.csv")
  expect_equal(robustnessAUC, robustnessAUCRead)
  unlink("tmp.csv")
  
  # Check the same when there is only one metric.
  ingroupVec <- c(0.1, 0.2, 0.3)
  names(ingroupVec) <- c(0.4, 0.5, 0.6)
  robustnessAUC <- methods::new("FERRET_ROC_AUC", auc = c(Metric1 = 0.6), 
                                roc = list("Metric1" = methods::new("FERRET_Similarities", ingroup = ingroupVec, outgroup = ingroupVec / 10)))
  WriteRobustnessAUC(robustnessAUC, "tmp.csv")
  robustnessWritten <- read.csv("tmp.csv", row.names = 1)
  expect_equal(colnames(robustnessWritten), c("Metric1", "cutoff", "inOrOut"))
  expect_equal(rownames(robustnessWritten), c("AUC", "IngroupROC_0.4", "IngroupROC_0.5", "IngroupROC_0.6",
                                             "OutgroupROC_0.4", "OutgroupROC_0.5", "OutgroupROC_0.6"))
  expect_true(is.na(robustnessWritten["AUC", "cutoff"]))
  expect_true(is.na(robustnessWritten["AUC", "inOrOut"]))
  expect_equal(robustnessWritten[5:7,"Metric1"], robustnessWritten[2:4,"Metric1"] / 10)
  expect_equal(robustnessWritten[2:7,"cutoff"], rep(c(0.4, 0.5, 0.6), 2))
  expect_equal(robustnessWritten[2:7,"inOrOut"], c(rep("Ingroup", 3), rep("Outgroup", 3)))
  robustnessAUCRead <- ReadRobustnessAUC("tmp.csv")
  expect_equal(robustnessAUC, robustnessAUCRead)
  unlink("tmp.csv")
  
  # Check the same when there is only one metric and one cutoff.
  ingroupVec <- 0.1
  names(ingroupVec) <- 0.4
  robustnessAUC <- methods::new("FERRET_ROC_AUC", auc = c(Metric1 = 0.6), 
                                roc = list("Metric1" = methods::new("FERRET_Similarities", ingroup = ingroupVec, outgroup = ingroupVec / 10)))
  WriteRobustnessAUC(robustnessAUC, "tmp.csv")
  robustnessWritten <- read.csv("tmp.csv", row.names = 1)
  expect_equal(colnames(robustnessWritten), c("Metric1", "cutoff", "inOrOut"))
  expect_equal(rownames(robustnessWritten), c("AUC", "IngroupROC_0.4", "OutgroupROC_0.4"))
  expect_true(is.na(robustnessWritten["AUC", "cutoff"]))
  expect_true(is.na(robustnessWritten["AUC", "inOrOut"]))
  expect_equal(robustnessWritten[3,"Metric1"], robustnessWritten[2,"Metric1"] / 10)
  expect_equal(robustnessWritten[2:3,"cutoff"], rep(0.4, 2))
  expect_equal(robustnessWritten[2:3,"inOrOut"], c("Ingroup", "Outgroup"))
  robustnessAUCRead <- ReadRobustnessAUC("tmp.csv")
  expect_equal(robustnessAUC, robustnessAUCRead)
  unlink("tmp.csv")
  
  suppressWarnings({
    # Check that an error is returned when the wrong type of object is passed.
    expect_error(WriteRobustnessAUC("hi", "tmp.csv"), "Error: Input must be an object of type FERRET_ROC_AUC")
  
    # Check that an error is thrown when trying to write to an invalid file.
    expect_error(WriteRobustnessAUC(robustnessAUC, "../aFolderThatDoesntExist/tmp.csv"), "Cannot write to file ../aFolderThatDoesntExist/tmp.csv")
  })
  
})
test_that("[FERRET] GetDifferentialNetworks() function yields expected results",{
  # Check that difference is correct with an example.
  net1 <- data.frame(source = c("tf1", "tf2", "tf3"),
                     target = c("gene1", "gene2", "gene3"),
                     score = c(0.1, 1.5, 0.1))
  rownames(net1) <- paste(net1$source, net1$target, sep = "_")
  net2 <- net1
  net2$score <- net1$score / 2
  net2Neg <- net2
  net2Neg$score <- -1 * net2Neg$score
  net1Copy <- net1
  net2Copy <- net2
  resultList <- list(net1,net2,net1Copy, net2Copy)
  names(resultList) <- c("net1", "net2", "net1Copy", "net2Copy")
  FERRET_ResultsObj <- methods::new("FERRET_Results",results = resultList, directory = "tmp", interpretationOfNegative = "poor")
  FERRET_Comparisons <- methods::new("FERRET_Comparisons", ingroup = list(new("FERRET_Comparison", source = "net1", target = "net1Copy")),
                                     outgroup = list(new("FERRET_Comparison", source = "net1", target = "net2")))
  expect_equal(GetDifferentialNetworks(FERRET_ResultsObj, FERRET_Comparisons), net2)
   
  # Check the same when only one edge overlaps.
  netOneEdgeOverlap <- data.frame(source = c("tf1", "tfA", "tfB"),
                                          target = c("gene1", "geneA", "geneB"),
                                          score = c(0.1, 1.5, 0.1))
  rownames(netOneEdgeOverlap) <- paste(netOneEdgeOverlap$source, netOneEdgeOverlap$target, sep = "_")
  netOneEdgeExpected <- data.frame(source = c("tf2", "tf3"),
                                   target = c("gene2", "gene3"),
                                   score = c(1.5, 0.1))
  rownames(netOneEdgeExpected) <- paste(netOneEdgeExpected$source, netOneEdgeExpected$target, sep = "_")
  resultList <- list(net1,netOneEdgeOverlap, net1Copy)
  net1Copy <- net1
  names(resultList) <- c("net1", "netOneEdgeOverlap", "net1Copy")
  FERRET_ResultsObj <- methods::new("FERRET_Results",results = resultList, directory = "tmp", interpretationOfNegative = "poor")
  FERRET_Comparisons <- methods::new("FERRET_Comparisons", ingroup = list(new("FERRET_Comparison", source = "net1", target = "net1Copy")),
                                     outgroup = list(new("FERRET_Comparison", source = "net1", target = "netOneEdgeOverlap")))
  expect_equal(GetDifferentialNetworks(FERRET_ResultsObj, FERRET_Comparisons), netOneEdgeExpected)
  
  # Check that we get the expected result when there are both positives and negatives.
  net1WithNeg <- net1
  net1WithNeg[2, "score"] <- -1 * net1[2, "score"]
  rownames(net1WithNeg) <- paste(net1WithNeg$source, net1WithNeg$target, sep = "_")
  net2ExpDiff <- net2
  net2ExpDiff[2, "score"] <- -1 * net1[2, "score"]
  rownames(net2ExpDiff) <- paste(net2ExpDiff$source, net2ExpDiff$target, sep = "_")
  net2NegExpDiff <- net2Neg[2,]
  net2NegExpDiff[1, "score"] <- 0.75
  rownames(net2NegExpDiff) <- paste(net2NegExpDiff$source, net2NegExpDiff$target, sep = "_")
  net1WithNegCopy <- net1WithNeg
  resultList <- list(net1WithNeg,net2,net1WithNegCopy,net2Copy)
  names(resultList) <- c("net1WithNeg", "net2", "net1WithNegCopy","net2Copy")
  FERRET_ResultsObj <- methods::new("FERRET_Results",results = resultList, directory = "tmp", interpretationOfNegative = "inhibitory")
  FERRET_Comparisons <- methods::new("FERRET_Comparisons", ingroup = c(new("FERRET_Comparison", source = "net1WithNeg", target = "net1WithNegCopy")),
                                     outgroup = c(new("FERRET_Comparison", source = "net1WithNeg", target = "net2")))
  expect_equal(GetDifferentialNetworks(FERRET_ResultsObj, FERRET_Comparisons), net2ExpDiff)
  FERRET_Comparisons <- methods::new("FERRET_Comparisons", ingroup = c(new("FERRET_Comparison", source = "net2", target = "net2Copy")),
                                     outgroup = c(new("FERRET_Comparison", source = "net2", target = "net1WithNeg")))
  expect_equal(GetDifferentialNetworks(FERRET_ResultsObj, FERRET_Comparisons), net2NegExpDiff)
  net2ExpDiff <- net1WithNeg[c(1,3),] 
  net2ExpDiff[,"score"] <- c(0.05, 0.05)
  net2NegExpDiff <- net2[2,]
  net2NegExpDiff[,"score"] <- net2NegExpDiff[,"score"] + 1.5
  FERRET_ResultsObj <- methods::new("FERRET_Results",results = resultList, directory = "tmp", interpretationOfNegative = "poor")
  FERRET_Comparisons <- methods::new("FERRET_Comparisons", ingroup = c(new("FERRET_Comparison", source = "net1WithNeg", target = "net1WithNegCopy")),
                                     outgroup = c(new("FERRET_Comparison", source = "net1WithNeg", target = "net2")))
  expect_equal(GetDifferentialNetworks(FERRET_ResultsObj, FERRET_Comparisons), net2ExpDiff)
  FERRET_Comparisons <- methods::new("FERRET_Comparisons", ingroup = c(new("FERRET_Comparison", source = "net2", target = "net2Copy")),
                                     outgroup = c(new("FERRET_Comparison", source = "net2", target = "net1WithNeg")))
  expect_equal(GetDifferentialNetworks(FERRET_ResultsObj, FERRET_Comparisons), net2NegExpDiff)
  
  # Check that the difference from a network to itself has all 0 edges.
  emptyDf <- as.data.frame(matrix(nrow = 0,ncol = 3))
  colnames(emptyDf) <- colnames(net1)
  resultList <- list(net1,net1Copy)
  names(resultList) <- c("net1","net1Copy")
  FERRET_ResultsObj <- methods::new("FERRET_Results",results = resultList, directory = "tmp")
  FERRET_Comparisons <- methods::new("FERRET_Comparisons", ingroup = c(new("FERRET_Comparison", source = "net1", target = "net1Copy")),
                                     outgroup = c(new("FERRET_Comparison", source = "net1", target = "net1Copy")))
  expect_equal(nrow(GetDifferentialNetworks(FERRET_ResultsObj, FERRET_Comparisons)), nrow(emptyDf))
  expect_equal(colnames(GetDifferentialNetworks(FERRET_ResultsObj, FERRET_Comparisons)), colnames(emptyDf))
  
  # Check that the difference from a network to the empty network returns the original network.
  suppressWarnings({
    resultList <- list(net1, emptyDf, net1Copy)
    names(resultList) <- c("net1", "emptyDf", "net1Copy")
    FERRET_ResultsObj <- methods::new("FERRET_Results",results = resultList, directory = "tmp")
    FERRET_Comparisons <- methods::new("FERRET_Comparisons", ingroup = c(new("FERRET_Comparison", source = "net1", target = "net1Copy")),
                                       outgroup = c(new("FERRET_Comparison", source = "net1", target = "emptyDf")))
    expect_equal(GetDifferentialNetworks(FERRET_ResultsObj, FERRET_Comparisons), net1)
    
    # Check that the difference between two empty networks is the empty network.
    FERRET_Comparisons <- methods::new("FERRET_Comparisons", ingroup = c(new("FERRET_Comparison", source = "emptyDf", target = "emptyDf")),
                                       outgroup = c(new("FERRET_Comparison", source = "emptyDf", target = "emptyDf")))
    expect_equal(nrow(GetDifferentialNetworks(FERRET_ResultsObj, FERRET_Comparisons)), nrow(emptyDf))
    expect_equal(colnames(GetDifferentialNetworks(FERRET_ResultsObj, FERRET_Comparisons)), colnames(emptyDf))
  })
})
test_that("[FERRET] GetCommonNetwork() function yields expected results",{
  # Check that the common network is correct with an example.
  net1 <- data.frame(source = c("tf1", "tf2", "tf3"),
                     target = c("gene1", "gene2", "gene3"),
                     score = c(0.1, 1.5, 0.1))
  rownames(net1) <- paste(net1$source, net1$target, sep = "_")
  net2 <- net1
  net2$score <- net1$score / 2
  net2Neg <- net2
  net2Neg$score <- -1 * net2Neg$score
  net1Copy <- net1
  net2Copy <- net2
  resultList <- list(net1, net2, net1Copy, net2Copy)
  names(resultList) <- c("net1", "net2", "net1Copy", "net2Copy")
  FERRET_ResultsObj <- methods::new("FERRET_Results",results = resultList, directory = "tmp")
  FERRET_Comparisons <- methods::new("FERRET_Comparisons", ingroup = c(new("FERRET_Comparison", source = "net1", target = "net1")),
                                     outgroup = c(new("FERRET_Comparison", source = "net1", target = "net2")))
  expect_equal(GetCommonNetwork(FERRET_ResultsObj, FERRET_Comparisons), net2)
  FERRET_Comparisons <- methods::new("FERRET_Comparisons", ingroup = c(new("FERRET_Comparison", source = "net2", target = "net2")),
                                     outgroup = c(new("FERRET_Comparison", source = "net2", target = "net1")))
  expect_equal(GetCommonNetwork(FERRET_ResultsObj, FERRET_Comparisons), net2)
  
  # Check the same when only one edge overlaps.
  netOneEdgeOverlap <- data.frame(source = c("tf1", "tfA", "tfB"),
                                  target = c("gene1", "geneA", "geneB"),
                                  score = c(0.1, 1.5, 0.1))
  rownames(netOneEdgeOverlap) <- paste(netOneEdgeOverlap$source, netOneEdgeOverlap$target, sep = "_")
  netOneEdgeExpected <- data.frame(source = "tf1",
                                   target = "gene1",
                                   score = 0.1)
  rownames(netOneEdgeExpected) <- paste(netOneEdgeExpected$source, netOneEdgeExpected$target, sep = "_")
  resultList <- list(net1, netOneEdgeOverlap, net1Copy)
  names(resultList) <- c("net1", "netOneEdgeOverlap", "net1Copy")
  FERRET_ResultsObj <- methods::new("FERRET_Results",results = resultList, directory = "tmp")
  FERRET_Comparisons <- methods::new("FERRET_Comparisons", ingroup = c(new("FERRET_Comparison", source = "net1", target = "net1Copy")),
                                     outgroup = c(new("FERRET_Comparison", source = "net1", target = "netOneEdgeOverlap")))
  expect_equal(GetCommonNetwork(FERRET_ResultsObj, FERRET_Comparisons), netOneEdgeExpected)
  
  # Check that the shared network between a network and itself is the same network.
  FERRET_Comparisons <- methods::new("FERRET_Comparisons", ingroup = c(new("FERRET_Comparison", source = "net1", target = "net1Copy")),
                                     outgroup = c(new("FERRET_Comparison", source = "net1", target = "net1Copy")))
  expect_equal(GetCommonNetwork(FERRET_ResultsObj, FERRET_Comparisons), net1)
  
  # Check that the network is correct with both inhibitory and poor when there are negatives.
  net1WithNeg <- net1
  net1WithNeg[2, "score"] <- -1 * net1[2, "score"]
  rownames(net1WithNeg) <- paste(net1WithNeg$source, net1WithNeg$target, sep = "_")
  net2ExpDiff <- net1WithNeg[c(1,3),] 
  net2ExpDiff[,"score"] <- c(0.05, 0.05)
  rownames(net2ExpDiff) <- paste(net2ExpDiff$source, net2ExpDiff$target, sep = "_")
  net1WithNegCopy <- net1WithNeg
  resultList <- list(net1WithNeg,net2,net1WithNegCopy,net2Copy)
  names(resultList) <- c("net1WithNeg", "net2", "net1WithNegCopy","net2Copy")
  FERRET_ResultsObj <- methods::new("FERRET_Results",results = resultList, directory = "tmp", interpretationOfNegative = "inhibitory")
  FERRET_Comparisons <- methods::new("FERRET_Comparisons", ingroup = c(new("FERRET_Comparison", source = "net1WithNeg", target = "net1WithNegCopy")),
                                     outgroup = c(new("FERRET_Comparison", source = "net1WithNeg", target = "net2")))
  expect_equal(GetCommonNetwork(FERRET_ResultsObj, FERRET_Comparisons), net2ExpDiff)
  FERRET_Comparisons <- methods::new("FERRET_Comparisons", ingroup = c(new("FERRET_Comparison", source = "net2", target = "net2Copy")),
                                     outgroup = c(new("FERRET_Comparison", source = "net2", target = "net1WithNeg")))
  expect_equal(GetCommonNetwork(FERRET_ResultsObj, FERRET_Comparisons), net2ExpDiff)
  net2ExpDiff <- net1WithNeg[c(1,3),] 
  net2ExpDiff[,"score"] <- c(1.55, 1.55)
  FERRET_ResultsObj <- methods::new("FERRET_Results",results = resultList, directory = "tmp", interpretationOfNegative = "poor")
  FERRET_Comparisons <- methods::new("FERRET_Comparisons", ingroup = c(new("FERRET_Comparison", source = "net1WithNeg", target = "net1WithNegCopy")),
                                     outgroup = c(new("FERRET_Comparison", source = "net1WithNeg", target = "net2")))
  expect_equal(GetCommonNetwork(FERRET_ResultsObj, FERRET_Comparisons), net2ExpDiff)
  FERRET_Comparisons <- methods::new("FERRET_Comparisons", ingroup = c(new("FERRET_Comparison", source = "net2", target = "net2Copy")),
                                     outgroup = c(new("FERRET_Comparison", source = "net2", target = "net1WithNeg")))
  expect_equal(GetCommonNetwork(FERRET_ResultsObj, FERRET_Comparisons), net2ExpDiff)
  
  # Check that the shared network with an empty network is empty.
  suppressWarnings({
    emptyDf <- as.data.frame(matrix(nrow = 0,ncol = 3))
    colnames(emptyDf) <- colnames(net1)
    resultList <- list(net1, emptyDf, net1Copy)
    names(resultList) <- c("net1", "emptyDf", "net1Copy")
    FERRET_ResultsObj <- methods::new("FERRET_Results",results = resultList, directory = "tmp")
    FERRET_Comparisons <- methods::new("FERRET_Comparisons", ingroup = c(new("FERRET_Comparison", source = "net1", target = "net1Copy")),
                                       outgroup = c(new("FERRET_Comparison", source = "net1", target = "emptyDf")))
    expect_equal(nrow(GetCommonNetwork(FERRET_ResultsObj, FERRET_Comparisons)), nrow(emptyDf))
    expect_equal(colnames(GetCommonNetwork(FERRET_ResultsObj, FERRET_Comparisons)), colnames(emptyDf))
    FERRET_Comparisons <- methods::new("FERRET_Comparisons", ingroup = c(new("FERRET_Comparison", source = "emptyDf", target = "emptyDf")),
                                       outgroup = c(new("FERRET_Comparison", source = "emptyDf", target = "emptyDf")))
    expect_equal(nrow(GetCommonNetwork(FERRET_ResultsObj, FERRET_Comparisons)), nrow(emptyDf))
    expect_equal(colnames(GetCommonNetwork(FERRET_ResultsObj, FERRET_Comparisons)), colnames(emptyDf))  
  })
})
test_that("[FERRET] GetOverlapWithPrior() function yields expected results",{
  # Test an example with partial overlap.
  net1 <- data.frame(source = c("tf1", "tf2", "tf3"),
                     target = c("gene1", "gene2", "gene3"),
                     score = c(0.1, 1.5, 0.1))
  prior <- data.frame(source = c("tf1", "tfA", "tfB"),
                     target = c("gene1", "gene2", "gene3"),
                     score = c(1, 1, 1))
  expect_equal(GetOverlapWithPrior(net1, prior), 0.1 / 1.7)

  # Test that result = 1 when the network and the prior match completely.
  prior <- data.frame(source = c("tf1", "tf2", "tf3"),
                      target = c("gene1", "gene2", "gene3"),
                      score = c(1, 1, 1))
  expect_equal(GetOverlapWithPrior(net1, prior), 1)
  
  # Test with no overlap.
  prior <- data.frame(source = c("tfC", "tfA", "tfB"),
                      target = c("gene1", "gene2", "gene3"),
                      score = c(1, 1, 1))
  expect_equal(GetOverlapWithPrior(net1, prior), 0)
  
  # Test with an empty prior.
  emptyDf <- as.data.frame(matrix(nrow = 0,ncol = 3))
  colnames(emptyDf) <- colnames(net1)
  emptyDf$source <- as.character(emptyDf$source)
  emptyDf$target <- as.character(emptyDf$target)
  emptyDf$score <- as.numeric(emptyDf$score)
  expect_error(GetOverlapWithPrior(net1, emptyDf), "Prior cannot be empty!")
  
  # Test with an empty network.
  expect_equal(GetOverlapWithPrior(emptyDf, prior), 0)
  
  # Test when both prior and network are empty.
  expect_error(GetOverlapWithPrior(emptyDf, emptyDf), "Prior cannot be empty!")
  
})
test_that("[FERRET] GetFullNetworkFromPCANetwork() function yields expected results",{
  # Test with an example.
  net1 <- data.frame(source = c("PC1", "PC2", "PC3"),
                     target = c("PC2", "PC3", "PC4"),
                     score = c(0.1, 1.5, 0.1))
  pca <- data.frame(PC1 = c(1, 0.5, 0.25, 0),
                    PC2 = c(0, 1, 0.5, 0.25),
                    PC3 = c(0.25, 0, 1, 0.5),
                    PC4 = c(0.5, 0.25, 0, 1))
  rownames(pca) <- c("gene1", "gene2", "gene3", "gene4")
  from1to2 <- 0.1 * 1 * 1 + 1.5 * 0 * 0 + 0.1 * 0.25 * 0.25
  from4to3 <- 0.1 * 0 * 0.5 + 1.5 * 0.25 * 1 + 0.1 * 0.5 * 0
  fullNet <- GetFullNetworkFromPCANetwork(net1, pca)
  expect_equal(fullNet["gene1_gene2", 3], from1to2)
  expect_equal(fullNet["gene4_gene3", 3], from4to3)
  expect_equal(max(fullNet[,3]), fullNet["gene2_gene3",3])
})