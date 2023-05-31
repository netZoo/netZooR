context("test SEAHORSE result")

test_that("seahorse function works", {
  # Simulate expression data
  expression_data = data.frame(matrix(rexp(1000, rate=.1), ncol=10, nrow = 100))
  rownames(expression_data) = paste("gene", 1:100, sep = "")
  colnames(expression_data) = paste("sample", 1:10, sep = "")
  
  # Simulate phenotypic data
  phenotype_data = data.frame(matrix(0, ncol=2, nrow = 10))
  colnames(phenotype_data) = c("sex", "height")
  rownames(phenotype_data) = colnames(expression_data)
  phenotype_data$sex = c(rep("male", nrow(phenotype_data)/2), rep("female", nrow(phenotype_data)/2))
  phenotype_data$height = 65 + sample.int(10, nrow(phenotype_data), replace = T)
  
  phenotype_dictionary = c("categorical", "numeric")
  
  # Create toy pathways
  pathways = list()
  pathways$pathway1 = sample(rownames(expression_data), 50)
  pathways$pathway2 = sample(rownames(expression_data), 30)
  pathways$pathway3 = sample(rownames(expression_data), 70)
  
  # Run seahorse
  results <- seahorse(expression_data, phenotype_data, phenotype_dictionary, pathways)
  })