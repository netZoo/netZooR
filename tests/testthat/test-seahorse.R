context("test SEAHORSE result")

test_that("seahorse function works", {
  expression_data <- read.csv("./GTEx_lung_expression_toydata.txt", sep="")
  phenotype_data <- read.csv("./GTEx_lung_phenotype_toydata.txt", sep="")
  phenotype_dictionary = c("categorical", "numeric")
  
  # Load KEGG pathways
  pathways_KEGG <- gmtPathways("./c2.cp.kegg.v2022.1.Hs.symbols.gmt")
  
  # Run seahorse
  results = seahorse(expression = expression_data, phenotype = phenotype_data, phenotype_dictionary, pathways = pathways_KEGG)
})