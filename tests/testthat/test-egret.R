context("test EGRET result")

test_that("EGRET function works", {
  system("curl -O  https://netzoo.s3.us-east-2.amazonaws.com/netZooR/unittest_datasets/EGRET/toy_qbic.txt")
  system("curl -O  https://netzoo.s3.us-east-2.amazonaws.com/netZooR/unittest_datasets/EGRET/toy_genotype.vcf")
  system("curl -O  https://netzoo.s3.us-east-2.amazonaws.com/netZooR/unittest_datasets/EGRET/toy_motif_prior.txt")
  system("curl -O  https://netzoo.s3.us-east-2.amazonaws.com/netZooR/unittest_datasets/EGRET/toy_expr.txt")
  system("curl -O  https://netzoo.s3.us-east-2.amazonaws.com/netZooR/unittest_datasets/EGRET/toy_ppi_prior.txt")
  system("curl -O  https://netzoo.s3.us-east-2.amazonaws.com/netZooR/unittest_datasets/EGRET/toy_eQTL.txt")
  system("curl -O  https://netzoo.s3.us-east-2.amazonaws.com/netZooR/unittest_datasets/EGRET/toy_map.txt")
  system("curl -O  https://netzoo.s3.us-east-2.amazonaws.com/netZooR/unittest_datasets/EGRET/EgretGT.RData")
  
  qbic <- read.table(file = "toy_qbic.txt", header = FALSE)
  vcf <- read.table("toy_genotype.vcf", header = FALSE, sep = "\t", stringsAsFactors = FALSE, colClasses = c("character", "numeric", "character", "character", "character", "character", "character", "character", "character", "character"))
  motif <- read.table("toy_motif_prior.txt", sep = "\t", header = FALSE)
  expr <- read.table("toy_expr.txt", header = FALSE, sep = "\t", row.names = 1)
  ppi <- read.table("toy_ppi_prior.txt", header = FALSE, sep = "\t")
  qtl <- read.table("toy_eQTL.txt", header = FALSE)
  nameGeneMap <- read.table("toy_map.txt", header = FALSE)
  tag <- "my_toy_egret_run"

  #runEgret(qtl,vcf,qbic,motif,expr,ppi,nameGeneMap,tag)
  #load("EgretGT.RData")
  #load("my_toy_egret_run_egret.RData")
  #load("my_toy_egret_run_panda.RData")
  #expect_equal(regnetE,regnetEServer,tolerance=1e-5)
  #expect_equal(regnetP,regnetPServer,tolerance=1e-5)
})