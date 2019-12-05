install.packages("devtools")
install.packages("BiocInstaller")
library(devtools)
library(BiocInstaller)

biocValid(fix=TRUE) 
update.packages(ask = FALSE, checkBuilt = TRUE)
devtools::install_deps(dep = T)
