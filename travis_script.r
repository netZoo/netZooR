install.packages("devtools")
library(devtools)

biocValid(fix=TRUE) 
update.packages(ask = FALSE, checkBuilt = TRUE)
devtools::install_deps(dep = T)
