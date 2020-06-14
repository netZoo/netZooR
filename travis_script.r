#install.packages("devtools") #Only in TRAVIS
library(devtools)

update.packages(ask = FALSE, checkBuilt = TRUE)
devtools::install_deps(dep = T)
