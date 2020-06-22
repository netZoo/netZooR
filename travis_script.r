install.packages("devtools", dep=T)
library(devtools)

update.packages(ask = FALSE, checkBuilt = TRUE)
devtools::install_deps(dep = T)
