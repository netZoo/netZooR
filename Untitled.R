setwd("~/Documents/GitHub/netZoo-draft/")



library(devtools)
library(roxygen2)

build()
document()

setwd("..")
install()






library(netZoo)
?runCytoscapePlot
?sourcePPI
