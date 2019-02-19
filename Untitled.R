setwd("~/Documents/GitHub/netZoo-draft/")



 library(devtools)
build()
setwd("..")
install()

library(roxygen2)


document()

library(netZoo)
?runCytoscapePlot
?sourcePPI
