setwd("~/Documents/GitHub/netZoo-devel/")



library(devtools)
library(roxygen2)

build()
document()

setwd("..")
install()

getwd()




library(netZoo)
?runCytoscapePlot
?sourcePPI
