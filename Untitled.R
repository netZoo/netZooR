setwd("~/Documents/GitHub/netZoo_devel/")



library(devtools)
library(roxygen2)

build()
document()

setwd("..")
install()



?runPanda
runLioness("~/Desktop/ToyData/ToyExpressionData.txt","~/Desktop/ToyData/ToyMotifData.txt","~/Desktop/ToyData/ToyPPIData.txt",rm_missing = F, )