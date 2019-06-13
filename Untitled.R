setwd("~/Documents/GitHub/netZoo_devel/")



library(devtools)
library(roxygen2)

build()
document()

setwd("..")
install()



?runPanda
