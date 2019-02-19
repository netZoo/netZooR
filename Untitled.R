setwd("~/Documents/GitHub/netZoo-draft/")


library(devtools)
build()
install()

library(roxygen2)

setwd("./cats")
document()

library(netZoo)
?run

