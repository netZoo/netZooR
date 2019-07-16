## ----echo=TRUE,message=FALSE,warning=FALSE-------------------------------
library(ALPACA)

simp.mat <- read.table('Example_2comm.txt',header=T)

simp.alp <- alpaca(simp.mat,NULL,verbose=F)
simp.alp2 <- simp.alp[[1]]
simp.memb <- as.vector(simp.alp2)
names(simp.memb) <- names(simp.alp2)

simp.memb



