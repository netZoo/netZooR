test_panda_r_vs_c <- function(){
  data("pandaToyData")
  panda.r <- with(pandaToyData, panda(motif, expression, ppi))
  reg.r <- reshape::melt.array(slot(panda.r,"regNet"))
  colnames(reg.r) <- c("TF", "Gene", "Z")
  reg.r <- reg.r[with(reg.r, order(TF, Gene)), ]
  
  data(pandaResultPairs)
  checkTrue( (1 - cor(pandaResultPairs$Z, reg.r$Z, method="spearman")) < 1e-10)
}
