test_plotZ <- function(){
  data("pandaToyData")
  panda.res1 <- with(pandaToyData, panda(motif, expression, ppi, hamming=1))
  panda.res2 <- with(pandaToyData, panda(motif, expression + 
                                           rnorm(prod(dim(expression)),sd=5), ppi, hamming=1))
  p <- plotZ(panda.res1, panda.res2, addLine=FALSE)
  checkTrue("ggplot" %in% class(p))
}
