test_testMotif <- function(){
  data(pandaToyData)
  data(pandaResult)
  regnet = slot(pandaResult,"regNet")
  p1 <- with(pandaToyData, testMotif(regnet, mode="augment", motif, expression, ppi, hamming=1))
  p2 <- with(pandaToyData, testMotif(regnet, mode="remove", motif, expression, ppi, hamming=1))
  checkTrue("ggplot" %in% class(p1) & "ggplot" %in% class(p2))
}
