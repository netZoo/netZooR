pandaObj <- setClass("panda", slots=c("regNet","coregNet","coopNet"))
setMethod("show","panda",function(object){print.panda(object)})
