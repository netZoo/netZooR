monsterObj <- setClass("monster", slots=c("tm","nullTM","numGenes","numSamples"))
setMethod("show","monster",function(object){print.monster(object)})
