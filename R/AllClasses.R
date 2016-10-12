monsterAnalysis <- setClass("monsterAnalysis", slots=c("tm","nullTM","numGenes","numSamples"))
setMethod("show","monsterAnalysis",function(object){print.monsterAnalysis(object)})