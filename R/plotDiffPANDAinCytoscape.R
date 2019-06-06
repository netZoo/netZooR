
plotDiffPANDAinCytoscape <- function(net1,net2,network.name="DiffPANDAnetwork"){
  # reshape PANDA networks
  merged.net <- merge(net1,net2,by=c("tf","gene"))
  # edge-weight differences
  merged.net$'x-y' <- merged.net$force.x-merged.net$force.y
  merged.net$'y-x' <- merged.net$force.y-merged.net$force.x
  # CDF
  fnx <- ecdf(merged.net$force.x)
  fny <- ecdf(merged.net$force.y)
  fnxy <- ecdf(merged.net$'x-y')
  fnyx <- ecdf(merged.net$'y-x')
  sub.net1 <- merged.net[fnx(merged.net$force.x)*(fnxy(merged.net$'x-y'))>0.8,]
  sub.net1$cond1 <- "1"
  sub.net1 <- sub.net1[,c(1,2,3,4,9)]
  colnames(sub.net1) <- c("tf","gene","motif","force","cond1")
  sub.net2 <- merged.net[fny(merged.net$force.y)*(fnyx(merged.net$'y-x'))>0.8,]
  sub.net2$cond1 <- "0"
  sub.net2 <- sub.net2[,c(1,2,3,6,9)]
  colnames(sub.net2) <- c("tf","gene","motif","force","cond1")
  merge.sub.net <- rbind(sub.net1,sub.net2)
  
  # plot 
  # launch Cytoscape 3.6.1 or greater
  cytoscapePing ()
  cytoscapeVersionInfo ()
  
  colnames(merge.sub.net) <- c("source","target","interaction","weight","cond.")
  merge.sub.net$weight_transformed <- rank(merge.sub.net$weight)/length(merge.sub.net$weight)*10000
  if(nrow(merge.sub.net)>6666){
    message("The maximum network objects (nodes+edges) that Cytoscape could support varies in different memories.
            0-20000 objects are suggested by 512M memory size for Cytoscape to view. 
            For better view function in Cytoscape please reduce the size of PANDA network imported")
  }
  #create nodes arg for creating a cytoscape plot
  panda_nodes <- data.frame(id=c(merge.sub.net$source,merge.sub.net$target), group=rep(c("TF","Gene"), each=nrow(merge.sub.net)),stringsAsFactors=FALSE)
  # creat cytoscape from DataFrames
  createNetworkFromDataFrames(nodes=panda_nodes,edges=merge.sub.net, title=network.name, collection="DataFrame Example")

  
}


