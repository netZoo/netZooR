
plotLIONESSinCytoscape <- function(lioness.net,sample=c(1,2),network.name="LIONESSnetwork"){
  net1 <- data.frame(cbind(lioness.net[,c(1,2,(sample[1]+2) )]))
  colnames(net1) <- c("tf","gene","force")
  net2 <- data.frame(cbind(lioness.net[,c(1,2,(sample[2]+2) )]))
  colnames(net2) <- c("tf","gene","force")
  # reshape LINOESS networks
  merged.net <- merge(net1,net2,by=c("tf","gene"))
  # edge-weight differences
  merged.net$'x-y' <- merged.net$force.x-merged.net$force.y
  merged.net$'y-x' <- merged.net$force.y-merged.net$force.x
  
  # CDF
  fnx <- ecdf(merged.net$force.x)
  fny <- ecdf(merged.net$force.y)
  fnxy <- ecdf(merged.net$'x-y')
  fnyx <- ecdf(merged.net$'y-x')
  #indx1 <- fnx(merged.net$force.x)*(fnxy(merged.net$'x-y'))>0.8
  sub.net1 <- merged.net[fnx(merged.net$force.x)*(fnxy(merged.net$'x-y'))>0.8,]
  sub.net1$cond <- "1"
  sub.net1 <- sub.net1[,c(1,2,3,7)]
  colnames(sub.net1) <- c("tf","gene","force","cond")
  sub.net2 <- merged.net[fny(merged.net$force.y)*(fnyx(merged.net$'y-x'))>0.8,]
  sub.net2$cond <- "2"
  sub.net2 <- sub.net2[,c(1,2,4,7)]
  colnames(sub.net2) <- c("tf","gene","force","cond")
  merge.sub.net <- rbind(sub.net1,sub.net2)
  
  # no significant difference edges
  left.net1 <- merged.net[fnx(merged.net$force.x)*(fnxy(merged.net$'x-y'))<=0.8,]
  left.net1$cond <- "0"
  left.net1 <- left.net1[,c(1,2,3,7)]
  colnames(left.net1) <- c("tf","gene","force","cond")
  
  left.net2 <- merged.net[fny(merged.net$force.y)*(fnyx(merged.net$'y-x'))<=0.8,]
  left.net2$cond <- "0"
  left.net2 <- left.net2[,c(1,2,4,7)]
  colnames(left.net2) <- c("tf","gene","force","cond")
  
  merge.left.net <- merge(left.net1,left.net2,by=c("tf","gene"))
  merge.left.net$force <- (merge.left.net$force.x+merge.left.net$force.y)/2
  merge.left.net <- merge.left.net[,c(1,2,7,4)]
  colnames(merge.left.net) <- c("tf","gene","force","cond")
  merge.sub.net <- rbind(merge.sub.net,merge.left.net)
  
  
  # plot 
  # launch Cytoscape 3.6.1 or greater
  cytoscapePing ()
  cytoscapeVersionInfo ()
  
  colnames(merge.sub.net) <- c("source","target","weight","cond.")
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

