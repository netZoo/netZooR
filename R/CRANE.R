#' Find the robust nodes in ALPACA community using CRANE
#'
#' @param input same input for alpaca: first column TF, second column Genes, third column edge weights from baseline condition, fourth column edge weights from disease condition.
#' @param alp alpca object in list format (output from alpaca package)
#' @param alpha alpha paramter perturbs each edge weights
#' @param beta beta parameter perturbs the strength of each node. Set this to 0 if you want nodes to have node strength identical to the orignal network.
#' @param iteration Number of CRANE distributions to create. Higher value leads to better ranking but longer runtime.
#' @param isParallel TRUE = use Multithread / FALSE = do not use Multithread
#'
#' @import foreach
#' @import doParallel
#' @import parallel
#'
#' @return list of data frames
#' @export
#' @examples
#' \dontrun{
#'
#' input=cbind(nonAng,ang[,3])
#' alp=alpaca(input,NULL,verbose = F)
#' alpListObject=alpacaCrane(input, alp, isParallel = T)
#'
#' }


alpacaCrane = function(input,alp,alpha=0.1,beta=0,iteration=30,isParallel=F){
  elist1=input[,1:3]
  tfcheck=paste(unique(elist1[,1]),"A",sep="_")
  gcheck=paste(unique(elist1[,2]),"B",sep="_")
  member=alpacaGetMember(alp)
  
  print("Converting to Adjecency Matrix...")
  A=elistToAdjMat(elist1,isBipartite=T)
  ridx=order(rownames(A))
  cidx=order(colnames(A))
  elist1=adjMatToElist(A[ridx,cidx])
  
  print("Generating CRANE Networks, this may take a long time (upto 10-20 minutes for Networks with 19000+ nodes) ...")
  if (isParallel){
    cores=detectCores()
    cl <- makeCluster(cores[1]-1)
    registerDoParallel(cl)
    exportFunctions=c("craneBipartite","elistRemoveTags","adjMatToElist","jutterDegree")
    elist.df=foreach(i=1:iteration,.export = exportFunctions, .combine='cbind', .packages=c("netZooR"),.verbose = F) %dopar% {
      craneBipartite(A,alpha=alpha,beta=beta)[,3]
    }
    stopCluster(cl)
  }else{
    elist.df=NULL
    for (i in 1:iteration){
      elist.df=cbind(elist.df,craneBipartite(A,alpha=alpha,beta=beta)[,3])
    }
  }
  gc()
  
  print("Creating Test Distribution...")
  elist1=elistAddTags(elist1)
  
  print("Detecting communities in control network...")
  ctrl.pos = elist1[elist1[,3]>=0,1:3]
  ctrl.elist = data.frame(red=ctrl.pos[,2],blue=ctrl.pos[,1],weights=ctrl.pos[,3])
  ctrl.condor = createCondorObject(ctrl.elist)
  ctrl.condor = condorCluster(ctrl.condor,project=F)
  ctrl.memb = c(ctrl.condor$red.memb[,2],ctrl.condor$blue.memb[,2])
  names(ctrl.memb) = c(as.character(ctrl.condor$red.memb[,1]),as.character(ctrl.condor$blue.memb[,1]))
  
  
  print("Comuting Differential Scores...")
  df=NULL
  for (i in 1:iteration){
    net.table=cbind(elist1,elist.df[,i])
    pos.table = net.table[intersect(which(net.table[,3]>=0),which(net.table[,4]>=0)),]
    pos.graph = graph.edgelist(as.matrix(pos.table[,1:2]),directed=T)
    
    if (length(setdiff(V(pos.graph)$name,names(ctrl.memb)))>0){
      uncounted <- setdiff(V(pos.graph)$name,names(ctrl.memb))
      unc.memb <- sample(1:max(ctrl.memb),length(uncounted),replace=T)
      names(unc.memb) <- uncounted
      ctrl.memb <- c(ctrl.memb,unc.memb)
    }
    if (i %% round(iteration/5)==0){
      print(paste("Computing Differential Score",(i/iteration)*100,"%"))
    }
    dwbm = computeDWBMmat.mscale(pos.table,ctrl.memb[V(pos.graph)$name])
    
    miss_tf=setdiff(tfcheck,rownames(dwbm))
    miss_g=setdiff(gcheck,colnames(dwbm))
    if (length(miss_tf)>0){
      temp=array(0,dim = c(length(miss_tf),ncol(dwbm)))
      rownames(temp)=miss_tf
      dwbm=rbind(dwbm,temp)
    }
    if (length(miss_g)>0){
      temp=array(0,dim = c(nrow(dwbm),length(miss_g)))
      colnames(temp)=miss_g
      dwbm=cbind(dwbm,temp)
    }
    
    temp=alpacaComputeDifferentialScoreFromDWBM(dwbm,member)
    
    if (i>1 & !identical(rownames(df),names(temp))){
      print("error")
      break
    }else{
      df=cbind(df,temp)
    }
  }
  
  isSig=c()
  tstat=c()
  if (identical(names(alp[[2]]),rownames(df))){
    for (i in 1:nrow(df)){
      if (rownames(df)[i]==""){
        isSig=c(isSig,NA)
        next
      }
      dist=df[i,]
      if (length(dist)>1){
        temp=t.test(dist, mu=alp[[2]][i], alternative = "less",conf.level=0.95)
        isSig=c(isSig,temp$p.value)
        temp2=temp$statistic[[1]]
        names(temp2)=rownames(df)[i]
        tstat=c(tstat,temp2)
      }else{
        temp.v=1
        names(temp.v)=rownames(df)[i]
        isSig=c(isSig,temp.v)
        temp.v[]=10e+200
        tstat=c(tstat,temp.v)
      }
    }
  }
  alp[[3]]=isSig
  alp[[4]]=tstat
  
  alp.out=alpacaObjectToDfList(alp)
  
  return(alp.out)
}


#' Pertrubs the bipartite network with fixed node strength
#'
#' @param df Adjacency Matrix or Edge list
#' @param alpha alpha paramter perturbs each edge weights
#' @param beta beta parameter perturbs the strength of each node. Set this to 0 if you want nodes to have node strength identical to the orignal network.
#' @param getAdj TRUE = this will return adjacency matrix instead of edge list
#' @param randomStart FALSE = initialize the first row with completely random edge weights.
#'
#' @return edge list
#' @export
#' @examples
#' \dontrun{
#'
#' # Using Edge list as input
#' elist=craneBipartite(nonAng)
#' elist=craneBipartite(nonAng,alpha=0.3)
#'
#' # Using Edge list as input and Adjcency Matrix as output
#' adjMatrix=craneBipartite(nonAng,alpha=0.1,getAdj=T)
#'
#' # Using Edge list as input and Adjcency Matrix as output
#' A=elistToAdjMat(nonAng)
#' elist=craneBipartite(A)
#'
#' }

craneBipartite= function(df,alpha=0.1,beta=0,getAdj=F,randomStart=F){
  if (isElist(df)){
    print("Converting Edgelist to Adj Matrix")
    A=elistToAdjMat(df,isBipartite=T)
  }else{
    A=df
  }
  
  print(paste("Applying Alpha =",alpha))
  rown=nrow(A)
  coln=ncol(A)
  
  # randomize nodes
  ridx=sample(1:rown,rown,replace = F)
  cidx=sample(1:coln,coln,replace = F)
  A=A[ridx,cidx]
  
  rlt_A=array(0,dim=dim(A))
  colnames(rlt_A)=colnames(A)
  rownames(rlt_A)=rownames(A)
  
  # beta implimentation: beta = 1 models subsample data
  if (beta!=0){
    print("Beta on TFs")
    tfD=rowSums(A)
    tfD.jit=jutterDegree(tfD,beta,beta_slope = T)
    tfD.delta=tfD.jit-tfD
    correction=tfD.delta/coln
    A=A+correction
    print("Beta on Genes")
    A=t(A)
    tfD=rowSums(A)
    tfD.jit=jutterDegree(tfD,beta,beta_slope = F)
    tfD.delta=tfD.jit-tfD
    correction=tfD.delta/rown
    A=t(A)
  }
  
  tfD=rowSums(A)
  gD=colSums(A)
  
  omin=min(A)
  omax=max(A)
  
  if (alpha >0){
    rlt_A[1,]=A[1,]
    if (randomStart){
      rlt_A[1,]=rnorm(length(A[1,]),mean(A[1,]),2)
    }
    print("Constructing Iteratuvely Perturbed Network")
    current_degree=rep(0,coln)
    for (i in 2:rown){
      if (i %% 100 ==0){
        print(i)
      }
      #Randomly add delta to previous the Rows by alpha
      temp=rlt_A[(i-1),]
      n=length(temp)
      avg=mean(temp)
      std=sd(temp)
      cur_min=min(temp)
      cur_max=max(temp)
      delta=rnorm(n,mean=0,sd=std)
      temp=temp+alpha*delta
      temp=(temp/sum(temp))*tfD[(i-1)]
      goal_degree=colSums(A[1:i,])
      k=0
      # Backtrack to Control min max
      while(1){
        if (k >150){
          break
        }
        temp_degree=current_degree+temp
        y=goal_degree-temp_degree
        
        lidx=which(y<omin)
        ridx=which(y>omax)
        idx=c(lidx,ridx)
        lcor=y[lidx]-omin
        rcor=y[ridx]-omax
        if(length(idx)==0){
          break
        }
        if (length(lidx)>0){
          temp[lidx]=temp[lidx]-abs(lcor)
        }
        if (length(ridx)>0){
          temp[ridx]=temp[ridx]+abs(rcor)
        }
        temp=(temp/sum(temp))*tfD[(i-1)]
        k=k+1
      }
      #Update Final Edges
      rlt_A[(i-1),]=temp
      current_degree=current_degree+temp
      y=goal_degree-current_degree
      rlt_A[i,]=y
    }
  }
  
  if(!isTRUE(all.equal(rowSums(rlt_A),tfD,check.attributes = F, use.names = F))){
    print("ROW ERROR Alpha limit has been reached")
    return(NULL)
  }
  if(!isTRUE(all.equal(colSums(rlt_A),gD,check.attributes = F, use.names = F))){
    print("COLUMN ERROR Alpha limit has been reached")
    return(NULL)
  }
  
  print("Sorting Nodes")
  ridx=order(rownames(rlt_A))
  cidx=order(colnames(rlt_A))
  if (getAdj){
    return(rlt_A[ridx,cidx])
  }
  elist=adjMatToElist(rlt_A[ridx,cidx])
  elist=elistRemoveTags(elist)
  return(elist)
}

#' Pertrubs the unipartite network with fixed node strength from adjacency matrix
#'
#' @param A Adjacency Matrix
#' @param alpha alpha paramter perturbs each edge weights
#' @param isSelfLoop TRUE/FALSE if self loop exists. co-expression matrix will have a self-loop of 1. Thus TRUE
#'
#' @return adjacency matrix
#' @export
#'
craneUnipartite = function(A,alpha=0.1,isSelfLoop=F){
  print(paste("Applying Alpha =",alpha))
  
  rown=nrow(A)
  coln=ncol(A)
  
  # randomize nodes
  ridx=sample(1:rown,rown,replace = F)
  A=A[ridx,ridx]
  
  rlt_A=array(0,dim=dim(A))
  colnames(rlt_A)=colnames(A)
  rownames(rlt_A)=rownames(A)
  
  Aorig=A
  oavg=mean(Aorig)
  diag(Aorig)=mean(Aorig)
  A[lower.tri(A)]=0
  
  rowD=rowSums(A)
  colD=colSums(A)
  gD=colSums(A)
  
  omin=min(A)
  oomin=min(A[A>0])
  omax=max(A)
  randn=rown-2
  if (alpha >0){
    rlt_A[1,]=A[1,]
    print("Constructing Iteratuvely Perturbed Network")
    current_degree=rep(0,coln)
    for (i in 2:(rown-1)){
      if (i %% 100 ==0){
        print(i)
      }
      #Randomly add delta to previous the Rows by alpha
      temp=rlt_A[(i-1),]
      n=randn
      avg=mean(temp[(i+1):rown])
      std=sd(Aorig[(i-1),])
      cur_min=min(temp[(i+1):rown])
      cur_max=max(temp[(i+1):rown])
      delta=rnorm(n,mean=0,sd=std)
      temp[(i+1):rown]=temp[(i+1):rown]+alpha*delta
      temp[(i+1):rown][is.na(temp[(i+1):rown])]=oomin
      if (any(temp[(i+1):rown]<omin)){
        temp[(i+1):rown]=temp[(i+1):rown]+abs(min(temp[(i+1):rown]))
      }
      temp[(i+1):rown]=(temp[(i+1):rown]/sum(temp[(i+1):rown]))*(rowD[(i-1)]-sum(temp[1:i]))
      goal_degree=colSums(A[1:i,])
      k=0
      # Backtrack to Control min max
      while(1){
        if (k >150){
          break
        }
        temp_degree=current_degree+temp
        y=goal_degree-temp_degree
        
        lidx=which(y<omin)
        ridx=which(y>omax)
        idx=c(lidx,ridx)
        lcor=y[lidx]-omin
        rcor=y[ridx]-omax
        if(length(idx)==0){
          break
        }
        if (length(lidx)>0){
          temp[lidx]=temp[lidx]-abs(lcor)
        }
        if (length(ridx)>0){
          temp[ridx]=temp[ridx]+abs(rcor)
        }
        temp[(i+1):rown][is.na(temp[(i+1):rown])]=oomin
        if (any(temp[(i+1):rown]<omin)){
          temp[(i+1):rown]=temp[(i+1):rown]+abs(min(temp[(i+1):rown]))+omin
        }
        temp[(i+1):rown]=(temp[(i+1):rown]/sum(temp[(i+1):rown]))*(rowD[(i-1)]-sum(temp[1:i]))
        k=k+1
      }
      #Final Update Last Edges
      rlt_A[(i-1),]=temp
      current_degree=current_degree+temp
      y=goal_degree-current_degree
      rlt_A[i,]=y
      randn=randn-1
    }
  }
  
  if (isSelfLoop){
    rlt_A[nrow(rlt_A),ncol(rlt_A)]=A[nrow(rlt_A),ncol(rlt_A)]
  }
  
  # Non-Iterative version of the algorithm (not recommended)
  if(!isTRUE(all.equal(rowSums(rlt_A),rowD,check.attributes = F, use.names = F))){
    print("ROW ERROR Alpha limit has been reached")
    return(NULL)
  }
  if(!isTRUE(all.equal(colSums(rlt_A),colD,check.attributes = F, use.names = F))){
    print("COLUMN ERROR Alpha limit has been reached")
    return(NULL)
  }
  
  for (i in 1:nrow(rlt_A)){
    for (j in 1:ncol(rlt_A)){
      rlt_A[j,i]=rlt_A[i,j]
    }
  }
  
  print("Sorting Nodes")
  ridx=order(rownames(rlt_A))
  cidx=order(colnames(rlt_A))
  
  return(rlt_A[ridx,cidx])
}

#' Compute Differential modularity score from differential modularity matrix
#'
#' This functions takes the precomputed differential modularity matrix and the genLouvain membership to compute the differential modularity score.
#' @param dwbm differential modularity matrix
#' @param louv.memb louvain community membership
#'
#' @return Vector of differntial modularity score
#'
alpacaComputeDifferentialScoreFromDWBM = function(dwbm,louv.memb){
  louv.Ascores <- NULL
  louv.Bscores <- NULL
  for (i in 1:max(louv.memb)){
    this.comm <- names(louv.memb)[louv.memb==i]
    this.tfs <- this.comm[grep("_A",this.comm)]
    this.genes <- this.comm[grep("_B",this.comm)]
    if (length(this.tfs)>1){
      tf.sums <- apply(as.matrix(dwbm[this.tfs,this.genes]),1,sum)
      gene.sums <- apply(as.matrix(dwbm[this.tfs,this.genes]),2,sum)
    } else {
      tf.sums <- sum(dwbm[this.tfs,this.genes])
      gene.sums <- dwbm[this.tfs,this.genes]
    }
    this.denom <- sum(dwbm[this.tfs,this.genes])
    louv.Ascores <- c(louv.Ascores,tf.sums/this.denom)
    louv.Bscores <- c(louv.Bscores,gene.sums/this.denom)
  }
  
  louv.scores <- c(louv.Ascores,louv.Bscores)
  
  return(louv.scores)
}

#' get the member vector from alpaca object
#'
#' @param alp alpaca object
#' @param target tf, gene, or all
#'
#' @return member vector
#'
alpacaGetMember = function(alp, target="all") {
  alp2 = alp[[1]]
  memb_vector = as.vector(alp2)
  names(memb_vector) = names(alp2)
  memb_vector=memb_vector[!is.na(memb_vector)]
  if (target =="tf"){
    return(memb_vector[grep("_A",names(memb_vector))])
  }else if (target == "gene"){
    return(memb_vector[grep("_B",names(memb_vector))])
  }else{
    return(memb_vector)
  }
}

#' creates condor object
#'
#' @param elist edge list
#'
#' @return condor object
condorCreateObject = function(elist){
  elist=elistAddTags(elist)
  elist=elist[elist[,3]>0,]
  elist=data.frame(red=elist[,2],blue=elist[,1],weights=elist[,3])
  elist=createCondorObject(elist,return.gcc=T)
  return(elist)
}

#' Run CONDOR clustering
#'
#' @param elist edge list
#' @param qscore TRUE = output qscore / FALSE = do not output qscore
#'
#' @return condor object
#'
condorRun = function(elist,qscore=F){
  cond.object=condorCreateObject(elist)
  cond.object=condorCluster(cond.object)
  if (qscore){
    cond.object=condorQscore(cond.object)
  }
  return(cond.object)
}

#' Adds "_A" to first column and "_B" to second column
#'
#' @param elist edge list
#'
#' @return edge list
#'
elistAddTags = function(elist){
  elist[,1]=paste(elist[,1],"A",sep="_")
  elist[,2]=paste(elist[,2],"B",sep="_")
  return(elist)
}

#' check if first two columns are identical
#'
#' @param elist1 edge list
#' @param elist2 edge list
#'
#' @return boolean
#'
elistIsEdgeOrderEqual=function(elist1,elist2){
  check1=identical(as.character(elist1[,1]),as.character(elist2[,1]))
  check2=identical(as.character(elist1[,2]),as.character(elist2[,2]))
  if (check1 & check2){
    return(TRUE)
  }else{
    return(FALSE)
  }
}

#' undo elistAddTags
#'
#' @param elist edge list
#'
#' @return edge list
#'
elistRemoveTags = function(elist){
  elist[,1]=gsub("_A","",elist[,1])
  elist[,2]=gsub("_B","",elist[,2])
  return(elist)
}

#' Sorts the edge list based on first two columns in alphabetical order
#'
#' @param elist edge list
#' @return edge list
#'
elistSort=function(elist){
  A=elistToAdjMat(elist,isBipartite=T)
  ridx=order(rownames(A))
  cidx=order(colnames(A))
  elist=adjMatToElist(A[ridx,cidx])
  return(elist)
}

#' Converts edge list to adjacency matrix
#'
#' @param elist edge list
#' @param isBipartite TRUE = for bipartite / FALSE = for unipartite
#'
#' @import igraph
#' @return Adjcency Matrix
#' @export
elistToAdjMat= function(elist,isBipartite=F){
  if (isTRUE(isBipartite)){
    elist=elistAddTags(elist)
  }
  colnames(elist)=c('from','to','weight')
  A_up = as.matrix(get.adjacency(graph.data.frame(elist),attr = 'weight')) # upper triangular matrix
  A = A_up[rownames(A_up) %in% elist[,1], colnames(A_up) %in% elist[,2]]
  return(A)
}

#' Converts alpaca output into list of data frames
#'
#' @param alp alpaca object
#'
#' @return list of data frames
#'
alpacaObjectToDfList= function(alp){
  alp.out=list()
  if (length(alp)<3){
    tf=alpacaGetMember(alp,target="tf")
    df=data.frame(tf,alp[[2]][names(tf)])
    colnames(df)=c("Membership","Differential Modularity Score")
    rownames(df)=gsub("_A","",rownames(df))
    alp.out[["TF"]]=df
    gene=alpacaGetMember(alp,target="gene")
    df=data.frame(gene,alp[[2]][names(gene)])
    colnames(df)=c("Membership","Differential Modularity Score")
    rownames(df)=gsub("_B","",rownames(df))
    alp.out[["Gene"]]=df
  }else{
    tf=alpacaGetMember(alp,target="tf")
    df=data.frame(tf,alp[[2]][names(tf)],alp[[3]][names(tf)],alp[[4]][names(tf)])
    colnames(df)=c("Membership","Differential Modularity Score","Pvalue","t-statistic")
    rownames(df)=gsub("_A","",rownames(df))
    alp.out[["TF"]]=df
    gene=alpacaGetMember(alp,target="gene")
    df=data.frame(gene,alp[[2]][names(gene)],alp[[3]][names(gene)],alp[[4]][names(gene)])
    colnames(df)=c("Membership","Differential Modularity Score","Pvalue","t-statistic")
    rownames(df)=gsub("_B","",rownames(df))
    alp.out[["Gene"]]=df
  }
  return(alp.out)
}

#' converts adjacency matrix to edge list
#'
#' @param adj_mat adjacency matrix
#'
#' @return edge list
#'
adjMatToElist = function(adj_mat){
  from=rep(row.names(adj_mat), ncol(adj_mat))
  to=rep(colnames(adj_mat), each=nrow(adj_mat))
  weight=as.numeric(unlist(c(adj_mat),use.names = F))
  elist=data.frame(from,to,weight)
  elist=elistRemoveTags(elist)
  return(elist)
}

#' CRANE Beta perturbation function. This function will add noice to the node strength sequence.
#'
#' @param nodeD Vector of node strength
#' @param beta beta
#' @param beta_slope TRUE=use predetermined slope to add noise / FALSE = use constant value for noise
#' @return vector with new strength distribution
#'
jutterDegree=function(nodeD,beta,beta_slope=T){
  if (beta==0){
    print("No BETA")
    return(nodeD)
  }
  print(paste("Applying Beta =",beta))
  if (beta<0){
    beta=abs(beta)
    correction=0.06
    const.sig=1
    yint=const.sig*beta
  }else{
    # Correction needed based on number of genes
    correction=0.025
    const.sig=400
    yint=const.sig*beta
  }
  if(beta_slope){
    beta.pos=beta*correction
    beta.neg=beta*correction*0.5
    print(paste("applying beta slope:",beta.pos))
    temp.sigma=nodeD
    temp.sigma[temp.sigma>=0]=abs(temp.sigma[temp.sigma>=0]*beta.pos)+yint
    temp.sigma[temp.sigma<0]=abs(temp.sigma[temp.sigma<0]*beta.neg)+yint
    temp.sigma=abs(temp.sigma)
  }else{
    print(paste("NO Slope Applied, Gene perturbation:", const.sig*beta))
    temp.sigma=const.sig*beta
  }
  nodeD=nodeD+rnorm(length(nodeD),0,temp.sigma)
  
  return(nodeD)
}

#' Check if data frame is an edge list
#'
#' @param df some data frame
#'
#' @return Boolean
#'
isElist = function(df){
  return(ncol(df)==3 & (is.factor(df[,1])|is.character(df[,1]))&(is.factor(df[,2])|is.character(df[,2]))&is.numeric(df[,3]) )
}
