LoadUtil <- function() {
  
  library(R.matlab)
  library("dplyr")
  library(ggplot2)
  #library(devtools)
  library(RColorBrewer)
  library(cluster)
  library(pvclust)
  library(xtable)
  library(limma)
  library(plyr)
  library(ggplot2)
  library(lattice)
  library(RColorBrewer)
  library(scales)
  library(factoextra)
  library(gridExtra)
  library('proxy')
  library("energy")
  library(RColorBrewer)
  library(cluster)
  library(pvclust)
  library(FactoMineR)
  library(xtable)
  library(limma)
  library(plyr)
  library(ggplot2)
  library(car)
  library(lattice)
  library(RColorBrewer)
  library(scales)
  library(factoextra)
  library(gridExtra)
  library(NMF)
  library(apcluster)
  library(R.matlab)
  library(fpc)
  require(stats)
  library(igraph)
  library("corrplot")
  library("madness")
  library("factoextra")
  library("Rtsne")
  library(rmatio)
  library(umap)
  library("fossil")
  library("SummarizedExperiment")

}


LoadPollenDataSet <- function() {
  file="D:/Codes/DropOut/data/Pollen.mat"
  dataset=read.mat(file)
  data=dataset$data
  genes=dataset$genes
  classes =dataset$labels
  
  colnames(data)<-classes
  rownames(data)<-genes

  return (data)
}

LoadKolodDataSet <- function() {
  file="D:/Codes/DropOut/data/Kolod.mat"
  #704*10685
  dataset=read.mat(file)
  data=t(dataset$in_X)
  genes=as.character(dataset$Genes)
  classes =dataset$true_labs
  
  colnames(data)<-classes
  rownames(data)<-genes
  
  return (data)
}

LoadUsoskinDataSet <- function() {
  file="D:/Codes/DropOut/data/Usoskin.mat"
  dataset=read.mat(file)
  data=dataset$in_X
  classes =dataset$true_labs
  data<-t(data)
  colnames(data)<-classes
  return (data)
  
}

LoadKolodziejczykDataSet <- function() {
  file="D:/Codes/DropOut/data/Kolodziejczyk.mat"
  dataset=read.mat(file)
  data=dataset$data
  classes =dataset$labels
  gene=dataset$genes

  data<-as.data.frame(data)
  colnames(data)<-unname(classes)
  rownames(data)<-unname(gene)
  
 # nrow(data)
  #ncol(data)
  #length(classes)
  #length(gene)
  #View(data[1:10,1:10])
  return (data)
  
}


LoadTest_2_Kolod<- function() {
 
  load("D:/Codes/DropOut/data/Test_2_Kolod.RData")
  dataset<-Test_2_Kolod
  data=dataset$in_X
  classes =(dataset$true_labs)
  classes<-unlist(classes)
  #data<-t(data)
  colnames(data)<-classes

  
  return (data)
  
}




SCC_filtration <- function(data,celltype)
{
  
  pca<-prcomp(t(data))
  pcadata<-pca$x
  dis<-matrix(nrow=ncol(data),ncol=2)
  dis[,1]<-pcadata[,1]
  dis[,2]<-pcadata[,2]
  neardis<-dist(dis,p=2)
  dismatrix<-matrix(nrow=ncol(data),ncol=ncol(data))
  k<-1
  for (i in 1:(ncol(data)-1))
  {
    for (j in (i+1):ncol(data))
    {
      dismatrix[j,i]<-neardis[k]
      dismatrix[i,j]<-neardis[k]
      k<-k+1
    }
  }
  for (i in 1:ncol(data)) dismatrix[i,i]<-1000000
  nearnode<-rep(0,ncol(data))
  nearnodedis<-rep(1000000,ncol(data))
  for (i in 1:ncol(data))
  {
    for (j in 1:ncol(data))
    {
      if (dismatrix[i,j]<nearnodedis[i])
      {
        nearnode[i]<-j
        nearnodedis[i]<-dismatrix[i,j]
      }
    }
  }
  sortdis<-sort(nearnodedis)
  q1<-floor(length(sortdis)/4)
  q3<-length(sortdis)-q1
  maxdis<-sortdis[q3]+1.5*(sortdis[q3]-sortdis[q1])
  removedCells<-array()
  remove=0
  for (i in 1:ncol(data))
  {
    if (nearnodedis[i]>maxdis)
    {
      data<-data[,-i]
      celltype<-celltype[-i]
      removedCells[remove]<-i
      remove=remove+1
    }
  }
  
 # result<-list(data,celltype)
 # result$data<-data
 # result$celltype<-celltype

  return (removedCells)
}


#load the rda file#
#load(file = "D:/Codes/DropOut/data/usoskin.rda")





scSimulator <- function(N=3, nDG=150, nMK=10, nNDG=8180, k=50,
                        seed=17,logmean=5.25,logsd=1, v=9.2){
  set.seed(seed)
  
  scmdSimulator <- function(N,logmean,logsd,nDG){
    trD <- array(0,dim=c(nDG,N))
    for (i in 1:N){
      trD[,i] <- rlnorm(nDG,meanlog = logmean,sdlog= logsd) -1
    }
    return(trD)
  }
  
  ## Simulating Differentially Expressed Genes
  
  DG <- scmdSimulator(N,logmean,logsd,nDG)
  
  ## Simulating Marker Genes
  mkSimulator <- function(N,logmean,logsd,nMK){
    y <- array(0,dim=c(nMK*N,N))
    for (i in 1:N){
      y[(i*nMK-nMK+1):(i*nMK) ,i] <- rlnorm(nMK,meanlog = logmean,sdlog = logsd) -1
    }
    return(y)
  }
  
  MK <- mkSimulator(N,logmean,logsd,nMK)
  
  DG <- rbind(DG,MK)
  
  ## Simulating the non-differentially expressed genes
  NDG0 <-  rlnorm(nNDG,meanlog = logmean,sdlog=logsd) -1
  
  NDG <- array(NA,dim=c(nNDG,N))
  for(i in 1:N){
    NDG[,i] <- NDG0
  }
  
  ## Expected Values Matrix
  EM <- rbind(DG,NDG)
  
  ## Simulating noise and drop-outs
  
  design <- rep(k,N)
  
  u <-0.7
  
  Simulator <- function(EM,design){
    a<-list() ## annotation
    b<-list() ## expression matrix - expected values
    c<-list() ## expression matrix - dropouts simulated
    d<-list() ## tags matrix
    TZ <- nrow(EM)
    for (i in 1:ncol(EM)){
      a[[i]]<- array(0,dim=c(TZ,design[i]))
      b[[i]]<- array(NA,dim=c(TZ,design[i]))
      c[[i]]<- array(NA,dim=c(TZ,design[i]))
      d[[i]]<- array(NA,dim=c(TZ,design[i]))
      
      a[[i]][EM[,i]>0,] <- 1
      
      ## probability function
      pi <- 1/(1 + exp(u * (log2(EM[,i]+1) - v)))
      for (j in 1:design[i]){
        b[[i]][,j] <- EM[,i]
        c[[i]][,j] <- EM[,i]
        for (k in 1:TZ){
          if(rbinom(1,1,pi[k])){
            c[[i]][k,j] <- max(rpois(1,1)-1,0)
            a[[i]][k,j] <- 2 * a[[i]][k,j]
          }
        }
        
        ## simulating noise
        nj<-rpois(TZ,lambda=round(c[[i]][,j]))*(a[[i]][,j]==1)
        nj <- nj+c[[i]][,j]*(a[[i]][,j]==2)
        nj[nj<0]<-0
        d[[i]][,j]<-nj
      }
    }
    
    A <- do.call(cbind,a)
    B <- do.call(cbind,b)
    C <- do.call(cbind,c)
    D <- do.call(cbind,d)
    
    G <-list(A,B,C,D)
    names(G) <-c("annotation","expectedValues","dropoutsSimulated","tags")
    return(G)
  }
  
  sData <- Simulator(EM,design)
  
  return(sData)
}