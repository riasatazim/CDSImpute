LoadUtil()
dataset<-LoadPollenDataSet() # Load Dataset
num <- array( c(1,2,3,4,5,6,7,8,9,10,11)) #Number of Clusters
NOC<-11   #Number of Clusters



dataset<-CDSAlgorithm(dataset)
transposed<-t(dataset)

wss <- (nrow(transposed)-1)*sum(apply(transposed,2,var))
for (i in 2:15) wss[i] <- sum(kmeans(transposed,centers=i)$withinss)
plot(1:15, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")

pcs <- prcomp(transposed, center = F, scale = F)
comp <- data.frame(pcs$x[,1:2])
k <- kmeans(comp, centers = NOC, nstart=25, iter.max=10000)
plot(comp, col=k$clust, pch=16)






#Functions
getPreprocessedData <- function (dataset)
{
  
  dataset<-dataset[complete.cases(dataset), ]
  dataset= apply(dataset[1:nrow(dataset),1:ncol(dataset)],c(1,2),function(x)  log2(1+as.numeric (x)))
  dataset=as.data.frame(dataset)
  dataset<-dataset[complete.cases(dataset), ]
  dataset=dataset[rowSums(dataset == 0) > 10, ]
  removedCells<-SCC_filtration(dataset,num)
  
  for (i in 1:length(removedCells)) 
  {
    dataset<-dataset[,-removedCells[i]]
    
  }
  
  return (dataset)
}


getSimilarIndexes <- function(cellIndex,pearson,spearman,eucledian,pearsonThreshold,spearmanThreshold,eucledianThreshold){
  truthTable_p<-(pearson[,cellIndex] > pearsonThreshold) 
  trueIndexes_p<-which(truthTable_p== TRUE) 
  truthTable_s<-(spearman[,cellIndex] > spearmanThreshold) 
  trueIndexes_s<-which(truthTable_s== TRUE) 
  truthTable_e<-(eucledian[,cellIndex] > eucledianThreshold) 
  trueIndexes_e<-which(truthTable_e== TRUE) 
  
  ## Use Intersection or Union. It is configuarable 
  trueIndexes<-( intersect(trueIndexes_p,trueIndexes_s))
  trueIndexes<-union(trueIndexes,trueIndexes_e)
  return(trueIndexes)
}


getErrorZero<-function(indexes,cellIndex,data,nearzero)
{
  count<-0
  sum<-0
  for (i in indexes )
  {
    if(dataset[cellIndex,i]>nearzero)
    {
      count<-count+1
      sum<-sum+dataset[cellIndex,i]
    }
  }
  
  totalItem<-(count*100)/length(indexes)
  
  # print(totalItem)
  # print("----------")
  if(totalItem > 0.50)
  {
    count=sum/count
    #print(count)
  }
  else{
    count=0
  }
  
  return (count)
}


CDSAlgorithm <- function (dataset)
{
  
  dataset<-getPreprocessedData(dataset)
  transposed<-t(dataset)
  sim1 <- negDistMat(transposed)
  
  m<-min(sim1)
  m<-abs(m)
  sim2<-sim1+m
  maxi<-max(sim2)
  
  eucledian<-sim2/maxi
  
  pearson<-cor(dataset, method = c("pearson"))
  spearman<-cor(dataset, method = c("spearman"))
  
  for (i in 1:nrow(dataset) ) {
    
    print(i)
    
    for (j in 1:ncol(dataset) ) {
      if(dataset[i,j] == 0)
      {
        pearsonThreshold<-quantile(pearson[,j], .95, na.rm=T)
        spearmanThreshold<-quantile(spearman[,j], .95, na.rm=T)
        eucledianThreshold<-quantile(eucledian[,j], .95, na.rm=T)
        
        similarIndexes<-getSimilarIndexes(j,pearson,spearman,eucledian,pearsonThreshold,spearmanThreshold,eucledianThreshold)
        #length(similarIndexes)
        dataset[i,j]<-getErrorZero(similarIndexes,i,dataset,nearzero=1)
        
      }
    }  
  }
  
  return (datset)
  
}


















