# This code is written for readability rather than efficiency.
df <- read.csv(file = "CD8final.csv")
df <- df[,c(1:4)]
df$AML <- df$AML+1


df$loc <- paste(df$x,",",df$y,",",df$z)
loc <- rep(NA,sum(df$AML))
loc[1:(df$AML[1])] <- df$loc[1]
for (i in 2:nrow(df)) {
  loc[(sum(df$AML[1:i-1])+1):(sum(df$AML[1:i]))]<- rep(df$loc[i],df$AML[i])
}

new <- data.frame(loc)



library(dplyr)
library(tidyr)
#library(tidyverse)
library(reshape2)
library(stringr)

DD <- data.frame(x=as.numeric(str_split_fixed(new$loc,",",3)[,1]),y=as.numeric(str_split_fixed(new$loc,",",3)[,2]),z=as.numeric(str_split_fixed(new$loc,",",3)[,3]))


# Compute DBSCAN using fpc package
set.seed(12345)
db <- fpc::dbscan(DD, eps =45, MinPts = (13+1)*7)
unique(db$cluster)



DD$cluster <- db$cluster

write.csv(DD,file ="AML_DBSCAN1.csv")




rm(list=ls())

df <- read.csv(file = "AML_DBSCAN1.csv")
df <- df[,-1]


df$loc <- paste(df$x,",",df$y,",",df$z)
loc <- unique(df$loc)
l <- c() #cell number
un <- c() #how many cell type in this loc
cluster <- c() 
for (i in 1:length(loc)) {
  a <- df[which(df$loc==loc[i]),]
  l[i] <- nrow(a)
  un[i] <- length(unique(a$cluster))
  cluster[i] <- a$cluster[1]
}


library(dplyr)
library(reshape2)
library(stringr)

DD <- data.frame(x=as.numeric(str_split_fixed(loc,",",3)[,1]),y=as.numeric(str_split_fixed(loc,",",3)[,2]),z=as.numeric(str_split_fixed(loc,",",3)[,3]))

DD$AML <- l
DD$cluster <- cluster


AMLcluster <- c()
for (i in 1:(length(unique(DD$cluster))-1)) {
  a <- DD[which(DD$cluster==i),]
  AMLcluster[i] <- sum(a$AML)
}

DD$AML <- DD$AML-1




df1 <- read.csv(file = "CD8final.csv")
dff <- left_join(df1,DD,by=c("x","y","z","AML"))
write.csv(dff,"AML_DBSCAN2.csv")



rm(list=ls())

DD <- read.csv("AML_DBSCAN2.csv")
DD <- DD[,-1]

distfunction <- function(a,b){
  distance <- data.frame(matrix(NA,nrow(a),nrow(b)))
  for (i in 1:nrow(a)) {
    for (j in 1:nrow(b)) {
      distance[i,j] <- sqrt((a$x[i]-b$x[j])^2+(a$y[i]-b$y[j])^2+(a$z[i]-b$z[j])^2)
    }
  }
  
  distance
}


newclusterfunction <- function(DD,grid){
  
  mergematrix <- data.frame(matrix(0,(length(unique(DD$cluster))-1),(length(unique(DD$cluster))-1)))
  for (i in 1:(length(unique(DD$cluster))-1)) {
    a <- DD[which(DD$cluster==i),]
    
    DDrange <- DD[which(DD$x>=min(a$x)-grid & DD$x<=max(a$x)+grid),]
    DDrange <- DDrange[which(DDrange$y>=min(a$y)-grid & DDrange$y<=max(a$y)+grid),]
    DDrange <- DDrange[which(DDrange$z>=min(a$z)-grid & DDrange$z<=max(a$z)+grid),]
    
    newbid <- c(unique(DDrange$cluster))
    bid <- newbid[-which(newbid==0 | newbid==i)]
    
    if(length(bid)>0){
      for (j in 1:length(bid)) {
        b <- DD[which(DD$cluster==bid[j]),]
        dist <- distfunction(a,b)
        if(any(dist==grid)==TRUE){
          mergematrix[i,bid[j]] <- 1
          mergematrix[bid[j],i] <- 1
        }
      }
    }
  }
  
  DD$newcluster <- DD$cluster
  
  for (i in 1:(length(unique(DD$cluster))-1)) {
    a1 <-  c(which(mergematrix[i,]>0),i)
    
    a <- a1
    if(length(a1)>1){
      for (j in 1:length(a1)) {
        a <- unique(c(a,which(mergematrix[a1[j],]>0)))
      }
    }
    diff <- setdiff(a,a1)
    
    while (length(diff)>0) {
      a1 <- a
      for (k in 1:length(diff)) {
        a <- unique(c(a,which(mergematrix[diff[k],]>0)))
      }
      diff <- setdiff(a,a1)
    }
    DD$newcluster[which(DD$cluster %in% a)] <- min(a)
  }
  
  DD$finalcluster <- DD$newcluster
  for (i in 1:nrow(DD)) {
    A=DD$newcluster
    B=order(unique(DD$newcluster))-1
    C=unique(A)
    DD$finalcluster[i] <- B[which(C==A[i])]
  }
  
  DD
  
}

newDD <- newclusterfunction(DD,45)
newDD <- newDD[,c(1:6,9)]
write.csv(newDD,"AML_DBSCAN3.csv")


rm(list=ls())

df <- read.csv("AML_DBSCAN3.csv")
df <- df[,-1]
length(unique(df$finalcluster))

AMLnum <- c()
for (i in 1:18) {
  AMLnum[i] <- sum(df$AML[which(df$finalcluster==i)])
}
cluster <- data.frame(finalcluster = c(1:18), AMLnum,order=rank(-AMLnum,ties.method = 'random'))

df$cluster <- 0
for (i in 1:18) {
  df$cluster[which(df$finalcluster==i)] <- cluster$order[i]
}


dfz <- df[,c(1:6,8)]
colnames(dfz) <- c("x","y","z","AML","MGK","TC","cluster")
write.csv(dfz,"final.csv")


