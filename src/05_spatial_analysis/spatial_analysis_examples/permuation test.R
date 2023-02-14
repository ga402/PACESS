df <- read.csv(file = "~/CD8_45_clusteralllevels_eps45_91.csv")
df <- df[,-1]

df$loc <- paste(df$x,",",df$y,",",df$z)
loc <- unique(df$loc)
l <- c() 
un <- c() 
cluster <- c() 
for (i in 1:length(loc)) {
  a <- df[which(df$loc==loc[i]),]
  l[i] <- nrow(a)
  un[i] <- length(unique(a$cluster))
  cluster[i] <- a$cluster[1]
}
table(un)
unique(cluster)

library(dplyr)
library(reshape2)
library(stringr)

DD <- data.frame(x=as.numeric(str_split_fixed(loc,",",3)[,1]),y=as.numeric(str_split_fixed(loc,",",3)[,2]),z=as.numeric(str_split_fixed(loc,",",3)[,3]))

DD$AML <- l
DD$cluster <- cluster
unique(DD$cluster)
head(DD)

AMLcluster <- c()
for (i in 1:(length(unique(DD$cluster))-1)) {
  a <- DD[which(DD$cluster==i),]
  AMLcluster[i] <- sum(a$AML)
}

AMLcluster

AMLcluster/sum(DD$AML)

sum(AMLcluster/sum(DD$AML) >0.01)

df1 <- read.csv(file = "~/CD8_45.csv")
df1 <- df1[,-1]
dff <- left_join(df1,DD,by=c("x","y","z","AML"))

unique(dff$cluster)
table(dff$cluster)

dff$cluster1 <- dff$cluster
dff$cluster1[is.na(dff$cluster)] <- 0

unique(dff$cluster1)
table(dff$cluster1)



PM <- function(treatment, outcome, n){
  original <- diff(tapply(outcome, treatment, mean))
  distribution=c()
  result=0
  for(i in 1:n){
    distribution[i]=diff(by(outcome, sample(treatment, length(treatment), FALSE), mean))
  }
  result=sum(abs(distribution) >= abs(original))/(n)
  return(result)
}




library(dplyr)
library(magrittr)
Rcpp::sourceCpp("~/Rcpp_distance_matrix.cpp")



aa <- c(1:43)
test <- c()
lengthx <- c()
lengthy <- c()
lengthz <- c()
for (i in aa) {
  a <- dff[which(dff$cluster==i),]
  minx <- min(a$x)-90
  maxx <- max(a$x)+90
  miny <- min(a$y)-90
  maxy <- max(a$y)+90
  minz <- min(a$z)-90
  maxz <- max(a$z)+90
  
  b <- dff[which(dff$x>minx & dff$x<maxx),]
  b <- b[which(b$y>miny & b$y<maxy),]
  b <- b[which(b$z>minz & b$z<maxz),]
  
  b1 <- b[which(b$cluster1==i),]
  b2 <- b[which(b$cluster1==0),]
  
  b2cluster <- c()
  for (k in 1:nrow(b2)) {
    bk <- rbind(b2[k,],b1)
    mat<- as.matrix(bk[, c('x', 'y', 'z')])
    dmat <- rcpp_distance3d(mat)
    b2cluster[k] <- sum(dmat[1,]>0&dmat[1,]<90)
  }
  b2change <- b2[which(b2cluster>0),]
  bb <- rbind(b2change,b1)
  bb$test <- 0
  bb$test[which(bb$cluster1!=i)] <- 1
  test[i] <- PM(factor(bb$test),bb$TC,100000) 
}

testresult <- data.frame(table(dff$cluster)[-1])
colnames(testresult) <- c("cluster","#cube")
testresult$cell <-AMLcluster/sum(DD$AML) 
testresult$test <- test
testresult$testTF <- (test<0.05)
testresult$perc <- AMLcluster/sum(DD$AML)


write.csv(testresult,file = "~/CD8_45_clustertestMean.csv")