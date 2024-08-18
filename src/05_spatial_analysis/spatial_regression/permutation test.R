# This code is written for readability rather than efficiency.

df <- read.csv(" ")
df <- df[,-1]


distfunction <- function(a,b){
  distance <- data.frame(matrix(NA,nrow(a),nrow(b)))
  for (i in 1:nrow(a)) {
    for (j in 1:nrow(b)) {
      distance[i,j] <- sqrt((a$x[i]-b$x[j])^2+(a$y[i]-b$y[j])^2+(a$z[i]-b$z[j])^2)
    }
  }
  
  distance
}

library(dplyr)
library(coin)

difffunction <- function(df){
  aa <- c(1:max(unique(df$cluster)))
  TCdiff <- c()
  MGKdiff <- c()
  
  for (i in aa) {
    a <- df[which(df$cluster==i),]
    minx <- min(a$x)-90
    maxx <- max(a$x)+90
    miny <- min(a$y)-90
    maxy <- max(a$y)+90
    minz <- min(a$z)-90
    maxz <- max(a$z)+90
    
    b <- df[which(df$x>minx & df$x<maxx),]
    b <- b[which(b$y>miny & b$y<maxy),]
    b <- b[which(b$z>minz & b$z<maxz),]
    
    
    b <- b[-which(b$cluster==i),]
    dist <- distfunction(a,b)==45
    out <- b[which(apply(dist, 2, sum)>0),]
    
    
    out$test <- 1
    a$test <- 0
    
    final <- rbind(a,out)
    TCdiff[i] <- mean(final$TC[which(final$cluster==0)])-mean(final$TC[which(final$cluster>0)])
    MGKdiff[i] <- mean(final$MGK[which(final$cluster==0)])-mean(final$MGK[which(final$cluster>0)])
    
  }
  
  testresult <- data.frame(TCdiff,MGKdiff)
}


diff <- difffunction(df)
write.csv(diff,"diff.csv")
