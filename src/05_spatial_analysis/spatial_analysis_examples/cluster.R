df <- read.csv(file = "~/CD8_45.csv")
df <- df[which(df$AML>0),c(2,3,4,5)]


df$loc <- paste(df$x,",",df$y,",",df$z)
loc <- rep(NA,sum(df$AML))
loc[1] <- df$loc[1]
for (i in 2:nrow(df)) {
  loc[(sum(df$AML[1:i-1])+1):(sum(df$AML[1:i]))]<- rep(df$loc[i],df$AML[i])
}

new <- data.frame(loc)



library(dplyr)
library(reshape2)
library(stringr)
library(fpc)
library(dbscan)

DD <- data.frame(x=as.numeric(str_split_fixed(new$loc,",",3)[,1]),y=as.numeric(str_split_fixed(new$loc,",",3)[,2]),z=as.numeric(str_split_fixed(new$loc,",",3)[,3]))


# Compute DBSCAN using fpc package
set.seed(12345)
db <- fpc::dbscan(DD, eps =45, MinPts = 13*7)
unique(db$cluster)



DD$cluster <- db$cluster


write.csv(DD,file = "~/CD8_45_clusteralllevels_eps45_91.csv")