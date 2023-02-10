df <- read.csv(file = "~/CD8_45_clusteralllevels_eps45_91.csv")
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




df1 <- read.csv(file = "~/CD8_45.csv")
df1 <- df1[,-1]
dff <- left_join(df1,DD,by=c("x","y","z","AML"))

dff$cluster1 <- dff$cluster
dff$cluster1[is.na(dff$cluster)] <- 0





meanR <- read.csv("~/CD8_45_clustertest.csv")


meanR$rank1 <- rank(meanR$cell,ties.method = "random")

or2 <- c(43:1)
meanR$rank <- 0
for (i in 1:nrow(meanR)) {
  meanR$rank[i] <- which(or2==meanR$rank1[i])
}

dff$clusterrank <- 0
for (i in 1:nrow(dff)) {
  if(dff$cluster1[i]==0){
    dff$clusterrank[i]=0
  }else{
    dff$clusterrank[i]=meanR$rank1[dff$cluster1[i]]
  }
  
}


dff$label <- 0
for (i in 1:nrow(dff)) {
  if(dff$cluster1[i]==0){
    dff$label[i]=0
  }else{
    dff$label[i]=meanR$rank[dff$cluster1[i]]
  }
  
}

library(reshape2) 
library(ggplot2)
library(grDevices) 
library(RColorBrewer)
library(directlabels) 
library(ggrepel)
colormap <- c("white","#CCCCCC","#999999","#666666")



a=0

df10 <- dff[which(dff$z==a),c(1,2,9,10)]
rank <- unique(df10$label[which(df10$label!=0)])
x <- c()
y <- c()
for (i in 1:length(rank)) {
  newdf <- df10[which(df10$label==rank[i]),]
  x[i] <- sort(newdf$x)[ceiling((nrow(newdf))/2)]
  y[i] <- sort(newdf$y)[ceiling((nrow(newdf))/2)]
}
df10$label=0
clusterrank2 <- rep(0,length(rank))
for (i in 1:length(rank)) {
  clusterrank2[i] <- which(or2==rank[i])
}
rankdata <- data.frame(x,y,clusterrank=clusterrank2,label=rank)
rankdata <- rankdata[which(rankdata$label<11),]
newdf10 <- rbind(df10,rankdata)



g1 <- ggplot(newdf10,aes(x=x,y=y,fill=clusterrank))+
  geom_tile()+
  scale_fill_continuous( )+
  scale_fill_gradientn(limits=c(0,43),colors=colormap)+
  labs(x="X",y="Y",fill="AML", title = expression(paste("AML cluster at 0 ", mu , "m depth")) ,parse = TRUE)+
  geom_text_repel(data = subset(newdf10, label > 0),aes(x,y,label = label),fontface="bold", color="red",size=4, min.segment.length = 0, box.padding = 0.3,arrow = arrow(length = unit(0.01, "npc"), type = "open", ends = "last"),segment.color="red",segment.size=0.2,nudge_y=1) +   #,arrow = arrow(length = unit(0.01, "npc"), type = "open", ends = "last"),segment.color="red",segment.size=0.2,nudge_y=1
  guides(fill=FALSE) +
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  scale_x_continuous(limits = c(400, 3200))+
  scale_y_continuous(limits = c(0, 5900))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(panel.background = element_rect(fill = 'black', colour = 'black'))+
  coord_fixed()   





a=45

df10 <- dff[which(dff$z==a),c(1,2,9,10)]
rank <- unique(df10$label[which(df10$label!=0)])
x <- c()
y <- c()
for (i in 1:length(rank)) {
  newdf <- df10[which(df10$label==rank[i]),]
  x[i] <- sort(newdf$x)[ceiling((nrow(newdf))/2)]
  y[i] <- sort(newdf$y)[ceiling((nrow(newdf))/2)]
}
df10$label=0
clusterrank2 <- rep(0,length(rank))
for (i in 1:length(rank)) {
  clusterrank2[i] <- which(or2==rank[i])
}
rankdata <- data.frame(x,y,clusterrank=clusterrank2,label=rank)
rankdata <- rankdata[which(rankdata$label<11),]
newdf10 <- rbind(df10,rankdata)


g2 <- ggplot(newdf10,aes(x=x,y=y,fill=clusterrank))+
  geom_tile()+
  scale_fill_continuous( )+
  scale_fill_gradientn(limits=c(0,43),colors=colormap)+
  labs(x="X",y="Y",fill="AML", title = expression(paste("AML cluster at 45 ", mu , "m depth")) ,parse = TRUE)+
  geom_text_repel(data = subset(newdf10, label > 0),aes(x,y,label = label),fontface="bold", color="red",size=4, min.segment.length = 0, box.padding = 0.3,arrow = arrow(length = unit(0.01, "npc"), type = "open", ends = "last"),segment.color="red",segment.size=0.2,nudge_y=1 )+ 
  guides(fill=FALSE) +
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  scale_x_continuous(limits = c(400, 3200))+
  scale_y_continuous(limits = c(0, 5900))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(panel.background = element_rect(fill = 'black', colour = 'black'))+
  coord_fixed()   


a=90

df10 <- dff[which(dff$z==a),c(1,2,9,10)]
rank <- unique(df10$label[which(df10$label!=0)])
x <- c()
y <- c()
for (i in 1:length(rank)) {
  newdf <- df10[which(df10$label==rank[i]),]
  x[i] <- sort(newdf$x)[ceiling((nrow(newdf))/2)]
  y[i] <- sort(newdf$y)[ceiling((nrow(newdf))/2)]
}
df10$label=0
clusterrank2 <- rep(0,length(rank))
for (i in 1:length(rank)) {
  clusterrank2[i] <- which(or2==rank[i])
}
rankdata <- data.frame(x,y,clusterrank=clusterrank2,label=rank) 
rankdata <- rankdata[which(rankdata$label<11),]
newdf10 <- rbind(df10,rankdata)




g3 <- ggplot(newdf10,aes(x=x,y=y,fill=clusterrank))+
  geom_tile()+
  scale_fill_continuous( )+
  scale_fill_gradientn(limits=c(0,43),colors=colormap)+
  labs(x="X",y="Y",fill="AML", title = expression(paste("AML cluster at 90 ", mu , "m depth")) ,parse = TRUE)+
  geom_text_repel(data = subset(newdf10, label > 0),aes(x,y,label = label),fontface="bold", color="red",size=4, min.segment.length = 0, box.padding = 0.3,arrow = arrow(length = unit(0.01, "npc"), type = "open", ends = "last"),segment.color="red",segment.size=0.2,nudge_y=1)+
  guides(fill=FALSE) +
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  scale_x_continuous(limits = c(400, 3200))+
  scale_y_continuous(limits = c(0, 5900))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(panel.background = element_rect(fill = 'black', colour = 'black'))+
  coord_fixed()   


a=135

df10 <- dff[which(dff$z==a),c(1,2,9,10)]
rank <- unique(df10$label[which(df10$label!=0)])
x <- c()
y <- c()
for (i in 1:length(rank)) {
  newdf <- df10[which(df10$label==rank[i]),]
  x[i] <- sort(newdf$x)[ceiling((nrow(newdf))/2)]
  y[i] <- sort(newdf$y)[ceiling((nrow(newdf))/2)]
}
df10$label=0
clusterrank2 <- rep(0,length(rank))
for (i in 1:length(rank)) {
  clusterrank2[i] <- which(or2==rank[i])
}
rankdata <- data.frame(x,y,clusterrank=clusterrank2,label=rank) 
rankdata <- rankdata[which(rankdata$label<11),]
newdf10 <- rbind(df10,rankdata)



g4 <- ggplot(newdf10,aes(x=x,y=y,fill=clusterrank))+
  geom_tile()+
  scale_fill_continuous( )+
  scale_fill_gradientn(limits=c(0,43),colors=colormap)+
  labs(x="X",y="Y",fill="AML", title = expression(paste("AML cluster at 135 ", mu , "m depth")) ,parse = TRUE)+
  geom_text_repel(data = subset(newdf10, label > 0),aes(x,y,label = label),fontface="bold", color="red",size=4, min.segment.length = 0, box.padding = 0.3,arrow = arrow(length = unit(0.01, "npc"), type = "open", ends = "last"),segment.color="red",segment.size=0.2,nudge_y=1)+ 
  guides(fill=FALSE) +
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  scale_x_continuous(limits = c(400, 3200))+
  scale_y_continuous(limits = c(0, 5900))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(panel.background = element_rect(fill = 'black', colour = 'black'))+
  coord_fixed()   

grid.arrange(g1, g2, g3,g4,ncol=4)






pdf("~/light_dark_grey_red10.pdf", width=12, height=4)  
grid.arrange(g1, g2, g3,g4,ncol=4)
dev.off()



