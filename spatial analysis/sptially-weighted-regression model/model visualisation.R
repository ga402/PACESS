df <- read.csv(file = "~/CD8_45.csv")
df <- df[,-1]

library(dplyr)
library(magrittr)
Rcpp::sourceCpp("~/Rcpp_distance_matrix.cpp")
Rcpp::sourceCpp("~/Rcpp_bikernel_function.cpp")
Rcpp::sourceCpp("~/Rcpp_weighted_coef.cpp")
Rcpp::sourceCpp("~/Rcpp_AICc_calculation.cpp")

X<-data.frame(df[,c("AML","MGK")])
Y<-df$TC %>% as.matrix()
Intercept<-rep(1,nrow(X))
X=cbind(Intercept,X) %>% as.matrix()

mat<- as.matrix(df[, c('x', 'y', 'z')])
dmat <- rcpp_distance3d(mat)




b <- 453

bikernel <- rcpp_bikernel((dmat), b)

out <- rcpp_calculate_coefficient(bikernel, X, Y)

AICc <- rcpp_AICc(X, bikernel, out$residuals %>% as.matrix())


localR <- function(bikernel,Y,hatY) {
  localR <- c()
  for(i in 1:nrow(bikernel)) {
    localR[i] <- 1-(bikernel[i,]%*%((Y-hatY)^2))/(bikernel[i,]%*%((Y-mean(Y))^2))
  }
  return(localR)
}

localR2 <- localR(bikernel,df$TC,out$yfitted)

result <- data.frame(out$coef,out$yfitted,out$residuals,localR2)
colnames(result) <- c("intercept","AMLbeta","MGKbeta","yhat","residuals","localR2")




df1 <- cbind(df,result)


library(reshape2) 
library(ggplot2)
library(grDevices) 
library(RColorBrewer)
library(directlabels) 
require(gridExtra)
colormap <- colorRampPalette(rev(brewer.pal(11,'Spectral')))(10)


a=0

df10 <- df1[which(df1$z==a),]
head(df10)


g1 <-  ggplot(df10,aes(x=x,y=y,fill=AMLbeta))+
  geom_tile()+
  scale_fill_continuous( )+
  scale_fill_gradientn(values = scales::rescale(c(-0.025,0,0.025)),limits=c(-0.025,0.2),breaks=seq(-0.025,0.2,0.05),colors=c("darkgreen","white","#FFB90F","#EEAD0E","#FFBC00","#FF7F00","#EE7600","#CD6600","#8B4500","#8B6914"))+ 
  labs(x="X",y="Y",fill="AML",title = expression(paste("(a) AML at 0 ", mu , "m depth")) ,parse = TRUE)+
  guides(fill=FALSE) +
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  scale_x_continuous(limits = c(400, 3200))+
  scale_y_continuous(limits = c(0, 5900))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(panel.background = element_rect(fill = 'black', colour = 'black'))+
  coord_fixed()   

g2 <-  ggplot(df10,aes(x=x,y=y,fill=(MGKbeta)))+
  geom_tile()+
  scale_fill_continuous( )+
  scale_fill_gradientn(values = scales::rescale(c(-0.25,0,0.25)),limits=c(-0.25,1),breaks=seq(-0.25,1,0.25),colors=c("blue","white","#FF3300","#FF0000","#CC0000","#990000"))+
  labs(x="X",y="Y",fill="MGK",title = expression(paste("(e) MGK at 0 ", mu , "m depth")) ,parse = TRUE)+
  guides(fill=FALSE) +
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  scale_x_continuous(limits = c(400, 3200))+
  scale_y_continuous(limits = c(0, 5900))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(panel.background = element_rect(fill = 'black', colour = 'black'))+
  coord_fixed()   


g3 <-  ggplot(df10,aes(x=x,y=y,fill=localR2))+
  geom_tile()+
  scale_fill_continuous(limits=c(0,1),breaks=seq(0,1,by=0.2))+
  scale_fill_gradientn(limits=c(0,1),colors=colormap,breaks=seq(0,1,by=0.2))+
  labs(x="X",y="Y",fill="localR2",title = expression(paste("(i) local ", R^2, " at 0", mu , "m depth")) ,parse = TRUE)+
  guides(fill=FALSE) +
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  scale_x_continuous(limits = c(400, 3200))+
  scale_y_continuous(limits = c(0, 5900))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(panel.background = element_rect(fill = 'black', colour = 'black'))+
  coord_fixed()   


a=45

df10 <- df1[which(df1$z==a),]
head(df10)


g145 <-  ggplot(df10,aes(x=x,y=y,fill=AMLbeta))+
  geom_tile()+
  scale_fill_continuous( )+
  scale_fill_gradientn(values = scales::rescale(c(-0.025,0,0.025)),limits=c(-0.025,0.2),breaks=seq(-0.025,0.2,0.05),colors=c("darkgreen","white","#FFB90F","#EEAD0E","#FFBC00","#FF7F00","#EE7600","#CD6600","#8B4500","#8B6914"))+ 
  labs(x="X",y="Y",fill="AML",title = expression(paste("(b) AML at 45 ", mu , "m depth")) ,parse = TRUE)+
  guides(fill=FALSE) +
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  scale_x_continuous(limits = c(400, 3200))+
  scale_y_continuous(limits = c(0, 5900))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(panel.background = element_rect(fill = 'black', colour = 'black'))+
  coord_fixed()   


g245 <-  ggplot(df10,aes(x=x,y=y,fill=(MGKbeta)))+
  geom_tile()+
  scale_fill_continuous( )+
  scale_fill_gradientn(values = scales::rescale(c(-0.25,0,0.25)),limits=c(-0.25,1),breaks=seq(-0.25,1,0.25),colors=c("blue","white","#FF3300","#FF0000","#CC0000","#990000"))+
  labs(x="X",y="Y",fill="MGK",title = expression(paste("(f) MGK at 45 ", mu , "m depth")) ,parse = TRUE)+
  guides(fill=FALSE) +
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  scale_x_continuous(limits = c(400, 3200))+
  scale_y_continuous(limits = c(0, 5900))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(panel.background = element_rect(fill = 'black', colour = 'black'))+
  coord_fixed()   


g345 <-  ggplot(df10,aes(x=x,y=y,fill=localR2))+
  geom_tile()+
  scale_fill_continuous(limits=c(0,1),breaks=seq(0,1,by=0.2))+
  scale_fill_gradientn(limits=c(0,1),colors=colormap,breaks=seq(0,1,by=0.2))+
  labs(x="X",y="Y",fill="localR2",title = expression(paste("(j) local ", R^2, " at 45", mu , "m depth")) ,parse = TRUE)+
  guides(fill=FALSE) +
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  scale_x_continuous(limits = c(400, 3200))+
  scale_y_continuous(limits = c(0, 5900))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(panel.background = element_rect(fill = 'black', colour = 'black'))+
  coord_fixed()   



a=90

df10 <- df1[which(df1$z==a),]
head(df10)


g190 <-  ggplot(df10,aes(x=x,y=y,fill=AMLbeta))+
  geom_tile()+
  scale_fill_continuous( )+
  scale_fill_gradientn(values = scales::rescale(c(-0.025,0,0.025)),limits=c(-0.025,0.2),breaks=seq(-0.025,0.2,0.05),colors=c("darkgreen","white","#FFB90F","#EEAD0E","#FFBC00","#FF7F00","#EE7600","#CD6600","#8B4500","#8B6914"))+ 
  labs(x="X",y="Y",fill="AML",title = expression(paste("(c) AML at 90 ", mu , "m depth")) ,parse = TRUE)+
  guides(fill=FALSE) +
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  scale_x_continuous(limits = c(400, 3200))+
  scale_y_continuous(limits = c(0, 5900))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(panel.background = element_rect(fill = 'black', colour = 'black'))+
  coord_fixed()   


g290 <-  ggplot(df10,aes(x=x,y=y,fill=(MGKbeta)))+
  geom_tile()+
  scale_fill_continuous( )+
  scale_fill_gradientn(values = scales::rescale(c(-0.25,0,0.25)),limits=c(-0.25,1),breaks=seq(-0.25,1,0.25),colors=c("blue","white","#FF3300","#FF0000","#CC0000","#990000"))+
  labs(x="X",y="Y",fill="MGK",title = expression(paste("(g) MGK at 90 ", mu , "m depth")) ,parse = TRUE)+
  guides(fill=FALSE) +
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  scale_x_continuous(limits = c(400, 3200))+
  scale_y_continuous(limits = c(0, 5900))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(panel.background = element_rect(fill = 'black', colour = 'black'))+
  coord_fixed()   


g390 <-  ggplot(df10,aes(x=x,y=y,fill=localR2))+
  geom_tile()+
  scale_fill_continuous(limits=c(0,1),breaks=seq(0,1,by=0.2))+
  scale_fill_gradientn(limits=c(0,1),colors=colormap,breaks=seq(0,1,by=0.2))+
  labs(x="X",y="Y",fill="localR2",title =expression(paste("(k) local ", R^2, " at 90", mu , "m depth")) ,parse = TRUE)+
  guides(fill=FALSE) +
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  scale_x_continuous(limits = c(400, 3200))+
  scale_y_continuous(limits = c(0, 5900))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(panel.background = element_rect(fill = 'black', colour = 'black'))+
  coord_fixed()   



a=135

df10 <- df1[which(df1$z==a),]
head(df10)


g1135 <-  ggplot(df10,aes(x=x,y=y,fill=AMLbeta))+
  geom_tile()+
  scale_fill_continuous( )+
  scale_fill_gradientn(values = scales::rescale(c(-0.025,0,0.025)),limits=c(-0.025,0.2),breaks=seq(-0.025,0.2,0.05),colors=c("darkgreen","white","#FFB90F","#EEAD0E","#FFBC00","#FF7F00","#EE7600","#CD6600","#8B4500","#8B6914"))+ 
  labs(x="X",y="Y",fill="AML",title = expression(paste("(d) AML at 135 ", mu , "m depth")) ,parse = TRUE)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  scale_x_continuous(limits = c(400, 3200))+
  scale_y_continuous(limits = c(0, 5900))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(panel.background = element_rect(fill = 'black', colour = 'black'))+
  coord_fixed()   


g2135 <-  ggplot(df10,aes(x=x,y=y,fill=(MGKbeta)))+
  geom_tile()+
  scale_fill_continuous( )+
  scale_fill_gradientn(values = scales::rescale(c(-0.25,0,0.25)),limits=c(-0.25,1),breaks=seq(-0.25,1,0.25),colors=c("blue","white","#FF3300","#FF0000","#CC0000","#990000"))+
  labs(x="X",y="Y",fill="MGK",title = expression(paste("(h) MGK at 135 ", mu , "m depth")) ,parse = TRUE)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  scale_x_continuous(limits = c(400, 3200))+
  scale_y_continuous(limits = c(0, 5900))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(panel.background = element_rect(fill = 'black', colour = 'black'))+
  coord_fixed()   


g3135 <-  ggplot(df10,aes(x=x,y=y,fill=localR2))+
  geom_tile()+
  scale_fill_continuous(limits=c(0,1),breaks=seq(0,1,by=0.2))+
  scale_fill_gradientn(limits=c(0,1),colors=colormap,breaks=seq(0,1,by=0.2))+
  labs(x="X",y="Y",fill=expression(paste("local ", R^2)),title = expression(paste("(m) local ", R^2, " at 135", mu , "m depth")) ,parse = TRUE)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  scale_x_continuous(limits = c(400, 3200))+
  scale_y_continuous(limits = c(0, 5900))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(panel.background = element_rect(fill = 'black', colour = 'black'))+
  coord_fixed()   




grid.arrange(g1,g145,g190,g1135,g2,g245,g290,g2135,g3,g345,g390,g3135,ncol=4)




pdf("~/CD8_45_image_model.pdf", width=12, height=10) 
grid.arrange(g1,g145,g190,g1135,g2,g245,g290,g2135,g3,g345,g390,g3135,ncol=4)
dev.off()

