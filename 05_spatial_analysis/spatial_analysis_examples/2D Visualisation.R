df <- read.csv(file = "~/CD8_45.csv")
df <- df[,-1]
head(df)



library(reshape2) 
library(ggplot2) 
library(grDevices) 
library(RColorBrewer)
library(directlabels) 
require(gridExtra)
colormap <- c("white",colorRampPalette(rev(brewer.pal(11,'Spectral')))(86))


a=0

df10 <- df[which(df$z==a),]
head(df10)


g1 <-  ggplot(df10,aes(x=x,y=y,fill=TC))+
  geom_tile()+
  scale_fill_continuous()+
  scale_fill_gradientn(limits = c(0,15),breaks=seq(0,15,by=3),colors = c("white","#5E4FA2","#5855A5","#535CA8","#4E63AB","#4969AE","#4470B1","#3E77B5","#397DB8","#3484BB","#358BBB","#3B92B8","#4199B5","#479FB3","#4DA6B0","#53ADAD"))+
  labs(x="X",y="Y",fill="T cell",title = expression(paste("T cells at 0 ", mu , "m depth")) ,parse = TRUE)+
  guides(fill=FALSE) + 
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  scale_x_continuous(limits = c(400, 3200))+
  scale_y_continuous(limits = c(0, 5900))+
  theme(plot.title = element_text(hjust = 0.5))+  
  theme(panel.background = element_rect(fill = 'black', colour = 'black'))+
  coord_fixed() 

g2 <-  ggplot(df10,aes(x=x,y=y,fill=(AML)))+
  geom_tile()+
  scale_fill_continuous()+
  scale_fill_gradientn(limits = c(0,86),breaks=seq(0,86,by=10),colors = colormap)+
  labs(x="X",y="Y",fill="AML",title = expression(paste("AML at 0 ", mu , "m depth")) ,parse = TRUE)+
  guides(fill=FALSE) +  
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  scale_x_continuous(limits = c(400, 3200))+
  scale_y_continuous(limits = c(0, 5900))+
  theme(plot.title = element_text(hjust = 0.5))+  
  theme(panel.background = element_rect(fill = 'black', colour = 'black'))+
  coord_fixed() 

g3 <- ggplot(df10,aes(x=x,y=y,fill=(MGK)))+
  geom_tile()+
  scale_fill_continuous()+
  scale_fill_gradientn(limits = c(0,6),breaks=seq(0,6,by=2),colors = c("white","#5E4FA2","#5855A5","#535CA8","#4E63AB","#4969AE","#4470B1"))+
  labs(x="X",y="Y",fill="MGK",title = expression(paste("MGK at 0 ", mu , "m depth")) ,parse = TRUE)+
  guides(fill=FALSE) +   
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  scale_x_continuous(limits = c(400, 3200))+
  scale_y_continuous(limits = c(0, 5900))+
  theme(plot.title = element_text(hjust = 0.5))+  
  theme(panel.background = element_rect(fill = 'black', colour = 'black'))+
  coord_fixed() 



a=45

df10 <- df[which(df$z==a),]
head(df10)


g145 <-  ggplot(df10,aes(x=x,y=y,fill=TC))+
  geom_tile()+
  scale_fill_continuous()+
  scale_fill_gradientn(limits = c(0,15),breaks=seq(0,15,by=3),colors = c("white","#5E4FA2","#5855A5","#535CA8","#4E63AB","#4969AE","#4470B1","#3E77B5","#397DB8","#3484BB","#358BBB","#3B92B8","#4199B5","#479FB3","#4DA6B0","#53ADAD"))+
  labs(x="X",y="Y",fill="T cell",title = expression(paste("T cells at 45 ", mu , "m depth")) ,parse = TRUE)+
  guides(fill=FALSE) + 
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  scale_x_continuous(limits = c(400, 3200))+
  scale_y_continuous(limits = c(0, 5900))+
  theme(plot.title = element_text(hjust = 0.5))+  
  theme(panel.background = element_rect(fill = 'black', colour = 'black'))+
  coord_fixed() 

g245 <-  ggplot(df10,aes(x=x,y=y,fill=(AML)))+
  geom_tile()+
  scale_fill_continuous()+
  scale_fill_gradientn(limits = c(0,86),breaks=seq(0,86,by=10),colors = colormap)+
  labs(x="X",y="Y",fill="AML",title = expression(paste("AML at 45 ", mu , "m depth")) ,parse = TRUE)+
  guides(fill=FALSE) +  
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  scale_x_continuous(limits = c(400, 3200))+
  scale_y_continuous(limits = c(0, 5900))+
  theme(plot.title = element_text(hjust = 0.5))+  
  theme(panel.background = element_rect(fill = 'black', colour = 'black'))+
  coord_fixed() 

g345 <- ggplot(df10,aes(x=x,y=y,fill=(MGK)))+
  geom_tile()+
  scale_fill_continuous()+
  scale_fill_gradientn(limits = c(0,6),breaks=seq(0,6,by=2),colors = c("white","#5E4FA2","#5855A5","#535CA8","#4E63AB","#4969AE","#4470B1"))+
  labs(x="X",y="Y",fill="MGK",title = expression(paste("MGK at 45 ", mu , "m depth")) ,parse = TRUE)+
  guides(fill=FALSE) +   
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  scale_x_continuous(limits = c(400, 3200))+
  scale_y_continuous(limits = c(0, 5900))+
  theme(plot.title = element_text(hjust = 0.5))+  
  theme(panel.background = element_rect(fill = 'black', colour = 'black'))+
  coord_fixed() 




a=90

df10 <- df[which(df$z==a),]
head(df10)


g19 <-  ggplot(df10,aes(x=x,y=y,fill=TC))+
  geom_tile()+
  scale_fill_continuous()+
  scale_fill_gradientn(limits = c(0,15),breaks=seq(0,15,by=3),colors = c("white","#5E4FA2","#5855A5","#535CA8","#4E63AB","#4969AE","#4470B1","#3E77B5","#397DB8","#3484BB","#358BBB","#3B92B8","#4199B5","#479FB3","#4DA6B0","#53ADAD"))+
  labs(x="X",y="Y",fill="T cell",title = expression(paste("T cells at 90 ", mu , "m depth")) ,parse = TRUE)+
  guides(fill=FALSE) + 
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  scale_x_continuous(limits = c(400, 3200))+
  scale_y_continuous(limits = c(0, 5900))+
  theme(plot.title = element_text(hjust = 0.5))+  
  theme(panel.background = element_rect(fill = 'black', colour = 'black'))+
  coord_fixed() 

g29 <-  ggplot(df10,aes(x=x,y=y,fill=(AML)))+
  geom_tile()+
  scale_fill_continuous()+
  scale_fill_gradientn(limits = c(0,86),breaks=seq(0,86,by=10),colors = colormap)+
  labs(x="X",y="Y",fill="AML",title = expression(paste("AML at 90 ", mu , "m depth")) ,parse = TRUE)+
  guides(fill=FALSE) +  
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  scale_x_continuous(limits = c(400, 3200))+
  scale_y_continuous(limits = c(0, 5900))+
  theme(plot.title = element_text(hjust = 0.5))+  
  theme(panel.background = element_rect(fill = 'black', colour = 'black'))+
  coord_fixed() 

g39 <- ggplot(df10,aes(x=x,y=y,fill=(MGK)))+
  geom_tile()+
  scale_fill_continuous()+
  scale_fill_gradientn(limits = c(0,6),breaks=seq(0,6,by=2),colors = c("white","#5E4FA2","#5855A5","#535CA8","#4E63AB","#4969AE","#4470B1"))+
  labs(x="X",y="Y",fill="MGK",title = expression(paste("MGK at 90 ", mu , "m depth")) ,parse = TRUE)+
  guides(fill=FALSE) +   
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  scale_x_continuous(limits = c(400, 3200))+
  scale_y_continuous(limits = c(0, 5900))+
  theme(plot.title = element_text(hjust = 0.5))+  
  theme(panel.background = element_rect(fill = 'black', colour = 'black'))+
  coord_fixed() 





a=135

df10 <- df[which(df$z==a),]
head(df10)


g1135 <-  ggplot(df10,aes(x=x,y=y,fill=TC))+
  geom_tile()+
  scale_fill_continuous()+
  scale_fill_gradientn(limits = c(0,15),breaks=seq(0,15,by=3),colors = c("white","#5E4FA2","#5855A5","#535CA8","#4E63AB","#4969AE","#4470B1","#3E77B5","#397DB8","#3484BB","#358BBB","#3B92B8","#4199B5","#479FB3","#4DA6B0","#53ADAD"))+
  labs(x="X",y="Y",fill="T cell",title = expression(paste("T cells at 135 ", mu , "m depth")) ,parse = TRUE)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  scale_x_continuous(limits = c(400, 3200))+
  scale_y_continuous(limits = c(0, 5900))+
  theme(plot.title = element_text(hjust = 0.5))+  
  theme(panel.background = element_rect(fill = 'black', colour = 'black'))+
  coord_fixed() 

g2135 <-  ggplot(df10,aes(x=x,y=y,fill=(AML)))+
  geom_tile()+
  scale_fill_continuous()+
  scale_fill_gradientn(limits = c(0,86),breaks=seq(0,86,by=10),colors = colormap)+
  labs(x="X",y="Y",fill="AML",title = expression(paste("AML at 135 ", mu , "m depth")) ,parse = TRUE)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  scale_x_continuous(limits = c(400, 3200))+
  scale_y_continuous(limits = c(0, 5900))+
  theme(plot.title = element_text(hjust = 0.5))+  
  theme(panel.background = element_rect(fill = 'black', colour = 'black'))+
  coord_fixed() 

g3135 <- ggplot(df10,aes(x=x,y=y,fill=(MGK)))+
  geom_tile()+
  scale_fill_continuous()+
  scale_fill_gradientn(limits = c(0,6),breaks=seq(0,6,by=2),colors = c("white","#5E4FA2","#5855A5","#535CA8","#4E63AB","#4969AE","#4470B1"))+
  labs(x="X",y="Y",fill="MGK",title = expression(paste("MGK at 135 ", mu , "m depth")) ,parse = TRUE)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  scale_x_continuous(limits = c(400, 3200))+
  scale_y_continuous(limits = c(0, 5900))+
  theme(plot.title = element_text(hjust = 0.5))+  
  theme(panel.background = element_rect(fill = 'black', colour = 'black'))+
  coord_fixed() 







gg <- grid.arrange(g2,g245,g29,g2135,g1,g145,g19,g1135,g3,g345,g39,g3135,ncol=4)




pdf("~/AML.pdf", width=12, height=4)  
grid.arrange(g2,g245,g29,g2135,ncol=4)
dev.off()



pdf("~/TC.pdf", width=12, height=4)  
grid.arrange(g1,g145,g19,g1135,ncol=4)
dev.off()




pdf("~/MGK.pdf", width=12, height=4)  
grid.arrange(g3,g345,g39,g3135,ncol=4)
dev.off()