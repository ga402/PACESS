df <- read.csv(file = "~/CD8_45.csv")
df <- df[,-1]
head(df)



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




library(lava)
min <- 0
max <- 6054

brange <- seq(min,max,by=1)
AICc <- rep(NA,length(brange))

for (j in 1:length(brange)) {
  
  b <- brange[j]
  
  dmat <- rcpp_distance3d(mat)
  
  bikernel <- rcpp_bikernel((dmat), b)
  
  rm(dmat)
  gc()
  
  out    <- try(rcpp_calculate_coefficient(bikernel, X, Y),silent=FALSE)
  if('try-error' %in% class(out))          
  {
    AICc[j] <- NA                               
  }else{
    AICc[j]<-rcpp_AICc(X, bikernel, out$residuals %>% as.matrix())
  }
  
  rm(bikernel)
  rm(out)
  gc()
}


AICcdf <- data.frame(b=brange,AICc=AICc)
library(ggplot2)
ggplot(AICcdf,aes(x=b,y=AICc))+geom_point()+geom_line()




write.csv(AICcdf,"~/CD8AICc_45.csv")


