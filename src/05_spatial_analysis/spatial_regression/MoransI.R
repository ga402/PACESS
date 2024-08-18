# This code is written for readability rather than efficiency.
df <- read.csv(file = " ")
colnames(df) <- c("x","y","z","AML","MGK","Tcell")


library(dplyr)
library(magrittr)
Rcpp::sourceCpp("Rcpp_distance_matrix.cpp")
Rcpp::sourceCpp("Rcpp_bikernel_function.cpp")


mat<- as.matrix(df[, c('x', 'y', 'z')])
dmat <- rcpp_distance3d(mat)
rm(mat)
gc()


MoransI_function <- function(a,dmat,df){
  
  dmat1 <- ifelse(dmat==45,1,0)
  
  m <- rep(NA,nrow(df))
  n <- rep(NA,nrow(df))
  n1 <- rep(NA,nrow(df))
  w <- c()
  
  for (i in 1:nrow(df)) {
    for (j in 1:nrow(df)) {
      if(dmat1[i,j]==0){
        m[j] <- 0
        w[j] <- 0
      }else{
        m[j] <- (a[i]-mean(a))*(a[j]-mean(a))
        w[j] <- 1
      }
    }
    n[i] <- sum(m)
    n1[i] <- (a[i]-mean(a))^2
  }
  
  I <- (nrow(df)*(sum(n)))/(sum(dmat1)*sum(n1))
  
  
  I
}


a <- MoransI_function(df$Tcell,dmat,df)
b <- MoransI_function(df$AML,dmat,df)
c <- MoransI_function(df$MGK,dmat,df)

I <- c(a,b,c)
names(I) <- c("Tcell","AML","MGK")
write.csv(I,"MoransI.csv")