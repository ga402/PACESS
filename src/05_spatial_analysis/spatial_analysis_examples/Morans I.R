df <- read.csv(file = "~/CD8_45.csv")
df <- df[,-1]



library(dplyr)
library(magrittr)
Rcpp::sourceCpp("~/Rcpp_distance_matrix.cpp")
Rcpp::sourceCpp("~/Rcpp_bikernel_function.cpp")


mat<- as.matrix(df[, c('x', 'y', 'z')])
dmat <- rcpp_distance3d(mat)
rm(mat)
gc()


a <- df$TC

dmat1 <- ifelse(dmat==45,1,0)

m <- rep(NA,nrow(df))
n <- rep(NA,nrow(df))
n1 <- rep(NA,nrow(df))
s0 <- sum(dmat1)
s1 <- c()
s2 <- c()
w <- c()
d1 <- c()
d2 <- c()
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
  s1[i] <- 2*sum(w)
  s2[i] <- (2*sum(dmat1[i,]))^2
  d1[i] <- (a[i]-mean(a))^4
}

I <- (nrow(df)*(sum(n)))/(sum(dmat1)*sum(n1))

S1 <- sum(s1)
S2 <- sum(s2)
A=nrow(df)*((nrow(df)^2-3*nrow(df)+3)*S1-nrow(df)*S2+3*(s0^2))
D=sum(d1)/(sum(n1)^2)
B=D*((nrow(df)^2-nrow(df))*S1-2*nrow(df)*S2+6*(s0^2))
C=(nrow(df)-1)*(nrow(df)-2)*(nrow(df)-3)*(s0^2)

EI2 <- (A-B)/C

EI= (-1)/(nrow(df)-1)
z <- (I-EI)/sqrt(EI2-EI^2)
p <- 2*pnorm(z, lower.tail=FALSE)






a <- df$AML

dmat1 <- ifelse(dmat==45,1,0)

m <- rep(NA,nrow(df))
n <- rep(NA,nrow(df))
n1 <- rep(NA,nrow(df))
s0 <- sum(dmat1)
s1 <- c()
s2 <- c()
w <- c()
d1 <- c()
d2 <- c()
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
  s1[i] <- 2*sum(w)
  s2[i] <- (2*sum(dmat1[i,]))^2
  d1[i] <- (a[i]-mean(a))^4
}

I <- (nrow(df)*(sum(n)))/(sum(dmat1)*sum(n1))

S1 <- sum(s1)
S2 <- sum(s2)
A=nrow(df)*((nrow(df)^2-3*nrow(df)+3)*S1-nrow(df)*S2+3*(s0^2))
D=sum(d1)/(sum(n1)^2)
B=D*((nrow(df)^2-nrow(df))*S1-2*nrow(df)*S2+6*(s0^2))
C=(nrow(df)-1)*(nrow(df)-2)*(nrow(df)-3)*(s0^2)

EI2 <- (A-B)/C

EI= (-1)/(nrow(df)-1)
z <- (I-EI)/sqrt(EI2-EI^2)
p <- 2*pnorm(z, lower.tail=FALSE)

I
z
p






a <- df$MGK

dmat1 <- ifelse(dmat==45,1,0)

m <- rep(NA,nrow(df))
n <- rep(NA,nrow(df))
n1 <- rep(NA,nrow(df))
s0 <- sum(dmat1)
s1 <- c()
s2 <- c()
w <- c()
d1 <- c()
d2 <- c()
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
  s1[i] <- 2*sum(w)
  s2[i] <- (2*sum(dmat1[i,]))^2
  d1[i] <- (a[i]-mean(a))^4
}

I <- (nrow(df)*(sum(n)))/(sum(dmat1)*sum(n1))

S1 <- sum(s1)
S2 <- sum(s2)
A=nrow(df)*((nrow(df)^2-3*nrow(df)+3)*S1-nrow(df)*S2+3*(s0^2))
D=sum(d1)/(sum(n1)^2)
B=D*((nrow(df)^2-nrow(df))*S1-2*nrow(df)*S2+6*(s0^2))
C=(nrow(df)-1)*(nrow(df)-2)*(nrow(df)-3)*(s0^2)

EI2 <- (A-B)/C

EI= (-1)/(nrow(df)-1)
z <- (I-EI)/sqrt(EI2-EI^2)
p <- 2*pnorm(z, lower.tail=FALSE)

I
z
p
