# This code is written for readability rather than efficiency.

df <- read.csv("final.csv")
df <- df[,-1]
head(df)
summary(df)

df$c <- 0
df$c[df$cluster>0] <- 1

df$c <- factor(df$c)

df$TC_q <- 0
df$TC_q[which(df$TC>0 &df$MGK==0)] <- 1

df$MGK_q <- 0
df$MGK_q[which(df$MGK>0 &df$TC==0)] <-1

df$TC_q <- factor(df$TC_q)
df$MGK_q <- factor(df$MGK_q)

df$TC_M <- 0
df$TC_M[which(df$MGK>0 &df$TC>0)] <- 1

df$TC_M <- factor(df$TC_M)

lm <- glm(c ~ TC_q + MGK_q+TC_M, family = binomial, data = df)
summary(lm)
exp(coef(lm))