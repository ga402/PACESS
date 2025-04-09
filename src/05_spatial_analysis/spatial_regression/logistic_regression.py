
# This code is written for readability rather than efficiency.
# 
# - It is aimed to give a set of logical steps explaining how to run the analysis
# in the manner we ran it for this study. 
# 




# Read and inspect data
df <- read.csv("your_file.csv")
df <- df[ , -1]  # Drop the first column (e.g., index)
head(df)
summary(df)

# Binary outcome: c = 1 if cluster > 0, else 0
df$c <- ifelse(df$cluster > 0, 1, 0)
df$c <- factor(df$c)

# Define TC_q: 1 if TC > 0 and MGK == 0
df$TC_q <- ifelse(df$TC > 0 & df$MGK == 0, 1, 0)
df$TC_q <- factor(df$TC_q)

# Define MGK_q: 1 if MGK > 0 and TC == 0
df$MGK_q <- ifelse(df$MGK > 0 & df$TC == 0, 1, 0)
df$MGK_q <- factor(df$MGK_q)

# Define TC_M: 1 if both MGK > 0 and TC > 0
df$TC_M <- ifelse(df$MGK > 0 & df$TC > 0, 1, 0)
df$TC_M <- factor(df$TC_M)

# Logistic regression
model <- glm(c ~ TC_q + MGK_q + TC_M, family = binomial, data = df)
summary(model)

# Exponentiated coefficients (odds ratios)
exp(coef(model))
