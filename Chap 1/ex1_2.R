rm(list=ls())

setwd("/Users/namnguyen/Documents/GitHub/Applied-Bayesian-in-R/Chap 1")
set.seed(123)
# Data set up ####
library(readxl)
inflation <- read_excel("inflation.xls")
colnames(inflation)[1] = "period"

t = nrow(inflation)
lag0 = inflation[3:t,]
lag1 = inflation[2:(t-1),]
lag2 = inflation[1:(t-2),]
x_0 = rep(1, length(lag0))
colnames(lag1)[2] = "lag1"
colnames(lag2)[2] = "lag2"
Y= as.matrix(lag0[,2])
X = as.matrix(cbind(x_0, lag1[,2],lag2[,2]))
t = nrow(X)
#Set priors and starting values ####

## Priors for B 

B0 = (c(0,0,0))
sigma0 = diag(1,3,3)

## Priors for Sigma2 

T0 = 1 
D0 = 0.1

## Starting values

B = B0
sigma2 =1
reps = 15000
burn = 12000
out1 = matrix(0,reps,3)
out2 = c()
out3 = c()
#out1 and out2 are empty matrices that will save the
#draws of B and sigma2 respectively

for(i in 1:reps){
  
  M = solve( solve(sigma0) + drop((1 / sigma2)) * (t(X) %*% X)) %*% ( solve(sigma0) %*% B0 + drop((1 / sigma2)) *( t(X) %*% Y) )
  V = solve(solve(sigma0) + drop((1 / sigma2))*(t(X) %*% X))
  chck= -1
  
  B = M + t((rnorm(3)) %*% chol(V))
  
  b= rbind(c(B[2],B[3]), c(1,0) )
  
  ee = max(abs(eigen(b)$values))
  
  resid = Y - X %*% B
  
  T1 = t + T0
  D1 = D0 + t(resid)%*% resid
  
  # draw from IG
  
  z0 = rnorm(T1)
  z0_z0 = t(z0) %*% z0
  
  sigma2 = D1 /z0_z0
  out1[i,]=t(B)
  out2[i]=sigma2
  
  
} 

# compute forecast for 1 years

yhat=matrix(0,14,reps)
yhat[1:2,] = X[249:250,1]
cfactor = sqrt(sigma2) 
for( i  in 1:reps){
  
  for(m in 3:14){
    yhat[m,i] = (c(1, yhat[(m-1),i], yhat[(m-2),i])) %*% (out1[i,]) +(rnorm(1) * cfactor)
    
  }
}


quants <- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)
EB1=apply( yhat[3:14,] , 1 , quantile , probs = quants , na.rm = TRUE )
EB1

# Plot series

if (!require('reshape2')) install.packages('reshape2'); library('reshape2')
library("ggplot2")
#if(!"devtools" %in% rownames(installed.packages())) install.packages("devtools")
#devtools::install_github("crsh/papaja")
#library("papaja") #APA journal graph template

##Combine series
T_combined = nrow(Y)+12
Y_combined=matrix(NA,T_combined,10)
Y_combined[1:nrow(Y),1]=Y
Y_combined[(nrow(Y)+1):(nrow(Y)+12),2:10]=t(EB1)
##Convert to data frame with date column
dates = seq(as.Date("1948-01-01"), by = "quarter", length.out = 262)
dates = as.character(dates)

Y_combined = cbind(Y_combined,dates)
Y_combined = as.data.frame(Y_combined)


for (i in 1:10){
  Y_combined[,i]=as.numeric(Y_combined[,i])
  Y_combined[nrow(Y),i] = Y_combined[nrow(Y),1]
}

Y_combined[,11]=as.Date(Y_combined[,11])

Y_combined=subset(Y_combined, dates>as.Date("2000-01-01"))
test_data_long <- reshape2::melt(Y_combined, id="dates")  # convert to long format

p <- ggplot(data=test_data_long,
            aes(x=dates, y=value, colour=variable)) +
  geom_line() +
  labs(x = "dates", y = "inflation",
               title = "Forecast of inflation")
p
#p + theme_apa()
