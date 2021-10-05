rm(list=ls())

setwd("D:/GitHub/Applied-Bayesian-in-R/Chap 1")
set.seed(123)
# Data set up ####
library(readxl)
inflation <- read_excel("inflation.xls")
colnames(inflation)[1] = "period"

t = nrow(inflation)
lag0 = inflation[3:t,]
lag1 = inflation[2:(t-1),]
lag2 = inflation[1:(t-2),]
X_0 = rep(1, length(lag0))
colnames(lag1)[2] = "lag1"
colnames(lag2)[2] = "lag2"
Y= as.matrix(lag0[,2])
X = as.matrix(cbind(X_0, lag1[,2],lag2[,2]))
t = nrow(X)

#Set priors and starting values ####

## Priors for B 
#MODEL IS Y=ALPHA+B1*Y(T-1)+B2*Y(T-2)+ET
B0 = (c(0,0,0))
sigma0 = diag(1,3,3)

## Priors for Sigma2 

T0 = 1 
D0 = 0.1

## Priors for rho
rho0 = 0
Sigma0r =1

## Starting values

B = B0
rho=rho0
sigma2 =1
reps = 500 #samples to be saved for inference
burn = 0 #throwing off
out1 = matrix(0,14,reps)
out2 = matrix(0,reps,3)
out3 = c()
out4 = c()
xstar = X[1:(nrow(X)-1),]
# Loop ####
for(i in 1:reps){
  
  ystar = Y-lag1[,2]*rho
  ystar = ystar[2:nrow(ystar),]
  xstar[,1] = X[1:(nrow(X)-1),1]-X[2:nrow(X),1]*rho
  xstar[,2] = X[1:(nrow(X)-1),2]-X[2:nrow(X),2]*rho
  xstar[,3] = X[1:(nrow(X)-1),3]-X[2:nrow(X),3]*rho
  
  M = solve( solve(sigma0) + drop((1 / sigma2)) * (t(xstar) %*% xstar)) %*% ( solve(sigma0) %*% B0 + drop((1 / sigma2)) *( t(xstar) %*% ystar) )
  V = solve(solve(sigma0) + drop((1 / sigma2))*(t(xstar) %*% xstar))

  chck= -1
  while (chck < 0) {
    B = M + t((rnorm(3)) %*% chol(V))
    b= rbind(c(B[2],B[3]), c(1,0) )
    ee = max(abs(eigen(b)$values))
    if (ee<=1){
      chck=1
    }
  }
  
  # Step 3 compute rho
  y = Y - X %*% B
  x = y[2:nrow(y),]
  y = y[1:(nrow(y)-1),]
  MM = solve( solve(Sigma0r) + drop((1 / sigma2)) * (t(x) %*% x)) %*% ( solve(Sigma0r) %*% rho0 + drop((1 / sigma2)) *( t(x) %*% y) )
  VV = solve( solve(Sigma0r) + drop((1 / sigma2)) * (t(x) %*% x)) 
  
  
  chck= -1
  while (chck < 0) {
    rho = MM + t(rnorm(1)%*% chol(VV))
    ee = abs(rho)
    if (ee<=1){
      chck=1
    }
  }
  
  #step 3 sample sigma2 conditional on B from IG(T1,D1);
  #compute residuals
  resid = ystar - xstar %*% B
  
  T1 = t + T0
  D1 = D0 + t(resid)%*% resid
  # draw from IG
  
  z0 = rnorm(T1)
  z0_z0 = t(z0) %*% z0
  
  sigma2 = D1 /z0_z0
  
  cfactor = sqrt(sigma2) 
  yhat=matrix(0,14,reps)
  vhat=matrix(0,14,reps)
  
  out2[i,]=t(B)
  out3[i]=rho
  out4[i]=sigma2
  if( i  > burn){

    yhat[1:2,] = Y[249:250,1]
    for(m in 3:14){
      vhat[m,i] = vhat[(m-1),i]%*%rho+rnorm(1)*cfactor
      yhat[m,i] = (c(1, yhat[(m-1),i], yhat[(m-2),i])) %*% (out2[i,]) +vhat[m,i]
      out1[m,i] = yhat[m,i] 
    }
    
  }
}

#plot marginal posterior distributions
par(mfrow=c(1,1))
hist(out2[burn:reps,1], main="Posterior distribution of intercept",xlab="constant",col="blue")

hist(out2[burn:reps,2], main="Posterior distribution of AR(1)",xlab="ar(1)",col="blue")
hist(out2[burn:reps,3], main="Posterior distribution of AR(2)",xlab="ar(2)",col="blue")
hist(out4[burn:reps], main="Posterior distribution of sigma2",xlab="Variance",col="blue")

hist(out3[burn:reps], main="Posterior distribution of rho",xlab="rho",col="blue")

# 
# #compute mean of the marginal posterior distribution of intercept
# MB=colMeans(out2[burn:reps,1])
# #Compute standard error
# VB=apply(out2[burn:reps,1],2,sd)
# #compute percentile
# quants <- c(0.05,0.95)
# EB=apply( out2[burn:reps,1] , 2 , quantile , probs = quants , na.rm = TRUE )
# 
# MB
# VB
# EB
# 
# #compute mean of the marginal posterior distribution of AR(1)
# MB=colMeans(out2[burn:reps,2])
# #Compute standard error
# VB=apply(out2[burn:reps,2],2,sd)
# #compute percentile
# quants <- c(0.05,0.95)
# EB=apply( out2[burn:reps,2] , 2 , quantile , probs = quants , na.rm = TRUE )
# 
# MB
# VB
# EB
# 
# #compute mean of the marginal posterior distribution of AR(2)
# MB=colMeans(out2[burn:reps,3])
# #Compute standard error
# VB=apply(out2[burn:reps,3],2,sd)
# #compute percentile
# quants <- c(0.05,0.95)
# EB=apply( out2[burn:reps,3] , 2 , quantile , probs = quants , na.rm = TRUE )
# 
# MB
# VB
# EB
# 
# #compute mean of the marginal posterior distribution of rho
# MB=colMeans(out3[burn:reps])
# #Compute standard error
# VB=apply(out3[burn:reps],2,sd)
# #compute percentile
# quants <- c(0.05,0.95)
# EB=apply( out3[burn:reps] , 2 , quantile , probs = quants , na.rm = TRUE )
# 
# MB
# VB
# EB
# 
# #compute mean of the marginal posterior distribution of sigma
# MB=colMeans(out4[burn:reps])
# #Compute standard error
# VB=apply(out4[burn:reps],2,sd)
# #compute percentile
# quants <- c(0.05,0.95)
# EB=apply( out4[burn:reps] , 2 , quantile , probs = quants , na.rm = TRUE )
# 
# MB
# VB
# EB

#Forecast
quants <- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)
EB1=apply( out1[3:14,burn:reps] , 1 , quantile , probs = quants , na.rm = TRUE )
EB1


# Plot series
### Rearrange and combine
T_combined = nrow(Y)+12
Y_combined=matrix(NA,T_combined,10)
Y_combined[1:nrow(Y),1]=Y
Y_combined[(nrow(Y)+1):(nrow(Y)+12),2:10]=t(EB1)
##Convert to time series
##Set date and variable

# Example codes: https://www.statmethods.net/advstats/timeseries.html
# # from Jan 2009 to Dec 2014 as a time series object
# myts <- ts(myvector, start=c(2009, 1), end=c(2014, 12), frequency=12)
# 
# # subset the time series (June 2014 to December 2014)
# myts2 <- window(myts, start=c(2014, 6), end=c(2014, 12))
# 
# # plot series
# plot(myts)

dates = seq(as.Date("1948-01-01"), by = "quarter", length.out = 262)
dates = as.character(dates)

Y_combined = cbind(Y_combined,dates)
Y_combined = as.data.frame(Y_combined)
myts <- ts(Y_combined, start=c(1948, 1), frequency=12)

# Plot series
library("reshape2")
library("ggplot2")
#if(!"devtools" %in% rownames(installed.packages())) install.packages("devtools")
#devtools::install_github("crsh/papaja")
library("papaja") #APA journal graph template

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
}

Y_combined[,11]=as.Date(Y_combined[,11])

Y_combined=subset(Y_combined, dates>as.Date("2000-01-01"))
test_data_long <- melt(Y_combined, id="dates")  # convert to long format

p <- ggplot(data=test_data_long,
            aes(x=dates, y=value, colour=variable)) +
  geom_line() +
  labs(x = "dates", y = "inflation",
       title = "Forecast of inflation")
p
p + theme_apa()