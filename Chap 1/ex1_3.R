## Section 0 ----

rm(list=ls())

setwd("D:/GitHub/Applied-Bayesian-in-R/Chap 1")
# set.seed(123)
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

#Set priors and starting values 

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
reps = 15000 #samples to be saved for inference
burn = 12000 #throwing off
out1 = matrix(0,14,reps)
out2 = matrix(0,reps,3)
out3 = c()
out4 = c()
xstar = matrix(0,(nrow(X)-1),3)

xlag0 = as.matrix(X[2:t,])
xlag1 = as.matrix(X[1:(t-1),])
## Loop
for(i in 1:reps){
  
  ystar = Y-lag1[,2]*rho
  ystar = ystar[2:nrow(ystar),]
  xstar = xlag0 - xlag1*drop(rho)
  
  M = solve( solve(sigma0) + drop((1 / sigma2)) * (t(xstar) %*% xstar)) %*% ( solve(sigma0)   %*% B0 + drop((1 / sigma2)) *( t(xstar) %*% ystar) )
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
  
  #Test code for matlab
  # B =(c(    0.1508,  1.4647,  -0.5267))
  # Step 3 compute rho
  y = Y - X %*% B
  x = y[1:(nrow(y)-1),]
  y = y[2:nrow(y),]
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
  
  
  ## Test code for matlab
  # rho = -0.0834
  #step 3 sample sigma2 conditional on B from IG(T1,D1);
  #compute residuals
  resid = ystar - xstar %*% B
  
  T1 = t + T0
  D1 = D0 + t(resid)%*% resid
  # draw from IG
  
  z0 = rnorm(T1)
  
  ##Matlab
  # z0 =c(   -0.5936,    1.5303,    1.7322,    0.9999,   -1.2749,   -0.3055,   -1.4697,    0.3927,   -1.5489,    0.1936,   -1.6989,   -0.1246,    0.0902,   -1.4038,    3.6343,   -1.2333,   -1.0590,    0.3080,    0.2599,   -1.6080,    1.0916,    0.7229,   -0.8250,   -0.6928,    0.2588,   -0.2062,    1.7994,    0.6527,    1.1839,    0.3350,    1.2304,    0.4781,    1.4572,    1.2094,    0.8147,    0.3630,    0.7009,   -0.8755,   -0.0458,    0.1601,   -0.0296,   -1.0226,   -0.7903,   -0.4719,    0.8646,   -0.1233,   -1.4399,    0.3701,    1.2072,    0.9195,    1.6516,   -0.1547,   -1.2086,   -0.8235,   -0.5219,   -1.4479,   -0.3988,    1.3256,   -0.4551,    0.8785,   -0.0248,    0.6788,   -1.1214,   -0.3492,    1.5018,    0.7030,   -0.8644,    0.8592,   -1.4859,    1.1707,    0.6713,    1.1859,    0.3922,   -2.3568,   -1.4334,   -0.3466,    0.5229,    1.2816,    1.8146,    1.5016,   -1.0764,   -0.9949,   -1.4605,   -0.0106,    0.0233,    0.1593,    0.8595,   -0.8514,   -0.2743,    0.6674,    0.3800,   -0.9014,    0.0653,   -0.2986,   -1.5941,   -0.1463,   -0.2101,   -0.0166,   -0.5062,    0.2641,   -1.4320,    0.1459,   -0.1308,   -0.1970,   -0.3011,   -1.1426,   -0.9057,   -1.2300,    2.0263,    0.1335,   -1.3039,   -0.7322,    0.2548,    1.4804,    0.0208,    1.0943,    0.1111,    0.8229,    1.3053,    1.3733,   -1.5081,    0.6432,    2.1557,    0.1422,   -1.5035,   -0.5482,   -1.6047,    1.1985,   -1.2870,   -0.7872,   -1.1791,    0.6178,   -0.7486,    0.9578,   -0.2219,    0.5970,   -0.2036,    1.6074,   -1.5897,   -0.9965,    1.0965,    0.3077,    1.1045,   -0.3867,    0.4425,   -0.9712,   -0.2855,    1.2405,    0.0615,   -0.6303,    0.1215,   -0.1435,    0.6775,   -0.4554,   -1.4972,    1.5450,   -1.2301,   -1.1834,    0.0288,   -0.0482,   -0.5631,    0.8091,    2.1604,   -0.4100,    1.4811,   -0.1718,   -0.2636,    1.1911,   -0.2346,   -1.4435,   -1.4490,   -2.1132,    0.6347,    1.7541,   -2.0682,    0.0135,    0.2552,   -0.2383,   -1.2675,   -2.0064,    0.8377,   -2.2665,   -1.2533,    0.5703,   -2.6337,   -0.4275,    1.1147,    0.5720,   -0.5201,   -0.5637,    0.0201,   -2.3718,    1.4571,    0.6926,    0.1547,    2.4355,    0.3196,   -1.3583,   -0.9750,    0.0394,    0.5055,    1.4310,   -0.8004,   -0.0890,   -2.2015,   -0.7436,    0.8837,    0.0259,    0.1712,   -0.0216,   -1.5249,    2.1933,    0.3211,    0.9549,   -0.7606,    0.7286,    0.2614,    0.9762,    0.1412,   -0.0969,    0.5340,    0.6383,   -0.8954,   -0.8952,    0.0886,    0.3825,   -1.1989,   -0.1951,    0.3406,    0.3959,    0.0415,    0.7060,    3.0576,    0.9091,    1.2968,    1.0235,    0.6306,   -0.1708,    0.5785,   -0.0832,    0.2158,   -0.5434,   -0.2287,    0.9761,    1.1513,    0.3855,   -0.9997,   -1.9844,    0.3681,    1.0568,    1.1833)
  
  ## End matlab
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

## 2nd section ----

#plot marginal posterior distributions
par(mfrow=c(2,3))
hist(out2[burn:reps,1], main="Posterior distribution of intercept",xlab="constant",col="blue")
hist(out2[burn:reps,2], main="Posterior distribution of AR(1)",xlab="ar(1)",col="blue")
hist(out2[burn:reps,3], main="Posterior distribution of AR(2)",xlab="ar(2)",col="blue")
hist(out4[burn:reps], main="Posterior distribution of sigma2",xlab="Variance",col="blue")

hist(out3[burn:reps], main="Posterior distribution of rho",xlab="rho",col="blue")


#compute mean of the marginal posterior distribution of intercept
# MB=colMeans(out2[burn:reps,])
# MB
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
# # #compute mean of the marginal posterior distribution of rho
# # MB=colMeans(out3[burn:reps])
# # #Compute standard error
# # VB=apply(out3[burn:reps],2,sd)
# # #compute percentile
# # quants <- c(0.05,0.95)
# # EB=apply( out3[burn:reps] , 2 , quantile , probs = quants , na.rm = TRUE )
# # 
# # MB
# # VB
# # EB
# # 
# # #compute mean of the marginal posterior distribution of sigma
# # MB=colMeans(out4[burn:reps])
# # #Compute standard error
# # VB=apply(out4[burn:reps],2,sd)
# # #compute percentile
# # quants <- c(0.05,0.95)
# # EB=apply( out4[burn:reps] , 2 , quantile , probs = quants , na.rm = TRUE )
# # 
# # MB
# # VB
# # EB

#Forecast
quants <- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)
EB1=apply( out1[3:14,burn:reps] , 1 , quantile , probs = quants , na.rm = TRUE )



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
# library("papaja") #APA journal graph template

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
test_data_long <- reshape2::melt(Y_combined, id="dates")  # convert to long format

p <- ggplot(data=test_data_long,
            aes(x=dates, y=value, colour=variable)) +
  geom_line() +
  labs(x = "dates", y = "inflation",
       title = "Forecast of inflation")
p
# p + theme_apa()