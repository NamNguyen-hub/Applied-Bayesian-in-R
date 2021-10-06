## Setting data source ----
rm(list=ls())
setwd("D:/GitHub/Applied-Bayesian-in-R/Chap 1")
set.seed(123)

# generate artificial data
T = 100
X = cbind(matrix(1,T,1), rnorm(T))
btrue=c(1, 0.5)
sigmatrue = 0.2
Y = X%*%btrue+rnorm(T)*sqrt(sigmatrue)

#set priors
T0=3
D0=2.5

B0=matrix(0,2,1)
Sigma0=diag(4,2,2)

#analytical computation of the marginal likelihood
mlika=mlikols(B0, Sigma0)