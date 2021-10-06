## Setting data source ----
rm(list=ls())
setwd("D:/GitHub/Applied-Bayesian-in-R/Chap 1")
set.seed(123)
source("./functions/lag0.R") # custom function to create lag variable

## Generate argificial data 
T = 200
Y = matrix(0,T,1)
b0 = c(0.1, 1.1, -0.9)
r0 = 0.8
e = rnorm(T)
v = matrix(0,T,1)

for(i in 3:T){
  v[i]=v[i-1]*r0+e[i]
  Y[i]=b0[1]+t(as.matrix(c(Y[i-1], Y[i-2])))%*%(b0[2:3])+v[i]
}

X = cbind(matrix(1,T,1), lag0(Y,1), lag0(Y,2))
#remove missing obs
Y = as.matrix(Y[3:T])
X = as.matrix(X[3:T,])

#step 1 set priors and starting values
#priors for B
B0 = (c(0,0,0))
Sigma0=diag(1,3,3)
#prior for sigma2
T0=1
D0=0.1
#prior for rho
rho0=0
Sigma0r=1

#starting values
B=B0
rho=rho0
sigma2=1

#parameters for the Geweke diagnostic
Pa=0.1
Pb=0.5

reps = 12500
burn = 11500
out1=c()
out2 = matrix(0,reps,3)
out3=c()


## Loop ----
for(i in 1:reps){
  #step 2 Sample B conditional on sigma N(M*,V*)
  #remove serial correlation
  ystar=Y-lag0(Y,1)*drop(rho)
  xstar=X-lag0(X,1)*drop(rho)
  ystar=ystar[2:nrow(ystar)]
  xstar=xstar[2:nrow(xstar),]
  M = solve(solve(Sigma0) + drop((1 / sigma2)) * (t(xstar) %*% xstar)) %*% ( solve(Sigma0)   %*% B0 + drop((1 / sigma2)) *( t(xstar) %*% ystar) )
  V = solve(solve(Sigma0) + drop((1 / sigma2))*(t(xstar) %*% xstar))
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
  x = y[1:(nrow(y)-1),]
  y = y[2:nrow(y),]
  MM = solve( solve(Sigma0r) + drop((1 / sigma2)) * (t(x) %*% x)) %*% ( solve(Sigma0r) %*% rho0 + drop((1 / sigma2)) *( t(x) %*% y) )
  VV = solve( solve(Sigma0r) + drop((1 / sigma2)) * (t(x) %*% x)) 
  #draw rho but again ensure stationarity
  chck= -1
  while (chck < 0) {
    rho = MM + t(rnorm(1)%*% chol(VV))
    ee = abs(rho)
    if (ee<=1){
      chck=1
    }
  }
  
  #compute residuals
  resid = ystar - xstar %*% B
  #compute posterior df and scale matrix
  T1 = T + T0
  D1 = D0 + t(resid)%*% resid
  # draw from Inverse Gamme
  z0 = rnorm(T1)
  z0_z0 = t(z0) %*% z0
  sigma2 = D1 /z0_z0
  
  out1[i]=sigma2
  out2[i,]=t(B)
  out3[i]=rho
}

out1=as.matrix(out1[(burn+1):reps])
out2=as.matrix(out2[(burn+1):reps])
out3=as.matrix(out3[(burn+1):reps])


## Compute statistics ----
#Compute Geweke diagnostics for convergence
Na = ceiling(Pa*nrow(out1)) #first sample
Nb = ceiling(Pb*nrow(out1)) #second sample
Nstar = nrow(out1)-Nb

#Calculation for sigma
Ma = mean(out1[1:Na]) #mean of first sub-sample
Sa = spectrum(out1[1:Na],freq=0)$bandwidth #variance of first sub-sample
Mb = mean(out1[(Nstar+1):nrow(out1)]) #mean of second sub-sample
Sb = spectrum(out1[(Nstar+1):nrow(out1)],freq=0)$bandwidth #variance of second sub-sample
CD1 = (Ma-Mb)/sqrt((Sa/Na)+(Sb/Nb))
RNE1 = (sd(out1))^2/spectrum(out1,freq=0)$bandwidth

#Calculation for elements of coefficient matrix
CD2=c()
RNE2=c()
for (j in 1:ncol(out2)){
  Ma = mean(out2[1:Na,j]) #mean of first sub-sample
  Sa = spectrum(out2[1:Na,j],freq=0)$bandwidth #variance of first sub-sample
  Mb = mean(out2[(Nstar+1):nrow(out1),j]) #mean of second sub-sample
  Sb = spectrum(out2[(Nstar+1):nrow(out1),j],freq=0)$bandwidth #variance of second sub-sample
  CD= (Ma-Mb)/sqrt((Sa/Na)+(Sb/Nb))
  CD2 = cbind(CD2, CD)
  RNE2 = cbind(RNE2,(sd(out2[,j]))^2/spectrum(out2[,j],freq=0)$bandwidth)
}

#Calculation for rho
Ma = mean(out3[1:Na]) #mean of first sub-sample
Sa = spectrum(out3[1:Na],freq=0)$bandwidth #variance of first sub-sample
Mb = mean(out3[(Nstar+1):nrow(out1)]) #mean of second sub-sample
Sb = spectrum(out3[(Nstar+1):nrow(out1)],freq=0)$bandwidth #variance of second sub-sample
CD1 = (Ma-Mb)/sqrt((Sa/Na)+(Sb/Nb))
RNE1 = (sd(out3))^2/spectrum(out3,freq=0)$bandwidth
