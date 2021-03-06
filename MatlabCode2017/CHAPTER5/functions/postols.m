function [out,lik]=postols(Y,X,Gammaold,B0,Sigma0,T0,D0)
T=rows(Y);
beta=Gammaold(1:2);
sigma2=Gammaold(3);
resid=Y-X*beta;
    lik=-(T/2)*log(2*pi*sigma2)-0.5*(((resid)'*(resid))/sigma2); %likelihood function
       normalprior=log(mvnpdf(beta,B0,Sigma0)); %evaluate prior for B1 and B2
    gammaprior=gampdf1(T0,D0,1/sigma2); %evaluate prior for 1/sigma
    out=lik+normalprior+gammaprior; 