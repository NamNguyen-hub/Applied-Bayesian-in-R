function out=loglikx(b1,b2,sigma,y,x1,x2);
v=y-x1*b1-x2*b2;
sterm=0;
isigma=invpd(sigma);
dsigma=logdet(isigma);
T=size(y,1);
N=size(y,2);
for i=1:T
    sterm=sterm+(v(i,:)*isigma*v(i,:)');
end
% out=(-(T*N)/2)*log(2*pi)+(T/2)*dsigma-0.5*sterm;

out=(T/2)*dsigma-0.5*sterm;