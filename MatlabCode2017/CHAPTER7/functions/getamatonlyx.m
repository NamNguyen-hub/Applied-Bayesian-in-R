function [iamat,amat]=getamatonlyx(res,Sbigin,hlast,c0,p00_a,N);
amatx=[];
for j=2:N
%v2=-a1*v1+sqrt(exp(h2))*e2
ytemp=res(:,j);  %v2
xtemp=res(:,1:j-1)*(-1); %-v1
ytemp=ytemp./sqrt(hlast(1:rows(hlast),1).*Sbigin(j));
xtemp=xtemp./repmat(sqrt(hlast(1:rows(hlast),1).*Sbigin(j)),1,cols(xtemp));
%prior means and variances
alpha21_0=c0(j,1:j-1)';
v00_a=diag(abs(alpha21_0))*p00_a;

alpha21_2=getreg(ytemp,xtemp,alpha21_0,v00_a,1); %draw from Normal
amatx=[amatx alpha21_2'];
end
amat=chofac(N,amatx');
iamat=invpd(chofac(N,amatx'));%chofac converts to a lower triangular 