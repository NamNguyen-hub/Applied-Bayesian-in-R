function [ fload,rhoout,problemout,rho0out,res ] = ...
    getfloadfact( dataS,pmatw,pmatc,FLOAD0,PFLOAD0,rmatin,nfact,NN,...
    rhoin,r0,v0,rho0,index,NC )
LI=cols(rhoin);
T=rows(pmatw);
%draw factor loadings
fload=zeros(NN,nfact*2);
rhoout=zeros(NN,LI);
rho0out=rho0;
problemout=0;
emat=eye(cols(fload));
id=unique(index);
idx=vec(repmat(1:NC,nfact,1));
res=zeros(T,NN);
jx=1;
for ii=1:length(id)
    datax=dataS(:,index==id(ii)); %country i
    pmatx=pmatc(:,idx==id(ii));
    pmat=[pmatw pmatx];
    NNx=cols(datax);
    fload1=zeros(NNx,nfact*2);
rhoout1=zeros(NNx,LI);
rho01=rho0(index==id(ii),:);
rho0out1=rho01;
FLOAD01=FLOAD0(index==id(ii),:);
resx=zeros(T,NNx);
rhoin1=rhoin(index==id(ii),:);
rmatin1=rmatin(:,index==id(ii));
for j=1:nfact*2
    yy=datax(:,j)-pmat*emat(j,:)';
    yys=transformrho(yy,rhoin1(j,:));
         fload1(j,:)=emat(j,:)';
     resx(:,j)=yys;   
         %sample rho
     resid=yy;  %heteroscedasticity and serial correlation
     ye=resid;
     xe=zeros(rows(ye),LI);
    
     for i=1:LI
         xe(:,i)=lag0(ye,i);
     end
     yeh=ye./sqrt(rmatin1(2:end,j));
    xeh=xe./repmat(sqrt(rmatin1(2:end,j)),1,cols(xe));
     
     
     yeh=yeh(LI+1:end,:);
     xeh=xeh(LI+1:end,:);
     [bdraw,problem]=getARx(yeh,xeh,r0,v0,1);
     if problem
         rhoout1(j,:)=rho01(j,:);
     else
          rhoout1(j,:)=bdraw';
          rho0out1(j,:)=bdraw';
     end
     problemout=problemout+problem;
         
end
    
for j=(nfact*2)+1:NNx
    yy=datax(:,j);
    xx=pmat;
    
    %remove serial correlation
    yys=transformrho(yy,rhoin1(j,:));
    xxs=transformrho(xx,rhoin1(j,:));
    
    %remove heteroscedasticity
    yyh=yys./sqrt(rmatin1(2:end,j));
    xxh=xxs./repmat(sqrt(rmatin1(2:end,j)),1,cols(xx));
    
    
    yyh=yyh(LI+1:end,:);
    xxh=xxh(LI+1:end,:);
        
        
        fload0=FLOAD01(j,:)';
        pfload0=PFLOAD0;
   
    BB=getreg(yyh,xxh,fload0,pfload0,1);
   
         fload1(j,:)=BB';
         
         
         %sample rho
     resid=yy-xx*BB;  %heteroscedasticity and serial correlation
     ye=resid;
     xe=zeros(rows(ye),LI);
    
     for i=1:LI
         xe(:,i)=lag0(ye,i);
     end
     yeh=ye./sqrt(rmatin1(2:end,j));
    xeh=xe./repmat(sqrt(rmatin1(2:end,j)),1,cols(xe));
     
     
     yeh=yeh(LI+1:end,:);
     xeh=xeh(LI+1:end,:);
     [bdraw,problem]=getARx(yeh,xeh,r0,v0,1);
     if problem
         rhoout1(j,:)=rho01(j,:);
     else
          rhoout1(j,:)=bdraw';
          rho0out1(j,:)=bdraw';
     end
     problemout=problemout+problem;
    %remove serial correlation
    yyss=transformrho(yy,rhoout1(j,:));
    xxss=transformrho(xx,rhoout1(j,:));  
    resx(:,j)=yyss-xxss*BB;
end
fload(jx:jx+NNx-1,:)=fload1;
rhoout(jx:jx+NNx-1,:)=rhoout1;
rho0out(jx:jx+NNx-1,:)=rho0out1;
res(:,jx:jx+NNx-1)=resx;
jx=jx+NNx;
end