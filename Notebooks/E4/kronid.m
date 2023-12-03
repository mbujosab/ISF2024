function [table] = kronid(z,u,ikron,s)
%
% KRONID - Computes an display information criteria to determine
%          the Kronecker indices of the system.
%

table=[];
n=size(z,1);
conterm=(n/2)*(1+log(2*pi));
for i=1:size(ikron,1)
  mcmill=sum(ikron(i,:));
  if (u==[])
    [F,sF,Th,sTh,Q,innov]=varmaeci(z,[],[mcmill+1:mcmill+3],ikron(i,:)',s);
    npar=sum(sum(sF~=0))+sum(sum(sTh~=0));
  else
    [F,sF,Th,sTh,Q,G,sG,innov]=varmaeci(z,u,[mcmill:mcmill+3],ikron(i,:)',s);
    npar=sum(sum(sF~=0))+sum(sum(sTh~=0))+sum(sum(sG~=0));
  end
  ldet=log(det(Q));
  loglik=-(n/2)*ldet-conterm;
  AIC=n*ldet+2*npar;
  SBC=n*ldet+npar*log(n);
  table=[table;i loglik npar AIC SBC ikron(i,:)];
end
[AICS,index] = sort(table(:,4));
disp('Sorted table:')
disp('  logLik   npar      AIC      SBC')
table(index,2:size(table,2))
disp('')
