## usage: X = acumula (y,p,t)
##
## agregacion temporal de series temporales 
##
## y     < matriz de series (Txk) 
## p     < acumulacion de p en p periodos
## t     < tipo de agregacion 
##   t=1         : agrega p periodos sumando 
##   t=2         : selecciona el ultimo de cada p periodos 
##   t=3         : selecciona el dato del momento intermedio de cada p
##                  datos (o el inmediatamente posterior)
##   t=4         : selecciona el dato del momento intermedio de cada p
##                  datos (o el inmediatamente anterior)
##   t=5         : selecciona el primero de cada p periodos 
##   t=0 u otros : agrega p periodos promediando
##
## X     > matriz de series acumuladas (Nxk)
##
##
## Copyright (C) 2007, 2008 Marcos Bujosa 

## Author:  MB <marcos.bujosa@ccee.ucm.es>
## Description:  agregacion temporal de series temporales 

function X = acumula(y,p,t)

if nargin < 3; t=0; end

T=size(y,1); N=floor(T/p); n=1; a=[]; b=-1;

if t==1
  sel= ones(1,p); s=0;
elseif t==2
  sel=zeros(1,p); s=0; sel(end)=1;
elseif t==3
  sel=zeros(1,p); s=ceil((p+1)/2); sel(s)=1; b=rem(T,p);
elseif t==4
  sel=zeros(1,p); s=fix((p+1)/2); sel(s)=1; b=rem(T,p);
elseif t==5
  sel=zeros(1,p); s=1; sel(1)=1; b=rem(T,p);
else
  sel= ones(1,p); s=0; n=p;
end

if b >= s; a=p*N+s; end

if isempty(a)
  X=kron(eye(N),sel)*y(1:p*N,:)/n;
else
  X=[kron(eye(N),sel)*y(1:p*N,:)/n;y(a,:)];
end
