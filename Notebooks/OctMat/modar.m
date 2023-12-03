%% usage:  [AR,s2_e,coef] = modar (y,n_retardos,int)
%% 
%% Estimación del Modelo Autorregresivo de una serie temporal
%%
%% y          > serie temporal 
%% n_retardos > orden del modelo AR
%% int        > matriz de datos con intervenciones
%% 
%% AR         > polinomio autorregresivo estimado
%% s2_e       > varianza residual del modelo AR
%% coef       > coeficientes asociados a las intervenciones
%%
%% Copyright (C) 1999, 2000 Marcos Bujosa

%% Author: Marcos Bujosa <marcos.bujosa@ccee.ucm.es>

function [AR,s2_e,coef]=modar(y,n_retardos,int)

  %% VALORES POR DEFECTO
  if nargin<2, error(' [AR, s2_e ] = modar (y, k)'), end
  if nargin<3, int=[]; end

  %% COMPROBACIONES PREVIAS
  [T,Nc]=size(y);
  if min(T,Nc)>1, error('y debe ser un vector'), end
  if T<Nc, y=y'; nc=T; T=Nc; end % garantizamos vector columna
  
  %% ESTIMACIÓN DEL MODELO AUTORREGRESIVO
  
  filas=length(y)-n_retardos;
  Y=y(1+n_retardos:filas+n_retardos); 
  
  X=ones(filas,n_retardos);
  for i=1:n_retardos
    X(:,i)=y(1+n_retardos-i:filas+n_retardos-i);
  end

  if ~isempty(int)
    [fint,nint]=size(int);
    int=int(1:filas);
  end
  
  X=[X,int];

  ar = pinv (X) * Y;		% realizamos la regresión por mco.
  
  r=Y-X*ar;			% calculamos el error
  s2_e=(r'*r)/(filas-n_retardos); % estimamos la varianza residual

  if ~isempty(int)
    coef=ar(n_retardos+1:n_retardos+nint);
  else
    coef=[];
  end
  
  ar=ar(1:n_retardos);
  AR=[1;-ar];

