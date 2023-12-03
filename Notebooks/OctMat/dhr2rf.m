%% Copyright (C) 2009 by Marcos Bujosa
%% 
%% This program is free software; you can redistribute it and/or modify
%% it under the terms of the GNU General Public License as published by
%% the Free Software Foundation; either version 2 of the License, or
%% (at your option) any later version.
%% 
%% This program is distributed in the hope that it will be useful,
%% but WITHOUT ANY WARRANTY; without even the implied warranty of
%% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%% GNU General Public License for more details.
%% 
%% You should have received a copy of the GNU General Public License
%% along with this program; if not, write to the Free Software
%% Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA

%% usage:  [tr,dr,lr,AR,MA,s2a,Phi,theta,din,lab,CAR,CMA] = dhr2rf (TVP,P,V,nvri,p)
%%
%% nvri=: varianza de innovaciones en la ecuación de obs. si
%% intervención en varianza Q_j=[nvri 0;0 V(j+1)]
%% p=:    periodicidad del   modelo

%% Author: Marcos Bujosa <marcos.bujosa@ccee.ucm.es>

function [tr,dr,lr,AR,MA,s2a,Phi,theta,din,lab,CAR,CMA,RAR,RMA,ctes] = dhr2rf (TVP,P,V,nvri,p)

  if nargin<3;V=[];end; if isempty(V);V=ones(1,size(P,2)+1);end
  if nargin<4;nvri=[];end; if isempty(nvri);nvri=0;end
  %% remove zeros from P, V, and TVP when not all harmonics in the DHR model
  TVP=TVP(:,P>0);V=V([1,find(P)+1]); P=P(P>0); 

  [theta,din,lab,CAR,CMA,RAR,RMA,ctes] = dhr2thd (TVP,P,V,nvri,p);
  
  
  %% MODELO ARMA DEL MODELO DHR COMPLETO
  [Phi, Gam, E, H, D, C, Q, S, R] = thd2ss(theta, din);
  [E,Q,U,iU,P] = sstoinn(Phi, E, H, C, Q, S, R);
  [Phi2, H2, T2, E2, Gam2, F2, Th2, G2] = echelon(size(Phi,1), Phi, H, E, Gam,D);
  AR=F2; MA=Th2; s2a=Q;

  %% CONVERTIR 'CASI' UNOS Y 'CASI' CEROS EN UNOS Y CEROS
  MA(abs(MA  )<10^-12)= 0;
  MA(abs(MA-1)<10^-12)= 1;
  MA(abs(MA+1)<10^-12)=-1;
  
  AR(abs(AR  )<10^-12)= 0;
  AR(abs(AR-1)<10^-12)= 1;
  AR(abs(AR+1)<10^-12)=-1;
  [tr,dr,lr] = arma2thd(F2(2:length(F2)),[],Th2(2:length(Th2)),[],Q,p);
  %%    prtmod(tr,dr,lr),   pause
