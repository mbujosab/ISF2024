%% Copyright (C) 2004, 2006 by Marcos Bujosa
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

%% usage:  [theta,din,lab,AR,MA,RAR,RMA,ctes] = dhr2thd (TVP,P,V,nvri,p)
%%
%% nvri=: varianza de innovaciones en la ecuación de obs. si
%% intervención en varianza Q_j=[nvri 0;0 V(j+1)]
%% p=:    periodicidad del   modelo

%% Author: Marcos Bujosa <marcos.bujosa@ccee.ucm.es>

function [theta,din,lab,CAR,CMA,RAR,RMA,ctes] = dhr2thd (TVP,P,V,nvri,p)

  %% Falta intervencion en varianza (mirar en cdhr2thd,m)
  
  if nargin<3;V=[];end; if isempty(V);V=ones(1,size(P,2)+1);end
  if nargin<4;nvri=[];end; if isempty(nvri);nvri=0;end
  %% remove zeros from P, V, and TVP when not all harmonics in the DHR model
  TVP=TVP(:,P>0);V=V([1,find(P)+1]); P=P(P>0); 

  [CAR,CMA,RAR,RMA,ctes]=dhr2arma(TVP,P);
  
  [t,d,l]=arma2thd([],[],[],[],V(1),p); % COMPONENTE IRREGULAR

  for j=[1:length(P)] 		% COMPONENTES DHR
    MA=CMA(:,j)';AR=CAR(:,j)'; MA=MA(1:max(find(MA))); AR=AR(1:max(find(AR)));
    
    V(j+1)=V(j+1)*MA(1)^2; MA=MA./MA(1); MA(1)=[]; AR(1)=[];

    [tj,dj,lj]=arma2thd([AR],[],[MA],[],V(j+1),p);

    %tj=[tj ones(size(tj))]; tj(tj(:,1)==V(j+1),2)=0;
    %tj=[tj,zeros(size(tj))]; tj(1:(size(AR,2)),2)=1;
    
    %% seguramente intervencion en varianza es meter nvri en alguna
    %% matriz de las de la representación ss de aqui abajo
    %%[Phi, Gam, E, H, D, C, Q, S, R] = thd2ss(tj, dj), pause   
    [t,d,l]=stackthd(t,d,tj,dj,l,lj);
  end
  [theta,din,lab]=comp2thd(t,d,l);  
  %% prtmod(theta,din,lab),pause
