%% Copyright (C) 2008 by Marcos Bujosa
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

%% usage:  [gain,phase] = esptarma (MA,AR,n,s2_e)
%%
%% Cálculo de la ganacia y la fase de un filtro lineal (ARMA)
%%
%% MA   > polinomio de la parte de media móvil 
%% AR   > polinomio de la parte autorregresiva
%% n    > frecuencias entre 0 y pi                      % (512)
%%          (si n escalar, w=vec(pi*(0:(n-1))'/n) 
%%           si n vector w=n ) 
%%
%% gain < ganacia del filtro
%% phase< fase del filtro

%% Author: Marcos Bujosa <marcos.bujosa@ccee.ucm.es> 
%% Bibliography: Pollock (1999) 'A handbook of time-series
%% analysis, signal processing and dynamics' pag 466

function [gain,phase,w] = gainandphase (MA,AR,n)

  ARG=['MA  ';
       'AR  ';
       'n   '];

  for a=nargin+1:size(ARG,1)
    eval(sprintf('%s=[];',ARG(a,:)))
  end

  %% Input arguments
  if isempty(MA), error('gainandphase: bad arguments'); end
  if isempty(AR), AR=1; end
  if isempty(n), n=512; end

  %% Frequencies
  if length(n)>1, 
    w=reshape (n, length(n), 1);
  else				% n points in the open interval (0, pi)
    w=reshape (pi*linspace (1/n,(n-1)/n,n), n, 1); 
  end

  W1=w*[0:(length(MA)-1)];
  DR=sum(cos(W1)*diag(MA),2);
  DI=sum(sin(W1)*diag(MA),2);
  
  W2=w*[0:(length(AR)-1)];
  GR=sum(cos(W2)*diag(AR),2);
  GI=sum(sin(W2)*diag(AR),2);
  
  PR=DR.*GR+DI.*GI;
  PI=DI.*GR-DR.*GI;

  n=PR.^2+PI.^2;
  d=GR.^2+GI.^2;

  gain=(sqrt(n)./d);
  phase=-angle(PR+i*PI);	    % en pollock grafico de phase invertido

  %% frequency-response 
  %% fr=(PR+i*PI)./d; 
