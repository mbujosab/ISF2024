%% Copyright (C) 1999, 2000, 2001 by Marcos Bujosa
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

%% usage:  [fy,w] = esptarma (MA,AR,n,s2_e)
%%
%% Cálculo del espectro de un modelo ARMA
%%
%% MA   > polinomio de la parte de media móvil 
%% AR   > polinomio de la parte autorregresiva
%% n    > frecuencias entre 0 y pi                      % (512)
%%          (si n escalar, w=vec(pi*(0:(n-1))'/n) 
%%           si n vector w=n ) 
%% s2_e > varianza de las innovaciones del proceso ARMA % (1)
%%
%% fy   < espectro
%% w    < el vector de frecuencias donde se ha evaluado fy. 

%% Author: Marcos Bujosa <marcos.bujosa@ccee.ucm.es>

function [fy,w] = esptarma (MA,AR,n,s2_e)

  ARG=['MA  ';
       'AR  ';
       'n   ';
       's2_e';];

  for a=nargin+1:size(ARG,1)
    eval(sprintf('%s=[];',ARG(a,:)));
  end

  %% Input arguments
  if isempty(MA), error('esptarma: bad arguments'); end
  if isempty(AR), AR=1; end
  if isempty(n), n=512; end
  if isempty(s2_e), s2_e=1; end

  %% Frequencies
  if length(n)>1, 
    w=reshape (n, length(n), 1);
  else				% n points in the open interval (0, pi)
    w=reshape (pi*linspace (1/n,(n-1)/n,n), n, 1); 
  end

  %% polynomial fraction Fourier transform
  z=exp(-i*w); E=polyval(MA,z)./polyval(AR,z);

  fy=E.*conj(E)*s2_e;		% Spectrum
