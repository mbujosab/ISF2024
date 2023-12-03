%% Copyright (C) 2000, 2001, 2003, 2004, 2006, 2008 by Marcos Bujosa 
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

%% -*- texinfo -*-
%% @deftypefn {Function File} {} aresp (@var{y})
%% @deftypefnx {Function File} {} aresp (@var{y}, @var{lag_range})
%% @deftypefnx {Function File} {} aresp (@var{y}, @var{lag_range}, @var{crit})
%% @deftypefnx {Function File} {} aresp (@var{y}, @var{lag_range}, @var{crit}, @var{IA})
%% @deftypefnx {Function File} {[@var{arspt}, @var{s}, @var{lags}, @var{AR}, @var{roots}, @var{norm}, @var{p}] =} aresp (@dots{})
%%
%% Called with a single vector argument @var{y}, calculates the
%% AR-spectrum of the time series
%%
%% The variable @var{lag_range} is the range of AR orders. By
%% default @code{@var{lag_range}=[1 length(@var{y})/3]} (from AR(1)
%% to AR(N/3)).
%%
%% The variable @var{crit} selects the information criterium used to
%% select the AR order. If given @code{@var{crit}=0} (default value)
%% the Akaike is used. If given @code{@var{crit}=1} then Bayesian.
%%
%% If @code{@var{IA}=2} a graph and the roots of the AR polynomial are
%% shown in the screen.
%%
%% Outputs
%%
%% @var{arspt}: AR-spectrum 
%%
%% @var{s}    : residual variance of AR model fitted to the series
%%
%% @table @var
%% @item arspt
%% AR-spectrum
%% @item s
%% residual variance of AR model fitted to the series
%% @item lags
%% order of the AR polynomial
%% @item roots
%% roots of the AR polynomial
%% @item norm
%% norm of the roots of the AR polynomial
%% @item p
%% periods associated to the roots of the AR polynomial
%% @code{2*pi/arg(@var{roots})}
%% @end table
%%
%% @seealso{cia, cib, esptarma}
%%
%% @end deftypefn

%% Marcos Bujosa <marcos.bujosa@ccee.ucm.es>

%% usage:  [arespt,s2_e,n_retardos,AR,raices,norma,periodos] = aresp...
%%    (y,lag_range,crit,IA,nfrec)
%%
%% Estimación del AR-espectro de una serie temporal. 
%%
%% y                > serie temporal
%% lag_range * > [orden mínimo retardos, orden máximo] % ([1 T/3])
%%                     (Si escalar => ese sera el orden del modelo AR)
%%                     (T = longitud de la muestra)
%% crit       * > Criterio de Información               % (0)
%%                     (=0 Akaike, =1 Bayesiano) 
%% IA             * > Información Adicional                 % (0)
%%                     (>0 calcula raíces, =2 información por pantalla) 
%% nfrec          * > frecuencias entre 0 y pi              % (512)
%%                     (numero o vector con las frecuencias) 
%% 
%% arespt           < función AR-espectro
%% s2_e             < varianza residual del modelo AR
%% n_retardos       < numero de retardos empleados en el modelo AR
%% AR               < polinomio autorregresivo estimado
%% raices           < raíces del polinomio AR
%% norma            < norma de las raíces
%% periodos         < periodos asociados a las raíces
%%   (los 3 últimos vectores ordenados según periodos de mayor a menor)
%%
%% Copyright (C) 1999, 2000, 2001, 2008  Marcos Bujosa

%% Author: Marcos Bujosa <marcos.bujosa@ccee.ucm.es>

function [arespt,s2_e,n_retardos,AR,raices,norma,periodos] = aresp (varargin)

% function [arespt,s2_e,n_retardos,AR,raices,norma,periodos]=aresp...
%       (y,lag_range,crit,IA,nfrec)
  
  ARG=['y        ';
       'lag_range';
       'crit     ';
       'IA       ';
       'nfrec    ';];	% nfrec:  vector de frecuencias donde el espectro f(w) es evaluado.
  
  for i=1:nargin
    eval(sprintf('%s=varargin{%f};',ARG(i,:),i))
  end
  
  for i=nargin+1:size(ARG,1)
    eval(sprintf('%s=[];',ARG(i,:)))
  end

  %% VALORES POR DEFECTO
  if isempty(lag_range); lag_range=1:length(y)/3; end
  if isempty(crit), crit=0; end
  if isempty(IA), IA=0; end
  if isempty(nfrec); nfrec=512; end

  %% CÁLCULO DEL AR-ESPECTRO
  if crit==0
    [AR,s2_e,n_retardos]=cia(y,lag_range,IA); % Crit. Inf. Akaike 
  elseif crit==1
    [AR,s2_e,n_retardos]=cib(y,lag_range,IA); % Crit. Inf. Bayesiano 
  end

  if nargin<2 | length(lag_range)==2 
    fprintf('\n\tAR(%1.f)\n\n',n_retardos); 
  end
  
  [arespt,w]=esptarma(1,AR,nfrec); % AR-Espectro.

  %% INFORMACIÓN ADICIONAL
  if IA>0 | nargout > 4		% Cálculo de las raíces
    [periodos,raices,norma] = rootsana (AR);
    if IA==2			% Se muestran las raíces por la pantalla
      fprintf(' \t RAICES \t    NORMA \t\t PERIODOS \n');
      raices_polinomio=[raices,real(norma),real(periodos)]
      figure; semilogy(w,arespt);
      set(gca,'XTick',0:pi/2:pi)
      set(gca,'XTickLabel',{'0','pi/2','pi'})
      set (gca(), "xgrid", "on");
      set(gca,'GridLineStyle','--')
      title(sprintf('f(w): Ar-Espectro de la serie. AR(%g)',n_retardos))
      xlabel('frecuencias w')
      ylabel('f(w) (Escala logarítmica)')
    end
    %pause
  end
 
