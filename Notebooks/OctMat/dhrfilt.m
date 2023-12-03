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

%% usage:  [trend,season,cycle,irreg,xhat,px,e] = dhrfilt ...
%%    (y,P,TVP,V,p,filter,g,terminal,output,missing) 

%% Author: Marcos Bujosa <marcos.bujosa@ccee.ucm.es>
 
function [trend,season,cycle,irreg,thetat,dint,xhat,px,e] = dhrfilt ...
      (y,P,TVP,V,p,filter,g,missing,int_var,exog,coef) 
  
  ARG=['y       ';
       'P       ';
       'TVP     ';
       'V       ';
       'p       ';
       'filter  ';
       'g       ';
       'missing ';
       'int_var ';
       'exog    ';
       'coef    ';];

  for i=nargin+1:size(ARG,1)
    eval(sprintf('%s=[];',ARG(i,:)))
  end

  %% VALORES POR DEFECTO

  if isempty(filter),filter=0;end % filtrado (defecto e4trend con suavizado)
  if isempty(g), g=0; end	% Por defecto no muestra los gráficos
  
  terminal = [];output = [];

  global E4OPTION		% hay que inicializar la toolbox e4
  if isempty(E4OPTION) || exist ('E4OPTION')==0
    e4init
  end
  
  while size(TVP,2)<size(P,2)
    TVP=[TVP,TVP(:,size(TVP,2))];
  end

  Pseason=P(P<=p);
  Pcycle= P(P>p & P<Inf);

  if ~isempty(int_var)
    V2=V;
    %V2([3:8])=100000;
    [t2,d2,l2] = dhr2thd (TVP,P,V2,10000,p);
  end
  
  [t,d,l] = dhr2thd (TVP,P,V,[],p);

  if filter==1 %e4trend (toinnov=1, to obtain exact estimates)
    [trend,season,cycle,irreg,thetat,dint,ixmodes,xhat]=e4trend(t,d,y,1);
    %%   hh=(dif(trend));  %% curiosidad!!
    %%   cte=irreg(1)./hh(1);
    %%   plot([irreg(1:end-1),(dif(trend))*cte])  
    if any(cycle)
      trend=[trend+cycle,trend,cycle];
      titulo_trend=sprintf('Tendencia + ciclo');
    else
      titulo_trend=sprintf('Tendencia');
    end
  elseif filter==2
    [Phi, Gam, E, H, D, C, Q, S, R] = thd2ss(t,d);
    D=coef; Gam=zeros(size(Phi,1),size(D,2));
    [t,d,l]=ss2thd(Phi, Gam, E, H, D, C, Q, S, R);
    
    if any (isnan(y)) 
      if isempty(int_var)
	[zhat, Pz, xhat, px] = fismiss(t, d, [y,exog]);
	e=y-xhat*H';
	e(find(isnan (y)))=0;
      else
	[zhat, Pz, xhat, px] = fismiss2(t, d, [y,exog], t2, d2,int_var);
	e=y-xhat*H';
	e(find(isnan (y)))=0;       
      end
    else
      if isempty(int_var)
	[xhat, px, e] = fismod(t, d, [y,exog]);
      else
	[xhat, px, e] = fismod2(t, d, [y,exog], t2, d2,int_var);
      end
    end
    
    C=xhat(:,[find(H)]);
    
    sel=find(H);
    if ~isempty(find(P(find(TVP(1,:)))==inf))
      trend=xhat(:,[sel(find(P(find(TVP(1,:)))==inf))]);
    else
      trend=zeros(size(y)); 
    end
    titulo_trend=sprintf('Tendencia');
    
    if ~isempty(find(P(find(TVP(1,:)))>p+1e-13 & P(find(TVP(1,:)))<inf))
      cycle=xhat(:,[sel(find(P(find(TVP(1,:)))>p+1e-13 & P(find(TVP(1,:)))<inf))]);
      if size(cycle,2)>1; cycle=[sum(cycle')',cycle]; end
    else
      cycle=zeros(size(y));
    end
    
    if ~isempty(find(P(find(TVP(1,:)))<=p+1e-13))
      season=xhat(:,[(sel(find(P(find(TVP(1,:)))<=p+1e-13)))]);
      if size(season,2)>1; season=[sum(season')',season]; end
    else
      season=zeros(size(y));
    end
    irreg=e;    
  else %e4trend (toinnov=0, smoothed estimates)
    [trend,season,cycle,irreg,thetat,dint,ixmodes,xhat]=e4trend(t,d,y,0);
    if any(cycle)
      trend=[trend+cycle,trend,cycle];
      titulo_trend=sprintf('Tendencia + ciclo');
    else
      titulo_trend=sprintf('Tendencia');
    end
  end


%   trend(trend==0)=NaN;
%   trend(trend==0)=NaN;
%   trend(trend==0)=NaN;
%   trend(trend==0)=NaN;


  if g==1
    num_fig=[];%figure+1;
    figure;plot(irreg);title('Componente irregular');

    if any(cycle)
      if size(cycle,2)>2; 
	for i=2:size(cycle,2)
	  if isempty(Pcycle); Pcycle='                 '; end
	  figure;plot(cycle(:,i));title(sprintf('Componente del ciclo %g',Pcycle(i-1)));
	end 
	figure;plot(cycle(:,1));title('Componente de ciclo completo');
      else 
	figure;plot(cycle);title('Componente de ciclo');
      end
    end
    
    if any(season)
      if size(season,2)>2; 
	for i=2:size(season,2)
	  if isempty(Pseason); Pseason='                 '; end
	  figure;plot(season(:,i));title(sprintf('Componente estacional de periodo %g',Pseason(i-1)));
	end 
	figure;plot(season(:,i));title('Componente estacional completo');
      else 
	figure;plot(season),title('Componente estacional');
      end
    end

     if any(trend)
       figure;plot([trend(:,1),y]);title(titulo_trend);legend('Tendencia','Serie    ');
     end
   end
   drawnow ();
