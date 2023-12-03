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

%% usage:  [MA,AR] = dhr2arma (R,TVP)
%%
%% Calcula el modelo ARMA de los componentes DHR
%%
%% R   > matriz con las raíces AR de cada componente
%% TVP > modelo DHR de cada componente
%%
%% MA  > Matriz cuyas columnas son los polinomio MA de cada componente
%% AR  > Matriz cuyas columnas son los polinomio AR de cada componente

%% Author: Marcos Bujosa <marcos.bujosa@ccee.ucm.es>

function [AR,MA,RAR,RMA,ctes] = dhr2arma (TVP,P,p)

  W=zeros(size(P)); W(P~=Inf)=2*pi./[P(P~=Inf)];
  AR  =zeros(5,size(P,2));  MA  =zeros(3,size(P,2));
  RAR =zeros(4,size(P,2));  RMA =zeros(2,size(P,2));
  ctes=zeros(1,size(P,2));
  %% discriminacion de componentes
  [r1,TN]=find(P==Inf|P==2);	% Tendencia y Nyquist (0 y pi)
  [r2,Ci]=find(P <Inf&P >2);	% Componentes intermedios
  
  re = fix((cos (W(ones(size(TVP,1),1),:)) .* TVP)*10^16)/10^16;
  im = fix((sin (W(ones(size(TVP,1),1),:)) .* TVP)*10^16)/10^16;

%%% frecuencias 0 y pi
  w =   W(1,TN);
  r = TVP(:,TN);
  
  MA(1  ,TN) = 1;
  AR(1:3,TN) = [ones(1,length(TN)); [-sum(r)].*cos(w); prod(r)];
    
  RAR(1:size(TVP,1),TN)=TVP(:,TN).*cos(w(ones(size(TVP,1),1),:));
  ctes(TN)=1;

%%% componentes intermedios 
  w =   W(1,Ci);
  r = TVP(:,Ci);

  AR(:,Ci)=[ones(1,length(Ci));
	    -2*(sum(r)).*(fix(cos(w)*10^16)/10^16);
 	    sum(r.^2)+4*prod(r).*(fix(cos(w)*10^16)/10^16).^2;
 	    -2*( r(1,:).*(r(2,:)).^2 + r(2,:).*(r(1,:)).^2).*(fix(cos(w)*10^16)/10^16);
 	    prod(r.^2)];
  
  %% parte MA para componentes tipo AR(1) o RW 
  a  = @(r,w) sqrt((1+r(1,:).^2 + sqrt((1+r(1,:).^2).^2-4*(r(1,:).*(fix(cos(w)*10^16)/10^16)).^2 ))/2);
  b  = @(r,w) r(1,:).*(fix(cos(w)*10^16)/10^16)./a(r,w);		
  x  = @(r,w) [a(r,w); -b(r,w)];
  
  ci = Ci(find(sum(r>0)==1));	% seleccion componentes tipo AR(1) 
  w  = W(1,ci);   
  rci = TVP(:,ci);
  
  MA(1:2,ci)=x(rci,w);

  RAR(1:2,ci)=[re(1,ci)+im(1,ci)*i;re(1,ci)-im(1,ci)*i];
  RMA(1,ci)=b(rci,w)./a(rci,w);
  ctes(ci)=a(rci,w);

  %% parte MA para componentes tipo AR(2), SRW o IRW   
  a     = @(w)     (fix(cos(2*w)*10^16)/10^16);;
  b     = @(AB,w)  sum(AB).*(fix(cos(w)*10^16)/10^16);
  c     = @(AB)    2+prod(AB);

  sq    = @(AB,w)  sqrt( (b(AB,w)).^2 + 4.*a(w).*(2.*a(w)-c(AB)));
  delta = @(AB,w)  [(b(AB,w)+sq(AB,w))./(2.*a(w));(b(AB,w)-sq(AB,w))./(2.*a(w))];
  gamma = @(AB,w)  [(-delta(AB,w)-sqrt(delta(AB,w).^2-4))./2];

  CI  = Ci(find(sum(r>0)==2)); % seleccion componentes tipo AR(2) 
  w   = W(1,CI);
  rCI = TVP(:,CI);
  AB  = rCI+(1./rCI);

  G = gamma(AB,w);
  G(abs(G)>abs(1./G))=1./G(abs(G)>abs(1./G)); % raices menores que 1

  lambda=sqrt((prod(AB).*(fix(cos(2*w)*10^16)/10^16))./prod(G))/2;
  MA(1:3,CI)=[ones(1,length(CI));-sum(G);prod(G)].*lambda(ones(3,1),:);

  RAR(:,CI)  =[re(1,CI)+im(1,CI)*i;re(1,CI)-im(1,CI)*i;
	       re(2,CI)+im(2,CI)*i;re(2,CI)-im(2,CI)*i];
  RMA(1:2,CI)=G;
  ctes(CI)   =lambda;

  MA(abs(MA  )<10^-12)= 0;
  MA(abs(MA-1)<10^-12)= 1;
  MA(abs(MA+1)<10^-12)=-1;

  AR(abs(AR  )<10^-12)= 0;
  AR(abs(AR-1)<10^-12)= 1;
  AR(abs(AR+1)<10^-12)=-1;

