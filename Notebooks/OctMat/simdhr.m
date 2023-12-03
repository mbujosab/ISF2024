%% usage: [y,C] = simdhr (n,p,TVP,V,P)
%%
%% n=:    sample size
%% p=:    sampling frequency or samples per unit 
%%        (Examples: 1 for annual data; 4 for quarterly; 12 for monthly, etc.)

%% Copyright (C) 2004-2014 by Marcos Bujosa
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

%% Author: Marcos Bujosa <mbujosab@ucm.es>
%% Keywords: linear dynamic harmonic regression, dhr, ldhr 

function [y,C] = simdhr (T,p,TVP,V,P)

  global E4OPTION		% Initializes the global variables in the toolbox4
  if isempty(E4OPTION) || exist ('E4OPTION')==0
    e4init;
  end
  
  if nargin<5; P=p./[0:p/2];end;

  if size(TVP,2)<size(P,2)
    TVP = [TVP, TVP(:,size(TVP,2))*ones(1,size(P,2)-size(TVP,2))];
  end

  if size(V,2)<size(P,2)
    V   = [V, V(:,size(V,2))*ones(1,size(P,2)-size(V,2)+1)];
  end

  if nargout<2
    [theta,din,lab,CAR,CMA,RAR,RMA,ctes] = dhr2thd (TVP,P,V,[],p);
    s = simmod(theta,din,T+1000); y=s(1001:end);
  else
    [CAR,CMA,RAR,RMA,ctes]=dhr2arma(TVP,P);
    
    [t,d,l]=arma2thd([],[],[],[],V(1),p); % IRREGULAR
    irr = simmod(t,d,T);
    
    for j=[1:length(P)] 		  % DHR COMPONENTS
      MA=CMA(:,j)';AR=CAR(:,j)'; MA=MA(1:max(find(MA))); AR=AR(1:max(find(AR)));
      
      V(j+1)=V(j+1)*MA(1)^2; MA=MA./MA(1); MA(1)=[]; AR(1)=[];
      
      [tj,dj,lj]=arma2thd([AR],[],[MA],[],V(j+1),p);
      s =  simmod(tj,dj,T+1000);
      S(:,j) = s(1001:end);
    end
    C = [irr,S];
    y = sum(C,2);
  end

%!demo
%! n=150; p=12; TVP=[1 1;1 0]; V=[500, 1, 5]; P=p./[0:p/2];
%! [y,C] = simdhr (n,p,TVP,V,P);
%! figure(1), plot (y); figure(2), plot(C)
%! # you should now see the series and its components
