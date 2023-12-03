function [yhat,dethat,pval,thopt,iter,fnew,g,h]=outliers(thopt,din,y,tol,nodif)

% outliers -  Automatic detection of outliers and decomposition of a time
%             series into deterministic and stochastic components
%    [yhat,dethat,pval,thopt,it,fval,g,h]=outliers(thopt,din,y,tol)
% theta  > parameter vector.
% din    > matrix which stores a description of the model dynamics.
% y      > matrix of observable variables.
% tol    > tolerance for intervention, defined as the number of standard errors
%          required to consider any residual an outlier.
% nodif  > if nodif=0, intervention is restricted to impulse-type effects, if nodif=1 it
%          allows for both, impulse and step-type effects
% yhat   < estimate of the stochastic component of y
% dethat < estimate of the (additive) deterministic component of y, such that:
%          y=yhat + dethat
% pval   < percent probability (risk) of undue intervention
% thopt  > optimal parameter vector of the model for yhat
% iter   < number of iterations.
% fnew   < value of the objective function at the optimum.
% g      < gradient at the optimum.
% h      < hessian (or numerical approximation) at the optimum.
%
% 7/3/97

% Copyright (C) 1997 Jaime Terceiro
% 
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2, or (at your option)
% any later version.
% 
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details. 
% 
% You should have received a copy of the GNU General Public License
% along with this file.  If not, write to the Free Software Foundation,
% 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.

global E4OPTION

if nargin < 4,  e4error(3); end
if nargin == 4, nodif=1; end
z=transdif(y,1,nodif);
din=touser(din,'differ');

pval=200*(1-normcdf(tol,0,1));

numatip=1;
while numatip ~= 0
  [thopt,iter,fnew,g,h] = e4min('lfmiss', thopt,'', din,z);
  [ehat,vt]=residual(thopt,din,z);
  nanes = find(abs(vt./std(vt))>tol);
  numatip=size(nanes,1);
  if numatip > 0
    z(abs(vt./std(vt))>tol,1)=NaN*ones(numatip,1);
  end
end

[zhat,pz,xhat,px]=fismiss(thopt,din,z);

yhat=zeros(size(y)); yhat(1,:)=y(1,:);

for i=2:size(y,1)
   yhat(i,:)=yhat(i-1,:)+zhat(i-1,:);
end

dethat=y-yhat;
