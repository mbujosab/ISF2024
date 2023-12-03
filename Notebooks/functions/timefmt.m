## Copyright (C) 2000 -- 2009 by Marcos Bujosa 
## 
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2 of the License, or
## (at your option) any later version.
## 
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA

## -*- texinfo -*-
## @deftypefn {Function File} {} timefmt (@var{P}, @var{Y}, @var{M}, @var{N})
## @deftypefnx {Function File} {} timefmt (@var{P}, @var{Y}, @var{M}, @var{N}, @var{e})
##
## Returns date markers for time series data
##
## @var{P} > Periodicity (1 yearly, 4 quaterly, 12 monthly)
##
## @var{Y} > First year, ('1969')
##
## @var{M} > First period, ('3' for march or third quarter)
##
## @var{N} > Number of dates 
##
## @var{e} > month numbers (default value e=0)
##
##  if @var{e}<0 : equidistant points 1969.00  (first quarter), 1969.25  (second quarter)...
## 
##  if @var{e}>=0: e number of zeros before the # of the (month or quarter or ...). If e=1, 1969.01  (first quarter), 1969.02  (second quarter)...
##
## @end deftypefn

## 
## Author: Marcos Bujosa <marcos.bujosa@ccee.ucm.es>

function T = timefmt (P,Y,M,N,e)

  if (nargin < 4 || nargin > 5)
    print_usage ();
  elseif nargin == 4
    e = 0;
  endif

  if M>P,  error (" M>P \n"); endif

  if e<0
    p=[0:1/P:(P-1)/P]';
#   elseif e==2
#     p=[1:P]'./10;
  else
    p=[1:P]'./(10^e);
  endif
    
  p=vec(p(:,ones(1,ceil((N+P)/P)))); # periods

  A=[Y:Y+ceil(N/P)];
  A=vec(A(ones(P,1),:));	# years

  if e<0
    T=A+p;
#   elseif e==2
#     d=sum((10.^[0:9])./P<1)	# decimals
#     T=A+p/(10^d)
  else 
    d=sum((10.^[0:9])./P<1);	# decimals
    T=A+p/(10^d);
  endif

  T=T(M:N+M-1);
  
endfunction