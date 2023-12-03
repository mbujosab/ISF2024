## Copyright (C) 2008 Marcos Bujosa
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

## [ ret ] = tasavariacion (X)
##
## tasavariacion calcula a tasa de variación de las columnas de X

## Author: Marcos Bujosa <marcos.bujosa@ccee.ucm.es>
## 
## 2008-04-09 Marcos Bujosa <marcos.bujosa@ccee.ucm.es>
## * Initial revision

function [ ret ] = tasavariacion (X,n)

  if nargin==1; n=1; endif

  [f,c]=size(X);
  if f==1 & c>f;
    X=X(:);
  endif

  ret=(X(n+1:end,:)-X(1:end-n,:))./X(1:end-n,:);

endfunction
