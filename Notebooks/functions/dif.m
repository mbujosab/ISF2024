## Copyright (C) 1999, 2002 by Marcos Bujosa
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

## usage:  y = dif (x,retardo) 
## 
## Operador diferencia (1-B) ó (1-B^r), donde r=retardo 
## dif(X) nos da la serie diferenciada: 
## Y(t)=X(t)-X(t-1) 
## dif(X,r) nos da la diferencia respecto al retardo r: 
## Y(t)=X(t)-X(t-r) si r es positivo (RETARDO) 
## Y(t)=X(t)-X(t+r) si r es es negativo (ADELANTO) 
##
## Si X es una matriz, dif actua sobre cada columna

## Author: Marcos Bujosa <marcos.bujosa@ccee.ucm.es>

function y = dif (x,retardo) 
   
  if nargin<2, retardo=1; end
  [f,c]=size(x);
  if f==1 & c>f;
    x=x(:);
  end
  T=size(x,1);
  
  if retardo==0;
    y=zeros(size(x));		         # no hay retardo 
  elseif retardo>0; 
    y=x(1+retardo:T,:)-x(1:T-retardo,:); # retardo r 
  elseif retardo<0;
    y=x(1:T+retardo,:)-x(1-retardo:T,:); # adelanto r 
  end 

  