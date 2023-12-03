%% Copyright (C) 1999, 2000, 2001, 2002, 2003, 2004, 2006, 2008 by Marcos Bujosa
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

%% usage:  [P,TVP,Pac,TVPac,Paux,TVPaux] = iddhr (p,AR,PD,PaP,TVPaP,ar) 
%% 
%% Identifica el modelo DHR de una serie temporal a partir de su modelo AR
%% 
%% p        > Periodicidad de los datos: 1(Anual), 4(trimestral),... 
%% AR       > Polinomio AR del AR-espectro 
%% PD     * > Parámetros de Decisión              ([.05, .45, 36, 250, 2]) 
%%             PD=[1-.05, % Limite para decidir si una raíz es unitaria 
%%		   0.45,  % Limite para decidir si una raíz es cero 
%%		   36,	  % Limite para decidir si una raíz es de T 
%%		   250,	  % Limite para decidir si una raíz es de S 
%%		   2];    % Núm. parámetros de suavidad por comp. 
%% 
%% PaP    * > Periodos a Priori de los componentes DHR (p./(0:floor(p/2))) 
%% TVPap  * > Matriz parámetros suavidad a priori    ([1,...,1;1,...1,0]) 
%% ar     * > Polinomio AR con raices forzadas
%% 
%% P        < Periodos de los componentes DHR 
%% TVP      < Parámetros de los componentes DHR 
%% Pac      < Periodos de los componentes adicionales
%% TVPac    < Parámetros de los componentes adicionales
%% Paux     < Periodos de los componentes adic. (raices múltiples)
%% TVPaux   < Parámetros de los componentes adic. (raices múltiples)
 
%% Author: Marcos Bujosa <marcos.bujosa@ccee.ucm.es>
%% Description:  Identifica el modelo DHR correspondiente a un polinomio
%% AR ajustado a una serie temporal  
 
function [TVP,P] = iddhr (p,Ar,PD,PaP,TVPaP,aR)
                                                     
  %% VALORES POR DEFECTO 
  if nargin<6,aR=[];;end

  PaPaP = [PaP;PaP];
  WaP   = 2*pi./PaPaP; 
  nc    = size(PaP,2);
  
  rDHR         = roots(Ar);
  nDHR         = abs(rDHR);
  precision    = 10^15; 
  wDHR         = fix(precision*angle(rDHR))/precision;  
  pDHR         = inf*ones(size(wDHR));
  pDHR(wDHR~=0)= round(precision*2*pi./wDHR(wDHR~=0))/precision;
  wDHR         = 2*pi./pDHR;

  x    = fix((cos (wDHR) .* nDHR)*10^16)/10^16;
  y    = fix((sin (wDHR) .* nDHR)*10^16)/10^16;
  rDHR = x+y*i;

  [kk, ind] = sort (abs(pDHR),'descend');
  pDHR = pDHR(ind);
  rDHR = rDHR(ind);
  nDHR = nDHR(ind);
  wDHR = wDHR(ind);

  rDHR = rDHR(pDHR>0);
  nDHR = nDHR(pDHR>0);
  wDHR = wDHR(pDHR>0);
  pDHR = pDHR(pDHR>0);

  M    = [pDHR,rDHR,nDHR,wDHR];
  
  if aR
    rAC        = roots(aR);
    nAC        = abs(rAC);
    wAC        = fix(precision*angle(rAC))/precision;  
    pAC        = inf*ones(size(wAC));
    pAC(wAC~=0)= round(precision*2*pi./wAC(wAC~=0))/precision;
    wAC        = 2*pi./pAC;
    
    x   = fix((cos (wAC) .* nAC)*10^16)/10^16;
    y   = fix((sin (wAC) .* nAC)*10^16)/10^16;
    rAC = x+y*i;

    rAC = rAC(pAC>0);
    nAC = nAC(pAC>0);
    wAC = wAC(pAC>0);
    pAC = pAC(pAC>0);
    
    M   = [M;pAC,rAC,nAC,wAC];
  end
  
  periodos = M(:,1); 
  raices   = M(:,2);
  normas   = M(:,3);
  Wj       = M(:,4);
  
  if ~isempty(min(PaP(PaP<Inf & PaP>p))+2*pi/(PD(4)-1))
    %% periodos mayores que p (except. Inf) al ciclo (no a T).
    PD(3) = max(PD(3),min(PaP(PaP<Inf & PaP>p))+2*pi/(PD(4)-1)); % mirar PD(3) en autodhr !!!!
  end

  %% PARÁMETROS DE DECISIÓN 
  lim_uno  = PD(1);	      % Rango -+ PD(1) alrededor de TVP
  lim_cero = PD(2);	      % Limite para decidir si una raíz es cero 
  lim_T    = 2*pi/PD(3);      % Limite para decidir si una raíz es de T 
  lim_S    = 2*pi/PD(4);      % Limite para decidir si una raíz es de S 

  %% Identificación de raices correspondientes a componentes
  frecMOD  = 2*pi./PaP;
  FrecAP   = frecMOD(ones(length(Wj),1),:); % frecuencias a priori
  Frec     = Wj(:,ones(nc,1));		    % frecuencias raices AR
  Peri     = periodos(:,ones(nc,1)); % periodos raices AR
  
  AA       = xor(abs(Frec-FrecAP)<lim_T&Peri>=PD(3)*ones(size(Peri)),...
		 abs(Frec-FrecAP)<lim_S&Peri <PD(3)*ones(size(Peri)));

  for j=1:size(AA,1) % cada raíz solo puede pertenecer a un componente (el primero de la lista)
    AA(j,(min(find(AA(j,:),2))+1):size(AA,2)) = 0;
  end
 
  %%  AA,pause

  [f,c] = find(AA);		%raices AR en 2*pi./PaP 
  
  fdhr  = f(normas(f)>lim_cero);  
  cdhr  = c(normas(f)>lim_cero); % quito pequeñas
  
  %% raices que sobran por TVPaP
  NC       = [1:nc];
  CDHR     = [1:size(cdhr,1)]';
  D        = cdhr(:,ones(nc,1))==NC(ones(length(cdhr),1),:);
  [s,k]    = sort(D,1);   
  [ii,iii] = sort(k);  
  A        = sum(sort(D,1))-sum(TVPaP);                           A(A<0) = 0;
  s(sort(floor(A(ones(1,size(cdhr,1)),:)./CDHR(:,ones(1,size(A,2)))))>0) = 0;

  %%  disp('Aqui esta el problema'),iii,s,iii+vector(ones(1,size(s,1)),:),s(iii+vector(ones(1,size(s,1)),:))
  if size(s,1)==1
    vector = zeros(size(s));
  else
    vector = [0:size(s,1):size(s,1)*(size(s,2)-1)];
  end

  [a,b]    = find(s(iii+vector(ones(1,size(s,1)),:)));
  fdhr     = fdhr(a); 
  cdhr     = cdhr(a);
  m        = sum(cdhr(:,ones(nc,1))==NC(ones(length(cdhr),1),:)); % multiplicidad
  M        = max(2,max(m));

  %% modelo identificado
  TVP      = zeros(M,nc);
  P        = zeros(M,nc);
  W        = zeros(M,nc);
  R        = zeros(M,nc);
  RE       = zeros(M,nc);
  IM       = zeros(M,nc);

  TVP(floor([M:-1:1]'*m/(M))>0) = normas(fdhr);
  P(  floor([M:-1:1]'*m/(M))>0) = periodos(fdhr);
  W(  floor([M:-1:1]'*m/(M))>0) = Wj(fdhr);
  R(  floor([M:-1:1]'*m/(M))>0) = raices(fdhr);

  %% componentes adicionales
  periodos(fdhr) = 0;
  Pac            = periodos(periodos>0)';
  TVPac          = normas(periodos>0)';
  Wac            = Wj(periodos>0)';
  Rac            = raices(periodos>0)';

  %% completando la matriz TVP con nuevas filas de zeros
  %%if M<size(TVPaP,1), TVP(M+1:size(TVPaP,1),:)=zeros(size(TVPaP)-[M 0]); end

  %% FORZANDO EL MODELO
  WS          = [2*pi*[0:floor(p/2)]'/p];
  WaP(TVP==0) = NaN; 
  WaPv        = WaP(:);

  %% raices proximas a TVPaP
  S      = find((abs(TVPaP(:)-TVP(:))<=lim_uno&TVPaP(:)>0)&any(WS(:,ones(1,length(WaPv)))==WaPv(:,ones(1,length(WS)))')');
  TVP(S) = TVPaP(S);
  P(S)   = PaPaP(S);
  W(S)   = WaP(S);

%   RE(S)=cos(2*pi./PaPaP(S)).*TVP(S);
%   IM(S)=sqrt(TVPaP(S).^2-RE(S).^2)*j;
%   R(S)=[RE(S)+IM(S)];

  %% raices proximas a WaP
  WS     = [2*pi*[1:floor(p/2)]'/p];
  %S      = find((abs(WaP-2*pi./P)<=lim_S&TVPaP>0));
  
  %solo modificamos la frecuencia de los estacionales
  S      = (find(any(WS(:,ones(1,length(WaPv)))==WaPv(:,ones(1,length(WS)))')));
  P(S)   = PaPaP(S);

%   W(S)=WaP(S);

%   RE(S)=cos(2*pi./PaPaP(S)).*TVP(S);
%   IM(S)=sqrt(TVPaP(S).^2-RE(S).^2)*j;
%   R(S)=[RE(S)+IM(S)];

  P      = P(1,:);
  TVP    = sort(TVP,'descend');

