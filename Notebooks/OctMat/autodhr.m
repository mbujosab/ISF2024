%% Copyright (C) 1999 - 2009 by Marcos Bujosa 
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

%% usage: [VAR,P,TVP,oar,MCNN,NVR] = autodhr ...
%%      (y,p,PD,RO,PaP,TVPaP,fijo) 
%%
%% Identifica un modelo DHR y estima las varianzas de las innovaciones
%% de los componentes
%% 
%% y       > Serie temporal
%% p     * > Periodicidad de los datos: 1(Anual), 4(trimestral), ... % (1)
%% PD    * > Parametros de Decision (vector con 5 valores) 
%%            PD=[ lim_1, lim_0, lim_T, lim_S, num_param_suavidad]
%%              lim_1:  Limite para decidir si una raiz es como TVP
%%                       Si raiz\in TVP\pm lim_1 => raiz=TVP       % (0.05)
%%              lim_0:  Limite para decidir si una raiz es cero
%%                       Si |raiz|>lim_0 => raiz modulo 0          % (0.45)
%%              lim_T:  Limite para decidir si una raiz es de T
%%                       Si  Pj   >lim_T => raiz a Tendencia       %  (p*3)
%%              lim_Wj: Limite para decidir si una raiz es de S
%%                       Si |Pcomponente-Pj| < lim_S => raiz a S   %  (500)
%%              nps   : Num. parametros de suavidad: alpha_j, beta_j,...
%%                       excepto para T y Nyquist maximo 2         %    (2)
%% RO    * > Rango Ordenes para los polinomios AR, RO=[minimo,maximo] 
%%              Si p=1 =>  RO=[8,20] % ( p+2 , min(floor(length(y)/4),3*p))
%% PaP   * > Periodos a Priori de componentes DHR   %   (p./(0:floor(p/2)))
%% TVPap * > Matriz parametros suavidad a priori    % ([1,1,1...;1,1,1...])
%% fijo  * > Si fijo=1 estima el modelo: PaP + TVPaP% (0)
%%
%% VAR        < Varianzas de las innovaciones de los componentes
%% P          < Periodos de los componentes DHR
%% TVP        < Parametros de los componentes DHR
%% oar        < Orden del polinomio AR elegido
%% model      < estructura con los modelos identificados y estimados
%% MCNN       < Tipo de estimacion empleado: 0 MCO; 1 MCNN
%% NVR        < Hyperparametros NVR del modelo

%% Author: Marcos Bujosa <marcos.bujosa@ccee.ucm.es>
%% Description:  Identifica y estima el modelo DHR de una serie temporal

%function [VAR,P,TVP,oar,,results,MCNN,NVR] =  autodhr (y,p,PD,RO,PaP,TVPaP,fijo,exog,missing)

function [VAR,P,TVP,oar,results,MCNN,NVR] =  autodhr (varargin)
  
  if nargin < 2,
    error('[VAR,P,TVP,oar] = autodhr(y,p,PD,RO,PaP,TVPaP,fijo,exog,missing)');
  end
  
  global E4OPTION		% hay que inicializar la toolbox e4
  if isempty(E4OPTION) || exist ('E4OPTION')==0
    e4init
  end
    
  ARG=['y      ';
       'p      ';
       'PD     ';
       'RO     ';
       'PaP    ';
       'TVPaP  ';
       'fijo   ';
       'exog   ';
       'missing'];
  
  for i=1:size(ARG,1)
    eval(sprintf('%s=[];',ARG(i,:)));
  end
  
  for i=1:nargin
    if ~isempty(varargin{i}); eval(sprintf('%s=varargin{%f};',ARG(i,:),i)); end
  end


%   inp=repmat(struct(...
% 		   'y',varargin{1},...
% 		   'p',varargin{2},...
% 		   'PD',zeros(1,4),...
% 		   'RO',zeros(1,2),...
% 		   'PaP',zeros(1,floor(varargin{2}/2)+1),...
% 		   'TVPaP',zeros(2,floor(varargin{2}/2)),...
% 		   'fijo',0,...
% 		   'exog',0,...
% 		   'missing',0),...
% 	    1,1);

  if isempty(PaP)
    PaP=[inf, p ./ (1:floor(p/2))];
  end
  if isempty(RO)
    RO=[p+4 : min(floor(length(y)/3), ceil((p+12)*1.5))];
  end
  

  inp=repmat(struct(...
		    'y',y,...
		    'p',p,...
		    'RO',[p+4 : min(floor(length(y)/3), ceil((p+12)*1.5))],...
		    'PaP',[inf, p ./ (1:floor(p/2))],...
		    'PD',[.05,.45, (p+6)*2, min((length(y)-max(RO))*1.5, 600)],...
		    'TVPaP',[PaP>=2; PaP>2],...
		    'fijo',0,...
		    'exog',0,...
		    'missing',0),...
	     1,1);

  for i=1:nargin
    if ~isempty(varargin{i}); eval(sprintf('inp.%s=varargin{%f};',ARG(i,:),i)); end
  end

%   for i=1:size(ARG,1)
%     eval(sprintf('%s=inp.%s;',ARG(i,:),ARG(i,:)));
%   end

%   inp.RO    = [inp.p+4:min(floor(length(inp.y)/3), ceil((inp.p+12)*1.5))];
%   inp.PaP   = [inf, inp.p ./ (1:floor(inp.p/2))];
%   inp.PD    = [.05,.45, (inp.p+6)*2, min((length(inp.y)-max(inp.RO))*1.5, 600)];
%   inp.TVPaP = [inp.PaP>=2; inp.PaP>2];

  if any(inp.PaP~=Inf) && any(inp.PaP==Inf) 
    inp.PD(3) = max(inp.PD(3),3/2*max(inp.PaP(inp.PaP<Inf)));% deja espacio a un componente cíclico.
  end

%   for i=3:nargin
%     if ~isempty(varargin{i}); eval(sprintf('inp.%s=varargin{%f};',ARG(i,:),i)); end
%   end

%   RO, PD, inp.PD, pause

  for i=1:size(ARG,1); 
    eval(sprintf('%s=inp.%s;',ARG(i,:),ARG(i,:)));
    %eval(sprintf('%s=getfield(inp.%s);',ARG(i,:),ARG(i,:)));
  end

  %% Completamos la matriz TVP si incompleta
  if size(TVPaP,1)<2
    TVPaP = [TVPaP; zeros(size(TVPaP))];
  end
  if size(TVPaP,2)<size(PaP,2)
    TVPaP = [TVPaP, TVPaP(:,size(TVPaP,2))*ones(1,size(PaP,2)-size(TVPaP,2))];
  end
  %% ajustamos el orden minimo del polinomio AR ajustado a la serie y
  %% filtramos los datos por la parte AR (si el modelo está fijado)
  if fijo 
    A  = sum(sum(TVPaP(:,PaP==Inf|PaP==2)))+2*sum(sum(TVPaP(:,PaP<Inf&PaP>2)));
    B  = max(A,max(RO));
    P  = PaP; TVP = TVPaP;	% modelo fijado independientemente del orden AR
    [ARdhr,MAdhr,RARdhr,RMAdhr,ctes] = dhr2arma (TVP,P); 
    Ar = real(poly(RARdhr(RARdhr~=0)));  s1=size(Ar,2); s2=s1;
    x  = filter(Ar,1,y);	% filtrado de la serie
  else
    A = sum(sum(TVPaP(:,PaP==Inf|PaP==2)))+sum(sum(TVPaP(:,PaP<Inf&PaP>2)));
    B = max(A,max(RO));
  end
  RO = RO(RO>=A&RO<=B);		% ajustamos el orden minimo del polinomio AR 
  if isempty(RO); RO=A; end  	% ffprintf('\n AR(%u) fitted \n',A);

  OARok = zeros(B,1); OARok2 = zeros(B,1); mcnn = zeros(B,1); R2 = -ones(B,1)*10^1000;

  results = repmat(struct(...
			  'P',{zeros(size(PaP))},...
			  'TVP',{zeros(size(TVPaP))},...
			  'VAR',{zeros(1,size(TVPaP,2)+1)},...
			  'NVR',{zeros(size(PaP))},...
			  'proc',{'    '},...
			  'adjust',{-10^1000}),...
		   B,1);

% 			  'fixed',{fijo},...
% 			  'periodicity',{p},...
  
  RCAM = 0;
  %if PaP==inf&p==1,RCAM=1;else,RCAM=0;end % RC a la Agustin si HP

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%% AQUI EMPIEZA EL ALGORITMO de ESTIMACION E IDENTIFICACION 
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %% ESTIMACION DE LOS NVR Y LAS VARIANZAS PARA DISTINTOS ORDENES DE AR
  for o = RO			% Ordenes AR que evaluamos
    
    %% DETERMINACION DEL MODELO DHR Y COMPONENTES AUXILIARES
    if fijo %% Parte Ar del polinomio AR=conv(Ar,aR)      
      aR      = ar_mdcov(x(s1:length(y)),o-s2+1); 
    else %% IDENTIFICACION DEL MODELO DHR      
      if exog	   % estimacion del modelo AR con variables explicativas
	Ar    = modar(y,o,exog);
      else			% estimacion del modelo AR
	Ar    = ar_mdcov(y,o);
      end      
      TVPp    = 2*ones(size(TVPaP)); % distinto de TVP para iniciar el 'while'
      aR      = [];   
      iter    = 0;   
      PD2     = PD;
      PD2(1)  = min(.01,PD(1)); 
      PD2(3)  = max((p+6)*3,PD(3)); 
      PD2(4)  = min((length(y)-max(RO))*2,max(800,PD(4)));      
      [TVP,P] = iddhr (p,Ar,PD,PaP,TVPaP,aR); % identificacion inicial
      TVPp    = zeros(size(TVP));	% para iniciar el 'while'
      %% busqueda iterativa
      while  iter<1 || ( (~all(all(TVP==TVPaP)) &&  any(any(round(TVP*100)~=round(TVPp*100)))) && iter<5 ) 
	TVPp  = TVP;		% guardo modelo actual para comparar con el siguiente
	TVP(TVP~=TVPaP) = 0;	% sólo mantengo los parametros de TVP idénticos a TVPaP
	[ARdhr,MAdhr,RARdhr,RMAdhr,ctes] = dhr2arma (TVP,P); 
	Ar = real(poly(RARdhr(RARdhr~=0)));	
	x  = (filter(Ar,1,y)); s1 = size(Ar,2); s = s1;
	aR = ar_mdcov(x(s:length(y)),o-s+1); % resto del AR (parte aR)
	
	%% identificacion modelo DHR y comp aux
  	if iter==1		% inicialmente parámentros estrictos
 	  PD2(1) = min(.01,PD(1)); PD2(3) = max((p+6)*3,PD(3)); 
 	  PD2(4) = min((length(y)-max(RO))*2,max(800,PD(4)));
	elseif iter==2		% 2ª vuelta se relajan
	  PD2(1) = min(.02,PD(1)); PD2(3) = max((p+6)*2.5,PD(3));
	  PD2(4) = min((length(y)-max(RO))*1.75,max(700,PD(4)));
	else
	  PD2 = PD;
	end			% en las sucesivas uso de PD por defecto
	
	[TVP,P]  = iddhr (p,Ar,PD,PaP,TVPaP,aR);
	[ARdhr,MAdhr,RARdhr,RMAdhr,ctes] = dhr2arma (TVP,P);
	iter     = iter+1;
      end

      [ARdhr,MAdhr,RARdhr,RMAdhr,ctes] = dhr2arma (TVP,P);
      Ar = real(poly(RARdhr(RARdhr~=0)));  s1=size(Ar,2); s2=s1;
      x  = filter(Ar,1,y);	% filtrado de la serie
      aR      = ar_mdcov(x(s1:length(y)),o-s2+1); 
    end

    [Pac,RARac] = arAC2roots (aR);
    
    PP   = [P, Pac];
    RAR  = [RARdhr, RARac];
    RMA  = [RMAdhr, zeros(size(RMAdhr,1), size(RARac,2))];
    CTES = [ctes, ones(1,size(RARac,2))];
    
    lim      = .5; %% clasificación de raices    
    RARE     = RAR;        RARE(abs(1-abs(RAR))<lim)       = 0;
    RARNE    = RAR;        RARNE(abs(1-abs(RAR))>=lim)     = 0;    
    RARacNE  = RARac;      RARacNE(abs(1-abs(RARac))>=lim) = 0;
    racne    = RARacNE(:);                           racne = racne(racne~=0);
    polyacNE = real(poly(racne));
    
    re    = RARE(:);  re  = re(re~=0);
    polyE = real(poly(re));
    
    r     = RAR(:);   r   = r(r~=0);
    polyR = real(poly(r));
    
    rne   = RARNE(:); rne = rne(rne~=0);
    polyNE= real(poly(rne));
    
    xx=(filter(polyR,1,y)); xx=xx([size(polyR,2):length(y)]); % serie filtrada R
    
    rxx   = size(xx,1)-1;  
    fxx   = (abs(fft(xx - mean (xx)))).^2/rxx; fxx(1) = []; % periodograma serie filtrada
    nfrec = reshape(2*pi*linspace(1/rxx,(rxx-1)/rxx,rxx),rxx,1);
    
    fy    = fxx.*real(esptarma(1,polyE,nfrec)); % espectro
    
    S     = zeros(rxx,size(PP,2));
    z     = exp(-i*nfrec);
    preci = 10^16;
    
    for c = 1:size(PP,2)
      r_resto_ne = reshape(RARNE(:,1:size(PP,2)~=c),size(RAR,1)*(size(PP,2)-1),1);
      ma         = RMA(:,c);
      re         = RARE(:,c);
      num        = [ma(ma~=0);
	            r_resto_ne(r_resto_ne~=0)];
      Mnum(1:size(real(poly(num)),2),c) = CTES(c)*fix(preci.*real(poly(num))')./preci;
      den        = re(re~=0);
      Mden(1:size(real(poly(den)),2),c) = fix(preci.*real(poly(den))')./preci;
      E          = polyval(real(poly(num)),z)./polyval(real(poly(den)),z);
      %%      S(:,c)=[E.*conj(E)*CTES(c)];
      S(:,c)     = esptarma (real(poly(num)),real(poly(den)),nfrec,CTES(c).^2);
    end
    
    CirrT = real(esptarma(polyNE,1,nfrec)); % Componente irreg transf
    
    N = find(isnan(fy)); fy(N)=[]; fxx(N)=[]; nfrec(N)=[]; S(N,:)=[]; CirrT(N)=[];
    
    Sdhr  = S(:,1:size(P,2));	  % espectro componentes DHR transf
    Sac   = S(:,size(P,2)+1:end); % espectro componentes ac transf
    ncdhr = size(Sdhr,2);
    
    %% %  if p==1, marg=pi/p;                 % todas las freq
    marg = min(1,.8+p*2/length(y)) * pi/p;
    A    = 2*pi./P; 
    F    = [A,2*pi-fliplr(A(A<pi))];
    in   = F'-marg; 
    IN   = in(:,ones(1,rxx));
    su   = F'+marg;
    SU   = su(:,ones(1,rxx));
    w    = any(IN<=nfrec(:,ones(length(F),1))'&nfrec(:,ones(length(F),1))'<=SU);
    
    %% esto quitaria las frecuencias que estan 'entorno' a los 'polos'
    %% FFF=[2*pi./PaP , 2*pi-2*pi ./fliplr(PaP(PaP~=2))];
    %% w(find((abs(nfrec(:,ones(size(FFF,2),1))-FFF(ones(size(nfrec,1),1),:))<10^-6)))=[];
    
    %% ortogonalizando parte DHR de componentes aux
    coef     = pinv([Sac(w,:),CirrT(w)])*Sdhr(w,:);
    GSO(w,:) = Sdhr(w,:)-[Sac(w,:),CirrT(w)]*coef; 
    VAR      = pinv([GSO(w,:)])*fxx(w); % estimacion varianza componentes
    
    CirrE = fy(w,:)-[Sdhr(w,:)*VAR(1:ncdhr)]; % Espectro Comp. irr empirico
    
    TransF = esptarma(1,polyacNE,nfrec); % CirrT con raices ac / sin raices DHR
    CirrE  = fy-[Sdhr*VAR(1:ncdhr)];	 % Espectro Comp. irr empirico
    irrE   = CirrE.*TransF;
    irrT   = CirrT.*TransF; 
    
    s2_irr = max( pinv(irrT)*irrE, pinv(CirrT)*CirrE );
    
    if  any([VAR;s2_irr]<0) && exist ('lsqnonneg', 'file') %|s2_irr<max(VAR)*2
      NNN     = lsqnonneg([Sdhr,CirrT],fy);
      VAR     = NNN(1:ncdhr); 
      s2_irr  = NNN(ncdhr+1);
      mcnn(o) = 1;
    elseif  any([VAR;s2_irr]<0) && exist ('nnls', 'file') % & ~exist ('lsqnonneg', 'file') % | s2_irr<max(VAR)*2
      tol     = 10*eps*norm([Sdhr,CirrT],1)*max(size([Sdhr,CirrT]));
      NNN     = nnls([Sdhr,CirrT],fy,tol);
      VAR     = NNN(1:ncdhr); 
      s2_irr  = NNN(ncdhr+1);
      mcnn(o) = 1;
    end
    
    Var = VAR(1:ncdhr);
    
    [theta,din,lab,AR,MA,s2a] = dhr2rf (TVP,P,[s2_irr,VAR'],[],p);

    xx     = (filter(AR,1,y)); 
    xx     = xx(size(AR,2):length(y));
    rx     = size(xx,1);  
    nfrec  = [1:rx-1]/rx*2*pi;
    Pxx    = (abs(fft(xx - mean (xx)))).^2/rx; % Periodograma
    Pxx(1) = [];
    
    teo    = esptarma(MA,1,nfrec',s2a);
    r      = (log(Pxx)-log(teo));
    R2(o)  = 1-dot(r(w),r(w))/dot(log(teo(w)),log(teo(w)));

    %% ALMACENANDO LOS RESULTADOS EN LA ESTRUCTURA 'model'
    if mcnn(o)
      TIPO = 'NNLS';
    else
      TIPO = 'OLS ';
    end
    
    results(o).P      = P; 
    results(o).TVP    = TVP;
    results(o).VAR    = [s2_irr,Var'];
    if s2_irr; 
      results(o).NVR  = [Var'./s2_irr];
    end
    results(o).proc   = TIPO;
    results(o).adjust = R2(o);
    
    %% SELECCIONAMOS LOS CASOS DONDE LAS VARIANZAS SON POSITIVAS; y
    %% además NVRs < 1 con estimacion MCO
    if all([[s2_irr,Var']>=0,s2_irr>Var']) && ~mcnn(o)
      OARok(o)  = o; 
      OARok2(o) = o;
    elseif all(s2_irr>0&Var>=0) % o al menos varianzas positivas 
      OARok2(o) = o;
    end
  end
  
  if any(OARok) 		% si todo bien
    OAR = OARok;
  elseif any(OARok2)		% si NVR's > 1, o NNLS
    OAR = OARok2;
    % disp('\n NVRs mayores que 1 o estimacion por MCNN')
    % disp('Verifique que la periodicidad es correcta;')
    % disp('y que los parametros de decision y el modelo son adecuados.')
    % disp('Si es asi, amplie la horquilla de ordenes AR\n')
  end
  
  ff=[];
  if exist('OAR') 		% no hay grandes problemas
    f = find(OAR)'; v = [1:2]';	% 2 es num de raices componente
    %% si modelos identificados iguales (empleamos modelos en la mediana)
    if all(all(cat(1,results(f).TVP)==cat(1,results(f(1)).TVP(vec(v(:,ones(length(f),1)))',:)))) % iguales al primero
      VAR       = cat(1,results(f).VAR);  
      NVR       = VAR(:,2:length(PaP)+1)./(VAR(:,1)*ones(1,length(PaP)));
      m         = median(NVR,1);
      M         = m(ones(size(m,2),1),:); 
      me        = max(max(abs(M-M'))); 
      mi        = min(max(abs(M-M')));
      DD        = abs(NVR-m(ones(size(NVR,1),1),:)); 
      DD(DD==0) = inf;
      
      medianos  = find(all(NVR>=m(ones(size(NVR,1),1),:)-max(min(DD))-10*eps&NVR<=m(ones(size(NVR,1),1),:)+max(min(DD))+10*eps,2));

      if isempty(medianos)
	medianos   = find(all(NVR>=m(ones(size(NVR,1),1),:)-mi-10*eps&NVR<=m(ones(size(NVR,1),1),:)+mi+10*eps,2));
	if isempty(medianos)
	  medianos = find(all(NVR>=m(ones(size(NVR,1),1),:)-me-10*eps&NVR<=m(ones(size(NVR,1),1),:)+me+10*eps,2));
	end
      end
      ff  = f(medianos);
      if ~isempty(ff)
	f = ff; 
      end
    end
    
    %% Ordenamos por ajuste espectral
    [a,reord] = sort(-R2(f));  
    mejores   = f(reord); 
    oar       = mejores(1);

    f         = find(OARok2); 
    [a,reord] = sort(-R2(f));  
    mejores2  = f(reord);     
  else 				% si sólo estimaciones negativas
    
    [a,mejores2] = sort(-R2);	%mejores2=mejores2(1:min(6,(B-A)+1));
    oar          = mejores2(1);
    
    disp('\n ERROR... Estimacion NEGATIVA por MCO!');
    disp('Verifique que la periodicidad es correcta;');
    disp('y que los parametros de decision y el modelo son adecuados.');
    disp('Si es asi, amplie la horquilla de ordenes AR\n');
  end
  
  %if mcnn(oar); disp('\n NON-NEGATIVE LEAST SQUARES ESTIMATION\n'); end
  
  VAR  = results(oar).VAR; 	% Varianzas
  P    = results(oar).P;	% Periodos de los Componentes
  TVP  = results(oar).TVP; 	% Modelo DHR
  NVR  = results(oar).NVR; 	% Modelo DHR
  MCNN = mcnn(oar);
  
  %% meter modelo ARIMA implicito en estructura				  
  %%[theta,din,lab,AR,MA,s2a] = dhr2thd (TVP,P,VAR,[],p);
  
  if nargout > 6     
    %% dibujitos
    [theta,din,lab,AR,MA,s2a] =  dhrgraph (y,TVP,P,VAR,p,oar);
    
    %% mostramos los resultados de los mejores modelos    
    dhrshowmodels (PaP,mejores2,results,mcnn,R2,NVR,oar);
  end

 

%%%% example ('autodhr')
%%%% demo ('autodhr', 2)


%!demo
%! clf
%! x = 0:0.1:2*pi; 
%! y1 = sin (x);
%! y2 = exp (x - 1);
%! ax = plotyy (x, y1, x - 1, y2, @plot, @semilogy);
%! xlabel ('X');
%! ylabel (ax(1), 'Axis 1');
%! ylabel (ax(2), 'Axis 2');

%!demo
%! clf
%! x = linspace (-1, 1, 201);
%! subplot (2, 2, 1)
%! plotyy (x, sin(pi*x), x, 10*cos(pi*x))
%! subplot (2, 2, 2)
%! surf (peaks (25))
%! subplot (2, 2, 3)
%! contour (peaks (25))
%! subplot (2, 2, 4)
%! plotyy (x, 10*sin(2*pi*x), x, cos(2*pi*x))
%! axis square
