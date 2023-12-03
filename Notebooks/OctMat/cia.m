%% usage: [AR,s2_e,n_retardos] = cia (serie,rango_retardos,IA) 
%% 
%% Calcula el orden de un polinomio autorregresivo según el  
%% Criterio de Información de Akaike 
%%  
%% Serie          > serie temporal  
%% rango_retardos > [orden mínimo retardos, orden máximo] % ([1 T/3]) 
%%                    (Si escalar => ese sera el orden del modelo AR) 
%%                    (T = longitud de la muestra) 
%% IA             > Información Adicional                 % (0) 
%%                    (>0 calcula raíces, =2 información por pantalla)  
%% 
%% AR             > polinomio autorregresivo estimado 
%% s2_e           > varianza residual del modelo AR 
%% n_retardos     > numero de retardos empleados en el modelo AR 
%% 
%% Copyright (C) 1999, 2000, 2001 Marcos Bujosa. 
 
%% Author: Marcos Bujosa <marcos.bujosa@ccee.ucm.es> 
%% References: Box, G. E. P., Jenkins, G. M., and Reinsel, G. C. (1994).
%% Time series analysis: forecasting and control. Prentice Hall, Inc.,
 
function [AR,s2_e,n_retardos] = cia (y,rango_retardos,IA) 

  T=length(y); y=y(:); 
   
  ARG=['y             ';
       'rango_retardos';
       'IA            ';];
  
  for i=nargin+1:size(ARG,1)
    eval(sprintf('%s=[];',ARG(i,:)))
  end

  %% VALORES POR DEFECTO 
  if isempty(rango_retardos); rango_retardos=1:(T/3); end 
  if isempty(IA), IA=0; end

  rango_retardos=sort(rango_retardos);
  %% CÁLCULO DEL MODELO AR 
  if size(rango_retardos)==[1 1] % si se fija el número de retardos
    n_retardos=rango_retardos;
    [AR,s2_e]=modar(y,rango_retardos); % Modelo Autorregresivo        
  else				% si se da un rango de retardos 
    j=0; 
    for p=rango_retardos
      j=j+1; 
      [AR,s2_e]=modar(y,p);	% estima el Modelo Autorregresivo 
      AIC(j,1)=log(s2_e)+p*2/(T-p); % Crit. de Inf. de Akaike  
    end 
    [aic,posicion]=min(AIC);
    
    n_retardos=rango_retardos(posicion); % nº retardos óptimo
    [AR,s2_e]=modar(y,n_retardos); % estima el Modelo Autorregres. 
  end     

  %% INFORMACIÓN ADICIONAL 
  if IA>1 
    fprintf('\n\tAR(%1.f)\n\n',n_retardos); 
    if length(rango_retardos)>1 % si se da un rango de retardos 
      plot(rango_retardos,AIC);
      title(sprintf('Criterio de Información de Akaike. Mínimo en %g retardos',... 
       n_retardos))
      xlabel('Número de retardos')
      ylabel('Valor de la función AIC')
    end 
    if n_retardos>floor(T/3) 
      disp('Número de retardos mayor que el 33% del tamaño de la muestra') 
    end 
  end   
