function [indice, valor]=maximos(y,n) 
   
  %% usage:  [indice, valor] = maximos (y,n) 
  %% 
  %% localiza los máximos locales de un vector de datos 
  %% indice es la posición del máximo local,  
  %% valor es el valor que toma el vector en esa posición. 
  %% 
  %% Si se introduce el argumento "n", sólo localiza los n 
  %% máximos de mayor valor.
  %% 
  %% Marcos Bujosa 1998, 1999, 2000, 2002
 
  indice=[]; valor=[]; 
  T=length(y); 
   
  if T==1  
    indice=1; valor=y;  
  else  
 
    if y(1)>=y(2) 
      indice=1; valor=y(1); 
    end 
     
    for i=2:T-1 
      if y(i)-y(i-1)>=0 & y(i)-y(i+1)>=0 
	if isempty(indice)==1; 
	  indice=i; valor=y(i); 
	else 
	  indice=[indice;i]; valor=[valor; y(i)]; 
	end 
      end 
    end 
     
    if y(T)>=y(T-1) 
      if isempty(indice)==1; 
	indice=T; valor=y(T); 
      else 
	indice=[indice;T]; valor=[valor; y(T)]; 
      end 
    end   
  end  
  
  if nargin>1
    
    if length(valor)<n; 
      n=length(valor);
    end
    [s,i]  = sort(-valor);
    i=i(1:n);
    indice = indice(i);
    valor  = valor(i);
  end
  
 