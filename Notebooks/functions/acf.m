% Calcula (asintóticamente[ N, no N-k]) la función de autocovarianza o 
% autocovarianza cruzada. Con un sólo argumento calcula la autocovarianza, 
% con dos argumentos calcula la autocovarianza cruzada, el tercero indica el  
% número de retardos a calcular (por defecto N/3). 
% FACOVA(X) calcula con retardos k=0...N/3,(vector columna N/3 +1 X 1) 
% FACOVA(X,Y) calcula con retardos k=-N/3...N/3,(vector columna [2*N/3 +1 ,  1]  ) 
% FACOVA(X,Y,R) calcula con retardos k=-R...R 
% gamma=E[(x(t)-E(x))(y(t+k)-E(y(t+k))]/N-k 
 
function retval = acf(X,Y,h) 

 
  if nargin==1		% necesita la segunda serie y el número de retardos 
    Y=X; 
    h=length(X)/3; 
  elseif nargin==2	% necesita el número de retardos 
    if max(size(X))~=max(size(Y)); 
      error("input vectors must have the same length")
    endif
    if min(size(X))~=1 | min(size(Y))~=1; 
      error("bad arguments")
    endif
    h=length(X)/3; 
  endif
  
  X=X(:);Y=Y(:); 
  
  [n, c] = size (X);

  X = X-mean(X);
  Y = Y-mean(Y);

  retval = zeros (h + 1, 1);

  for i = 0 : h
    retval(i+1) = diag (X(i+1:n).' * conj (Y(1:n-i))).' / n;
  endfor

endfunction
