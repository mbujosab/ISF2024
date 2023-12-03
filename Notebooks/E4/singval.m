function [S] = singval(y, i, k, s)

% Funcion interna para el calculo de correlaciones canonicas y no canonicas
% S = singval2(y, i, k, s)
% k = 1 -- No canonica
% k = 0 -- Canonica

if nargin < 4, s = 1; end
ext = 1; % Siempre con matrices BHK extendidas
j = i;   % dim(pasado) = dim(futuro)

[N m] = size(y);    
Yi = blkhkel(y, i+j, s, ext);
N = min(size(Yi,2), size(y,1));
     
[Q, R] = qr(Yi', 0);
R = R';
ix = zeros(3,2);
ix(:,1) = [1; i*m+1; (i+1)*m+1];
ix(:,2) = [ix(2:3,1)-1; (i+j)*m];

if k
  % No canonica
   [U S V] = svd(R(ix(2,1):ix(3,2),ix(2,1):ix(3,2)) \ R(ix(2,1):ix(3,2),ix(1,1):ix(1,2)));
else   
  % Canonica
   Om = R(ix(2,1):ix(3,2),ix(1,1):ix(3,2))*R(ix(2,1):ix(3,2),ix(1,1):ix(3,2))';
   isOm = pinv(sqrtm(Om));
   [U S V] = svd(isOm*R(ix(2,1):ix(3,2),ix(1,1):ix(1,2)));       
end     
S = diag(S);