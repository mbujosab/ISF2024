function [f, var, beta, sbeta, z1] = fvarima(theta, din, z, opt1, opt2)
%
%	     [f, var, beta, sbeta, z1] = fvarima(theta, din, z, fown)
%   theta    > vector de parámetros
%   din      > matriz que describe la dinámica del modelo
%   z        > matriz de variables observables.
%   fown     > nombre de la función de usuario (sólo en modelos de usuario)
%   f       <  valor de la función de verosimilitud.
%   innov   <  (opcional) almacena la secuencia de innovaciones.
%   ssvect  <  (opcional) secuencia del vector de estado.
% 7/3/96
% (C@@)

global EEOPTION
saveinn = 0; if nargout >= 5, saveinn = 1; end

econd = EEOPTION(4);
zeps = EEOPTION(15);

if nargin < 4
   opt1 = [];
end

n=size(z,1);
r=size(z,2);

[Phi, Gam, E, H, D, C, Q, S, R] = thd2ee(theta, din, opt1);
m = size(H,1);

[Sigm, nonstat] = lyapunov(Phi, E*E');  % control de no estacionariedad
if nonstat > 0, e4error(26); end

l = size(Phi,1);
x0 = zeros(l,1); X0U = zeros(l,r-1);
WW  = zeros(l,l); WZ = zeros(l,1); WU = zeros(l,r-1);
Phib = (Phi - E*H);
Phibb0 = H;

if saveinn, N = eye(l); end

z1 = zeros(n,1);
u1 = zeros(n,r-1);

for t = 1:n
%
    z1(t)  = z(t,1)' - H*x0;
    x0  = Phi*x0 + E*z1(t);

    u1(t,:)  = z(t,2:r) - H*X0U;
    X0U  = Phi*X0U + E*u1(t,:);
    
    WW  = WW + Phibb0'*Phibb0;
    WZ  = WZ + Phibb0'*z1(t);
    WU  = WU + Phibb0'*u1(t,:);
    Phibb0 = Phibb0*Phib;
end

if any(isnan([WW WZ])) | any(isinf([WW WZ]))
   B1, WW, WZ
   e4error(25);
else
   [Ns S Ns] = svd(Sigm);
   k = find(diag(S) < zeps);
   if size(k,1) < l
      if size(k,1)
      %   
         N = Ns(:,1:k(1)-1);
         Sigm = S(1:k(1)-1,1:k(1)-1);
         WW = N'*WW*N;
         WZ = N'*WZ;
         WU = N'*WU;
      %
      end
      M = cholp(Sigm);
      T = cholp(eye(size(M,1))+M*WW*M');

      if econd == 2
         T2 = cholp(WW);
         TWU = T2'\WU;
         TWZ = T2'\WZ;
      else
         TWU = T'\(M*WU);
         TWZ = T'\(M*WZ);
      end

      UU = u1'*u1 - TWU'*TWU;
      UZ = u1'*z1 - TWU'*TWZ;
      iUU   = inv(UU);     
      beta  = iUU*UZ;
      sbeta = sqrt(diag(iUU));

      var  = (z1'*z1 - TWZ'*TWZ - beta'*UZ)/n;

      ff = n*sum(log(var))+n + 2*sum(log(diag(T)));

   end
   %
end

f = 0.5*(ff + n*m*log(2*pi));

if saveinn
%
   x0 = N*pinv(WW)*(WZ - WU*beta);

   for t = 1:n
       z1(t)  = z(t,1) - H*x0 - z(t,2:r)*beta;
       x0  = Phi*x0 + E*z1(t);
   end
%
end
