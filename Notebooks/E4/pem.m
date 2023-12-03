function f = pem(theta, din, z)
% FVARMAX - Función de verosimilitud exacta de un modelo ARMAX, TF, ESTR o cualquier
%    otro en el que se verifique que Q = R = S y C = I, suponiendo que el vector de
%    estado inicial es un parámetro a estimar, sin incertidumbre. Cuando se verifican
%    las condiciones expuestas, esta función proporciona exactamente el mismo
%    resultado que fvmod() con vcond = IDJONG, y es VALIDA tanto para sistemas
%    ESTACIONARIOS como NO ESTACIONARIOS.
%	     [f{, innov, ssvect}] = fvarmax(theta, din, z, fown)
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

if nargin < 3,  e4error(3); end

[Phi, Gam, E, H, D, C, Q, S, R] = thd2ss(theta, din);
n=size(z,1);
l = size(Phi,1);
m = size(H,1);
r=max([size(Gam,2), size(D,2)]);

x0 = zeros(l,1);

z1 = zeros(m,n);

for t = 1:n  % main loop
    if r
    %
       z1(:,t)  = z(t,1:m)' - H*x0 - D*z(t,m+1:m+r)';
       x0  = Phi*x0 + Gam*z(t,m+1:m+r)' + E*z1(:,t);
    %
    else
    %
       z1(:,t)  = z(t,:)' - H*x0;
       x0  = Phi*x0 + E*z1(:,t);
    %
    end
end  % t

f = log(det(z1*z1'/n));
if f == -inf, f = inf; end

