function [Phix,sPhix,Hx,sHx,Ex,sEx,Qx,Dx,sDx,innovx] = sidinti(y, u, i, ki, s)
%
% Ver sidint.m para sinopsis. La diferencia es que puede recibir en i un vector con
% posibles valores. Se seleccionará el valor de i que minimice la varianza del error de
% predicción

i = i(:);
r = size(u,2);
m = size(y,2);
n = sum(ki);

if nargin < 5, s = 1; end

fmax = inf;
for k=1:size(i,1)
   [Phi,sPhi,H,sH,E,sE,Q,D,sD,innov] = sidint(y, u, i(k), ki, s);
   [thet, di] = ee2thd(Phi, zeros(n,r), E+Phi(:,1:m), [H zeros(m,n-m)], D, [], Q);
   f = pem(thet,di,[y u]);
   if f < fmax
      fmax = f;
      Phix = Phi;
      sPhix = sPhi;
      Hx  = H;
      sHx = sH;
      Ex  = E;
      sEx = sE;
      Qx  = Q;
      Dx  = D;
      sDx = sD;
      innovx = innov;
   end
end
