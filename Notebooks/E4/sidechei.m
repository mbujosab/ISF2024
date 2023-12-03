function [Phix,sPhix,Hx,sHx,Ex,sEx,Qx,A1x,A2x,A3x,A4x,A5x] = sidechei(y, u, i, ki, s)
%
% Ver sidechel.m para sinopsis. La diferencia es que puede recibir en i un vector con
% posibles valores. Se seleccionará el valor de i que minimice la varianza del error de
% predicción

i = i(:);
r = size(u,2);
m = size(y,2);
n = sum(ki);
mb = sum(ki > 0);
ik = find(ki > 0);

if nargin < 5, s = 1; end

fmax = inf;
for k=1:size(i,1)
   [Phi,sPhi,H,sH,E,sE,Q,A1,A2,A3,A4,A5] = sidechel(y, u, i(k), ki, s);
   f = det(Q);

%   if r
%      Phi2 = Phi; E2 = E;
%      Phi2(:,1:mb) = Phi(:,1:mb)*H(ik,:);
%      E2(:,ik) = E(:,ik)+Phi2(:,1:mb);
%      [thet, di] = ee2thd(Phi, A1, E2, [H zeros(m,n-m)], A3, [], Q);
%      f = pem(thet,di,[y u]);
%   else
%      Phi2 = Phi; E2 = E;
%      Phi2(:,1:mb) = Phi(:,1:mb)*H(ik,:);
%      E2(:,ik) = E(:,ik)+Phi2(:,1:mb);
%      [thet, di] = ee2thd(Phi2, [], E2, [H zeros(m,n-mb)], [], [], Q);
%      f = pem(thet,di,y);
%   end
   if f < fmax
      fmax = f;
      Phix = Phi;
      sPhix = sPhi;
      Hx   = H;
      sHx   = sH;
      Ex   = E;
      sEx   = sE;
      Qx   = Q;
      A1x = A1;
      A2x = A2;
      A3x = A3;
      A4x = A4;
      A5x = A5;
   end
end
