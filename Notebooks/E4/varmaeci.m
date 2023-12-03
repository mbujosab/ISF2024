function [Fx,sFx,Thx,sThx,Sigmax,A1x,A2x,A3x] = varmaeci(y, u, i, ki, s)
%
% Ver varmaech.m para sinopsis. La diferencia es que puede recibir en i un vector con
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
   [F,sF,Th,sTh,Sigma,A1,A2,A3] = varmaech(y, u, i(k), ki, s);
   f = det(Sigma);

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
      Fx = F;
      sFx = sF;
      Thx   = Th;
      sThx   = sTh;
      Sigmax = Sigma;
      A1x = A1;
      A2x = A2;
      A3x = A3;
   end
end
