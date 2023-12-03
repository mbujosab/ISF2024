function [Phi, Gam, E, H, D, C, Q, S, R] = diff(th,dh)
%
% [Phi, Gam, E, H, D, C, Q, S, R] = diff(theta,din)
% 
% Devuelve la estructura del modelo resultado de aplicar una diferencia regular al modelo de entrada
%
%
zeps = .00001*1e-3;
[H_D, type, m, r, s, n, np,userflag, userf, innov] = e4gthead(dh);
dh = tomod(dh);
[Phi, Gam, E, H, D, C, Q] = thd2ss(th, dh);

if ~n
    Phi = [];
    E = [];
    Gam = [];
    H = [];
end
Phi = [Phi zeros(n,m);H zeros(m)];
if innov(1)
   E = [E;C];
   S = Q;
   R = Q;
else
   E = [E zeros(n,innov(3));zeros(m,innov(2)) C];
   Q = [Q zeros(innov(2),innov(3));zeros(innov(3),innov(2)) R];
   S = [S;eye(innov(3))];
end
 
if r, Gam = [Gam;D]; end

H = [H -eye(m)];

if r, B = [Gam, E]; else, B = E; end
   
[Phi,B,H] = minreal(Phi,B,H,zeros(m,size(B,2)),zeps);
if r
   Gam = B(:,1:r);
   E = B(:,r+1:size(B,2));
else
   E = B;
end

  