function [F, dF, A, dA, V, dV, G, dG] = estr_dv(theta,din,di)
% ESTR_DV  - Computes the partial derivatives of the matrices in a (ESTR) 
% structural model with respect to the i-th parameter of theta vector.
%   [F, dF, A, dA, V, dV, G, dG] = estr_dv(theta, din, i)
%
% 5/3/97
% Copyright (c) Jaime Terceiro, 1997

global EEOPTION

mtipo = din(1,1); np = din(1,6); di = fix(di);
if mtipo ~= 2, e4error(14); end
if np ~= size(theta,1), e4error(1); end
if any(size(di)~=1), e4error(2); end
m = din(1,2); r = din(1,3); s = din(1,4); k = din(1,5); diagv = din(1,7);
din2 = din(2:m+1,:); % Para simplificar
ardin  = din2(1:m,1:2*m);  madin  = din2(1:m,2*m+1:2*m+m);
inpdin = []; if r > 0, inpdin = din2(1:m,3*m+1:3*m+2*r); end
p  = max(max(ardin)); q  = max(max(madin));
if r, g = max(max(inpdin)); else g = -1; end
k1 = max([p, q, g]);
if (k1 ~= k), e4error(1); end

count = 1;
cc = zeros(m,m);
dtheta = zeros(np,1); dtheta(di) = 1;
F = eye(m, (1+k)*m); dF = zeros(m, (1+k)*m);
A = eye(m, (1+k)*m); dA = zeros(m, (1+k)*m);
G = zeros(m, (1+k)*r); dG = zeros(m, (1+k)*r);

Ft = cc;
dintmp = (ardin(:, 1:2:2*m) == 0);
dintmp(eye(m)) = zeros(m,1);
np = sum(sum(dintmp));
Ft(dintmp') = theta(count:count+np-1);
Ft(eye(m)) = ones(m,1);
F(:,1:m) = Ft';

Ft = cc;
Ft(dintmp') = dtheta(count:count+np-1);
dF(:,1:m) = Ft';
count = count+np;

if size(din,2) > 7
   if din(1,8) % VARMAX ECHELON STRUCTURE
      A(:,1:m)  = F(:,1:m); 
      dA(:,1:m) = Ft';
   end
end

for l = 1:p
    Ft = cc;
    dintmp = ( (ardin(:, 1:2:2*m) <= l) & (ardin(:, 2:2:2*m) >=l)  );
    np = sum(sum(dintmp));
    Ft(dintmp') = theta(count:count+np-1);
    F(:,l*m+1:(l+1)*m) = Ft';
    Ft(dintmp') = dtheta(count:count+np-1);
    dF(:,l*m+1:(l+1)*m) = Ft';
    count = count+np;
end

for l = 1:q
    Ft = cc;
    dintmp = (madin > 0);
    np = sum(sum(dintmp));
    Ft(dintmp') = theta(count:count+np-1);
    A(:,l*m+1:(l+1)*m) = Ft';
    Ft(dintmp') = dtheta(count:count+np-1);
    dA(:,l*m+1:(l+1)*m) = Ft';
    madin = madin - 1;
    count = count+np;
end

cc = zeros(r,m);
for l = 0:g
    Ft = cc;
    dintmp = ( (inpdin(:, 1:2:2*r) <= l) & (inpdin(:, 2:2:2*r) >=l)  );
    np = sum(sum(dintmp));
    Ft(dintmp') = theta(count:count+np-1);
    G(:,l*r+1:(l+1)*r) = Ft';
    Ft(dintmp') = dtheta(count:count+np-1);
    dG(:,l*r+1:(l+1)*r) = Ft';
    count = count+np;
end

V = zeros(m,m); dV = V;
if diagv
   V = diag(theta(count:count+m-1));
   dV = diag(dtheta(count:count+m-1));
else
   cc = ones(m,m);
   V(triu(ones(m,m))) = theta(count:count+(m*(m+1)/2)-1);
   dV(triu(ones(m,m))) = dtheta(count:count+(m*(m+1)/2)-1);
   if ~EEOPTION(5)
        V( tril(cc,-1) ) = V( tril(cc,-1)' );
        dV( tril(cc,-1) ) = dV( tril(cc,-1)' );
   else
        V = V';
        dV = dV';
   end
end

% Force V positive definite
if ~EEOPTION(5)
   if min(eig(V)) <= 0
      Vu = cholp(V);
      V  = Vu'*Vu;
   end
else
    dV = dV * V' + V * dV';
    V = V*V';
end