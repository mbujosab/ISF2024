function [F, A, V, G] = thd2estr(theta, din)
% THD2ESTR - Converts a THD model to the equivalent ESTR formulation.
%    [F, A, V, G] = thd2estr(theta, din)
% For the model:
%    Fb(B)y(t) = Gb(B)u(t) + Ab(B)e(t)
% returns the matrices:
% F = [F0 F1 ... Fk], A = [I A1 ... Ak], G = [G0 G1 ... Gk], V = V(e(t))
% k = max(p, q, g)
% The procedure also checks that V is positive definite.
%
% 7/3/97
% Copyright (c) Jaime Terceiro, 1997

global EEOPTION

mtipo = din(1,1); np = din(1,6);
if np ~= size(theta,1), e4error(1); end
m = din(1,2); r = din(1,3); s = din(1,4); k = din(1,5); diagv = din(1,7);
din2 = din(2:m+1,:);
ardin  = din2(1:m,1:2*m);  madin  = din2(1:m,2*m+1:2*m+m);
inpdin = []; if r > 0, inpdin = din2(1:m,3*m+1:3*m+2*r); end
p  = max(max(ardin)); q  = max(max(madin));
if r, g = max(max(inpdin)); else g = -1; end
k1 = max([p, q, g]);
if (k1 ~= k), e4error(5); end

count = 1;
cc = zeros(m,m);
F  = eye(m, (1+k)*m);
A  = eye(m, (1+k)*m);
G  = zeros(m, (1+k)*r);

Ft = cc;
dintmp = (ardin(:, 1:2:2*m) == 0);
dintmp(eye(m)) = zeros(m,1);
np = sum(sum(dintmp));
Ft(dintmp') = theta(count:count+np-1);
Ft(eye(m)) = ones(m,1);
F(:,1:m) = Ft';

if size(din,2) > 7
   if din(1,8) % VARMAX ECHELON STRUCTURE
      A(:,1:m) = Ft';
   end
end
     
count = count+np;

for l = 1:p
    Ft = cc;
    dintmp = ( (ardin(:, 1:2:2*m) <= l) & (ardin(:, 2:2:2*m) >=l)  );
    np = sum(sum(dintmp));
    Ft(dintmp') = theta(count:count+np-1);
    F(:,l*m+1:(l+1)*m) = Ft';
    count = count+np;
end

for l = 1:q
    Ft = cc;
    dintmp = (madin > 0);
    np = sum(sum(dintmp));
    Ft(dintmp') = theta(count:count+np-1);
    A(:,l*m+1:(l+1)*m) = Ft';
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
    count = count+np;
end

V = zeros(m,m);
if diagv
   V = diag(theta(count:count+m-1));
else
   cc = ones(m,m);
   V(triu(ones(m,m))) = theta(count:count+(m*(m+1)/2)-1);
   if ~EEOPTION(5), V( tril(cc,-1) ) = V( tril(cc,-1)' );
   else V = V'; end
end

% Force V positive definite
if ~EEOPTION(5)
   if min(eig(V)) <= 0
      Vu = cholp(V);
      V  = Vu'*Vu;
   end
else
   V = V*V';  
end