function [dPhi,dGam,dE,dH,dD,dQ] = eedvper(theta, din, i, userf)
% EEDVPER  - Computes the partial derivatives of the matrices of a periodic
%    model with respect to the i-th parameter in theta.
%    [dPhi, dGam, dE, dH, dD, dQ] = eedvper(theta, din, i, userf)
% theta > parameter vector.
% din   > matrix which stores a description of the model dynamics.
% i     > index of the parameter in the denominator of the partial derivative.
% userf > the first row contains the name of the user function and
%         the second row contains the name of the function which computes
%         the derivatives of the user function (only in user models).
%
% 5/3/97
% Copyright (c)

mtipo=din(1,1);

if mtipo >= 1000
    if nargin < 4, e4error(13); end
    if size(userf,1) < 2, e4error(13); end
    [dPhi,dGam,dE,dH,dD,dQ] = eedvpero(theta,din,i, Phi, E, D, userf(2,:));
    return;
end

if fix(din(1,1)/100) ~= 3, e4error(14); end
din(1,1) = din(1,1) - 300;

m = din(1,2); s = din(1,9); r = din(1,3); nmod = din(1,8);
k = din(1:m+1:nmod*(m+1),5);
asg = din(size(din,1),1:s);
kmax = max(max(k),1);
F = zeros(m*nmod,m*(kmax+1)); A = zeros(m*nmod,m*(kmax+1));
G = zeros(m*nmod,r*(kmax+1)); V = zeros(m*nmod,m);
iA0 = zeros(m*nmod,m);
   
h = 1:m+1:nmod*(m+1);
j = 1:m:s*m;
np = [0; cumsum(din(1:m+1:nmod*(m+1),6))];

im = find(np >= i);
im = im(1)-1;

for l=1:nmod
    if im ~= l
       [F(j(l):j(l)+m-1,1:(k(l)+1)*m), A(j(l):j(l)+m-1,1:(k(l)+1)*m), V(j(l):j(l)+m-1,:), G(j(l):j(l)+m-1,1:(k(l)+1)*r)] = thd2fr(theta(np(l)+1:np(l+1),1), din(h(l):h(l)+m,:));
       iA0(j(l):j(l)+m-1,:) = pinv(A(j(l):j(l)+m-1,1:m));
    else
       [F(j(l):j(l)+m-1,1:(k(l)+1)*m), dF, A(j(l):j(l)+m-1,1:(k(l)+1)*m), dA, V(j(l):j(l)+m-1,:), dV, ...
         G(j(l):j(l)+m-1,1:(k(l)+1)*r), dG] = fr_dv(theta(np(l)+1:np(l+1),1), din(h(l):h(l)+m,:), i-np(l));
       iA0(j(l):j(l)+m-1,:) = pinv(A(j(l):j(l)+m-1,1:m));
       diA0 = -iA0(j(l):j(l)+m-1,:)*dA(:,1:m)*iA0(j(l):j(l)+m-1,:);
    end
end

k(~k) = ones(sum(~k),1);
ksum = sum(k(asg(1:s)));

Phi = zeros(ksum*m,m);
E   = zeros(ksum*m,m);

dPhi = zeros(ksum*m,kmax*m);
dE   = zeros(ksum*m,m);
dQ   = zeros(m*s,m);
dC   = zeros(m*s,m);

if r
   Gam = zeros(ksum*m,r);
   D   = zeros(m*s,r);
   dGam = zeros(ksum*m,r);
   dD   = zeros(m*s,r);
end

kr = 0; kr2 = 0;
kc = m;
kn = zeros(s,1);
Im = eye(m);

idm = find(asg == im);

for i=1:s
    for h=1:kmax
        idx = rem(i+h-1,s)+1;
        if h <= k(asg(idx))
        %
           kr = kr+m;
           Phi(kr-m+1:kr,:) = -F(j(asg(idx)):j(asg(idx))+m-1,h*m+1:(h+1)*m);
           E(kr-m+1:kr,:) =  A(j(asg(idx)):j(asg(idx))+m-1,h*m+1:(h+1)*m);
           if r
              Gam(kr-m+1:kr,:) =  G(j(asg(idx)):j(asg(idx))+m-1,h*r+1:(h+1)*r);
           end

           if asg(idx) == im
              dPhi(kr-m+1:kr,1:m) = -dF(:,h*m+1:(h+1)*m);
              dE(kr-m+1:kr,:)     =  dA(:,h*m+1:(h+1)*m);
              if r
                 dGam(kr-m+1:kr,:) = dG(:,h*r+1:(h+1)*r);
              end
           end
        end
    end
    kn(i) = kr - kr2;
    if asg(i) == im
       dE(kr2+1:kr,:) = dE(kr2+1:kr,:)*iA0(j(asg(i)):j(asg(i))+m-1,1:m) + E(kr2+1:kr,:)*diA0 + dPhi(kr2+1:kr,1:m);
       if r
          dD(j(i):j(i)+m-1,:) = dG(:,1:r);
          dGam(kr2+1:kr,:) = dGam(kr2+1:kr,:)+dPhi(kr2+1:kr,1:m)*G(j(asg(i)):j(asg(i))+m-1,1:r)+Phi(kr2+1:kr,:)*dG(:,1:r);
       end
       dQ1 = dA(:,1:m)*V(j(asg(i)):j(asg(i))+m-1,:)*A(j(asg(i)):j(asg(i))+m-1,1:m)';
       dQ(j(i):j(i)+m-1,:) = dQ1 + dQ1' + A(j(asg(i)):j(asg(i))+m-1,1:m)*dV*A(j(asg(i)):j(asg(i))+m-1,1:m)';
    else
       dE(kr2+1:kr,:) = dE(kr2+1:kr,:)*iA0(j(asg(i)):j(asg(i))+m-1,1:m) + dPhi(kr2+1:kr,1:m);
       if r
          dGam(kr2+1:kr,:) = dGam(kr2+1:kr,:)+dPhi(kr2+1:kr,1:m)*G(j(asg(i)):j(asg(i))+m-1,1:r);
       end
    end
    kr2 = kr;
end   

dH = zeros(s*m,max(kn));
dPhi = dPhi(:,1:max(kn));