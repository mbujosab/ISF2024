function [thetat,dint,ixmodes] = blockdiag(theta,din,toinnov)
% BLOCKDIAG - Obtains the block-diagonal representation of a SS model
%    [thetat,dint,ixmodes] = blockdiag(theta,din,toinnov)
% theta   > parameter vector.
% din     > matrix which stores a description of the model dynamics.
% toinnov > logical flag:
%             toinnov=1, transforms the model to obtain exact estimates
%             toinnov=0 (default) not transforms
% thetat, dint < theta-din values corresponding to the block-diagonal SS model
% ixmodes < indexes for the different states:
%           1 trend states
%           2 seasonal states
%           3 cycle states
% 18/7/02
% Copyright (c) Casals, Jerez y Sotoca, 1999

global E4OPTION
zeps = .0001;
one = .99;
if nargin < 3, toinnov = 0; end

[H_D, type, m, r, s, n, np, userflag, userf, innov, szpriv] = e4gthead(din);
[Phi, Gam, E, H, D, C, Q, S, R] = thd2ss(theta,din);

if toinnov & ~innov(1),[E,Q] = sstoinn(Phi, E, H, C, Q, S, R); end
[Phit, iT, T, eigen] = bkdiag(Phi);
Et = T*E;
Ht = H*iT;
if r, Gamt = T*Gam; else, Gamt = []; end

if E4OPTION(5)
   Q = cholp(Q)';
   R = cholp(R)';
end

if innov(1) | toinnov
   [thetat,dint] = ss2thd(Phit,Gamt,Et,Ht,D,eye(m),Q);
else
   [thetat,dint] = ss2thd(Phit,Gamt,Et,Ht,D,C,Q,S,R);
end

if n
   per = zeros(n,1);
   a = real(eigen);
   b = abs(imag(eigen));
   k = abs(a) < zeps & b > zeps;
   per(k) = 4*ones(sum(k),1);
   k2 = b < zeps;
   k3 = sign(a) < 0 & k2;
   per(k3) = 2*ones(sum(k3),1);
   k3 = ~(k | k2);
   if sum(k3)
      at = atan(b(k3)./a(k3));
      k2 = at < 0;
      at(k2) = at(k2) + pi;
      per(k3) = 2*pi./at;
   end

   ktrend = real(eigen) > one & ~per;
   sper = s./(1:(floor(s/2)));
   kseason = ~ones(n,1);

   for i=1:size(sper,2)
       kseason = kseason | (per - zeps < sper(i) & per + zeps > sper(i));
   end

   kcycle = ~(kseason | ktrend);
else
   ktrend = []; kcycle = []; kseason = [];
end

ixmodes = zeros(size(Phi,1),1);
ixmodes(ktrend) = ones(sum(ktrend),1);
ixmodes(kseason) = ones(sum(kseason),1)*2;
ixmodes(~ixmodes) = ones(sum(~ixmodes),1)*3;
