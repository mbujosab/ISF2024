function [theta2, hess] = subes1(theta, din, z, garch)
%
% March 2002
% Copyright (c) Jaime Terceiro, 1997


global E4OPTION
E4OPTION(9) = 1000; % Maxiter

if nargin < 3,  e4error(3); end

% innovations?
[H_D, type, m, r, s, n, np, userflag, userf, innov] = e4gthead(din);
if ~innov(1)
    variance = E4OPTION(5);
    E4OPTION(5) = 1; % Cholesky factors
end

[Phi, Gam, E, H, D, C, Q, S, R] = thd2ss(theta, din);
N = size(z,1);
n = size(Phi,1);
m = size(H,1);
r = max([size(Gam,2), size(D,2)]);
if size(z,2) ~= m+r, e4error(11); end

if any(abs(eig(Phi)) >= E4OPTION(12)), stat = 0; else stat = 1; end
if stat | innov(1), offs = 0; else, offs = 1; end

i = max(round(log(N)),ceil(n/m)+2+offs);

if N < 2*i*(m+r+1-stat)-1, j = max(fix((N-i*(m+r+1)+1)/(m+r+1)),ceil(n/m)+1); else j = i; end
ij = i+j;

if N < ij*(m+r+1-stat)-1
   i = ceil(n/m)+1+offs; j = i; ij = i+j;
   if N < i*(m+r+1-stat)+j*(r+1-stat)-1, e4error(35); end
end

missing = any(any(isnan(z(:,1:m))));
if missing
   z2 = z;
   for h=1:m
       k = find(isnan(z(:,h)));
       sk = size(k,1);
       if sk
          k2 = find(k-[min(0,k(1)-2);k(1:sk-1)] > 1);
       
          if size(k2,1) <= 1
             if k(1) == 1, z1 = z(k(sk)+1,h);
             elseif k(sk) == N, z1 = z(k(1)-1,h);
             else, z1 = (z(k(sk)+1,h) + z(k(1)-1,h))/2; end
             z2(k,h) = z1*ones(sk,1);
          else
             if k(1) == 1, z1 = z(k(k2(2)-1)+1,h);
             else, z1 = (z(k(k2(2)-1)+1,h) + z(k(1)-1,h))/2; end
             z2(k(k2(1)):k(k2(2)-1),h) = z1*ones(k(k2(2)-1)-k(k2(1))+1,1);

             for l=2:size(k2,1)-1
                 z2(k(k2(l)):k(k2(l+1)-1),h) = (z(k(k2(l))-1,h)+z(k(k2(l+1)-1)+1,h))/2*ones(k(k2(l+1)-1)-k(k2(l))+1,1);
             end
 
             if k(sk) == N, z1 = z(k(k2(l+1))-1,h);
             else, z1 = (z(k(sk)+1,h) + z(k(k2(l+1))-1,h))/2; end
             z2(k(k2(l+1)):k(sk),h) = z1*ones(k(sk)-k(k2(l+1))+1,1);
          end
       end
   end
end
   
if r
   Yi = blkhkel(z(:,1:m), ij, 1, stat);
   N = min(size(Yi,2), size(z,1));
   if missing
      K  = any(isnan(Yi));
      Yi = blkhkel(z2(:,1:m)*N/(N-sum(K)/2), ij, 1, stat);
      Ui = blkhkel(z(:,m+1:m+r)*N/(N-sum(K)/2), ij, 1, stat);
      Yi(:,K) = Yi(:,K)/2;
      Ui(:,K) = Ui(:,K)/2;
   else
      Ui = blkhkel(z(:,m+1:m+r), ij, 1, stat);
   end
   [Q, R] = qr([Ui(r*i+1:r*ij,:);Ui(1:r*i,:);Yi]', 0);
   ix = zeros(5,2);
   ix(:,1) = [1; j*r+1; ij*r+1; ij*r+i*m+1; ij*r+(i+1)*m+1];
   ix(:,2) = [ix(2:5,1)-1; ij*(m+r)];
   R = R'/sqrt(N);
else
   Yi = blkhkel(z, ij, 1, stat);  %%% HEY seas!!
   N = min(size(Yi,2), size(z,1));
   if missing
      K  = any(isnan(Yi));
      Yi = blkhkel(z2(:,1:m)*N/(N-sum(K)/2), ij, 1, stat);
      Yi(:,K) = Yi(:,K)/2;
   end
   [Q, R] = qr(Yi', 0);
   R = R'/sqrt(N); 
   ix = zeros(3,2);
   if stat | innov(1)
      ix(:,1) = [1; i*m+1; (i+1)*m+1];
   else 
      ix(:,1) = [m+1; i*m+1; (i+1)*m+1];
   end
   ix(:,2) = [ix(2:3,1)-1; ij*m];
end

if size(R,1) > size(R,2)
   R = [R zeros(size(R,1), size(R,1)-size(R,2))];
   pond = 0;
else
   if stat | innov(1), pond = 1; else, pond = 0; end
end

sth = size(theta,1);
if size(theta,2) == 1, theta2 = [theta zeros(sth,1)]; else, theta2 = theta; end

index = e4ds(theta,din);
index = ~(~index);

k = find(theta2(:,2)==0);

if size(k,1) > 0
   verbose = E4OPTION(10); 
   E4OPTION(10) = 0; % verbose no
   if ~stat & ~r 
       theta2(k,1) = -.5*ones(size(k,1),1);  % VEr porque con phi=-1.01 pbs?
       [f e] = lffast(theta2,din,z(:,1:m));
       covs = vech(cov(e));
   else
       covs = vech(cov(z(:,1:m)));
   end
      theta2(k,1) = -.5*ones(size(k,1),1);   % VEr porque con phi=-1.01 pbs?
 if innov(1)
     theta2(find(~theta2(:,2) & index(:,2)))= covs;
 else
     theta2(find(~theta2(:,2) & index(:,2)))= covs(1:m,1:m)*ones(size(theta2(find(~theta2(:,2) & index(:,2)))));
 end
      [theta2, iterat, lfopt, grad, hess]  = e4min('subes2', theta2, '', din, R, ix, [i j N innov(1) garch]);
%  theta2 = e4min('subes2', theta2, '', din, R, ix, [i j N innov(1) garch]);   
end
if size(theta,2) == 1, theta2 = theta2(:,1); else, theta2(:,2) = theta(:,2); end

E4OPTION(10) = verbose; 
if ~innov(1) E4OPTION(5) = variance;, end
