function [Phi, Gam, E, H, D, C, Q, Phig, Gamg, Eg, Hg, Dg] = bekk(theta, din)

theta1 = theta(1:din(1,6),1);
theta2 = theta(din(1,6)+1:size(theta,1),1);
din = tomod(din);

n = din(1,2);

AB = reshape(theta2,n,2*n);

B = AB(:,n+1:2*n);
A = AB(:,1:n);

n2 = n*(n+1)/2;
k = tril(ones(n));
k = k(:);

V1 = kron(A,A);
V2 = kron(B,B);
V2 = V2(k,k);
V1 = V1(k,k) + V2;

[te,de] = arma2thd(-V1,[],-V2,[],zeros(n2,1),1);
[t,d] = garc2thd(theta1,din,te,de);
[Phi, Gam, E, H, D, C, Q, Phig, Gamg, Eg, Hg, Dg] = garch2ee(t, d);
