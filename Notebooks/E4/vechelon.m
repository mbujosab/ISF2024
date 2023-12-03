function [Phi, Gam, E, H, D, C, Q, S, R] = vechelon(th, dn)

[H_D, type, m, r, s, n] = e4gthead(dn);

nvar= n/m+1;
theta = th;
%theta(nvar*m^2+1:nvar*m^2+m^2,1) = th(1:m^2,1)-th(nvar*m^2+1:nvar*m^2+m^2,1);
theta(nvar*m^2+1:nvar*m^2+m^2,1) = thetat(1:m^2,1);
din = tomod(dn);
[Phi, Gam, E, H, D, C, Q, S, R] = thd2ss(theta, din);
