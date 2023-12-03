function [f, z1, x0] = lffast(theta, din, z)
% lffast - Fast evaluation of the exact likelihood function for any time-invariant SS model
%    [f, innov, ssvect] = lffast(theta, din, z, userf)
% Under these conditions, this function outputs the same result that lfmod with
% vcond = 'idej'. It is valid for stationary and/or nonstationary models.
% theta  > parameter vector.
% din    > matrix which stores a description of the model dynamics.
% z      > matrix of observable variables.
% userf  > contains the names of the user functions that provide the system
%         matrices (only for user models).
% f      < value of the likelihood function.
% innov  < (optional) stores the sequence of innovations.
% ssvect < (optional) stores the sequence of values of the state vector.
%
% 7/3/97

% Copyright (C) 1997 Jaime Terceiro
% 
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2, or (at your option)
% any later version.
% 
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details. 
% 
% You should have received a copy of the GNU General Public License
% along with this file.  If not, write to the Free Software Foundation,
% 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.

global E4OPTION

if nargin < 3,  e4error(3); end

saveinn = 0; if nargout > 1, saveinn = 1; end
scaleb  = E4OPTION(2);
econd   = E4OPTION(4);
zeps    = E4OPTION(15);
if econd == 2, MV = 1; else, MV = 0; end

[H_D, type, m, r, s, n, np, userflag, userf, innov] = e4gthead(din);
[Phi, Gam, E, H, D, C, Q, S, R] = thd2ss(theta, din);
n = size(z,1);
l = size(Phi,1);
m = size(H,1);
r = max([size(Gam,2), size(D,2)]);
if size(z,2) ~= m+r, e4error(11); end

if ~innov(1)
%
   iN = inv(C*R*C');
   ESCt = E*S*C';
   M = E*Q*E' - ESCt*iN*ESCt';
   if trace(M) < zeps
      P=zeros(l);
   else
      A = Phi - ESCt*iN*H;
      T = [ A' zeros(l); -M   eye(l) ];
      Y = [eye(l) H'*iN*H; zeros(l) A ];

      [W,d] = eig(T,Y);
      d = diag(d);
      [e,index] = sort(abs(d));
      WW = W(:,index(1:l));
      P = real(WW(l+1:2*l,:)*pinv(WW(1:l,:)));
   end
   
   Q = H*P*H'+C*R*C'
   if scaleb, U = cholp(Q, abs(Q));
   else       U = cholp(Q); end
   iU = eye(m)/U';
   E = (Phi*P*H' + ESCt)*iU'*iU
%
else
%
   Q = C*Q*C';
   E = E * pinv(C);

   if scaleb, U = cholp(Q, abs(Q));
   else       U = cholp(Q); end
   iU = eye(m)/U';
%
end

if r & (econd == 1 | econd == 4)
%
   if econd == 1
      u0 = mean(z(:,m+1:m+r))';
   else
      u0 = z(1,m+1:m+r)';
   end

   [x0, Sigm, iSigm, nonstat] =  djccl(Phi, E*Q*E', 0, Gam*u0);
%
else
   [x0, Sigm, iSigm, nonstat] =  djccl(Phi, E*Q*E', 0);
end

if ~nonstat, iSigm = pinv(Sigm); end

Phibb0 = eye(l);
WW  = zeros(l); WZ = zeros(l,1);
Phib = Phi - E*H;
iQ = iU'*iU;

for t = 1:l  % initial filter loop
%
    if r
       z1  = z(t,1:m)' - H*x0 - D*z(t,m+1:m+r)';
       x0  = Phi*x0 + Gam*z(t,m+1:m+r)' + E*z1;
    else
       z1  = z(t,:)' - H*x0;
       x0  = Phi*x0 + E*z1;
    end

    Y = H*Phibb0;
    WW  = WW + Y'*iQ*Y;
    WZ  = WZ + Y'*iQ*z1;
    Phibb0 = Phib*Phibb0;
%
end  % t

if any(isnan([WW WZ])) | any(isinf([WW WZ])), e4error(25); end

P1T = pinv(iSigm + WW);

% esta línea es la que está bien
% if ~MV, x1T = P1T*WZ; else, x1T = pinv(WW)*WZ; end
% Sigm = Phibb0*P1T*Phibb0';
% estas otras están provisionalmente para las pruebas bwana
if ~MV
   x1T = P1T*WZ;
   Sigm = Phibb0*P1T*Phibb0';
else
   x1T = pinv(WW)*WZ;
   Sigm = Phibb0*pinv(WW)*Phibb0';
end

iSigm = pinv(Sigm);
x0 = x0 + Phibb0*x1T;

WW  = zeros(l,l); WZ = zeros(l,1);
Phibb0 = iU*H;

if saveinn, N = eye(l); end

ff = 0.0;

z1 = zeros(m,n);

for t = l+1:n  % main loop
    if r
    %
       z1(:,t)  = z(t,1:m)' - H*x0 - D*z(t,m+1:m+r)';
       x0  = Phi*x0 + Gam*z(t,m+1:m+r)' + E*z1(:,t);
    %
    else
    %
       z1(:,t)  = z(t,:)' - H*x0;
       x0  = Phi*x0 + E*z1(:,t);
    %
    end
    
    WW  = WW + Phibb0'*Phibb0;
    WZ  = WZ + Phibb0'*iU*z1(:,t);
    Phibb0 = Phibb0*Phib;
%
end  % t

ff  =  2*(n-l)*sum(log(diag(U))) + sum(sum((iU*z1).^2));

if any(isnan([WW WZ])) | any(isinf([WW WZ])), e4error(25); end

[Ns S Ns] = svd(Sigm);
k = find(diag(S) < zeps);
if size(k,1) < l
   if size(k,1)
   %   
      N = Ns(:,1:k(1)-1);
      Sigm = S(1:k(1)-1,1:k(1)-1);
      WW = N'*WW*N;
      WZ = N'*WZ;
   %
   end
   M = cholp(Sigm);
   T = cholp(eye(size(M,1))+M*WW*M');
   ff = ff + 2*sum(log(diag(T)));
   ff = ff - sum((T'\(M*WZ)).^2);
end
%

f = 0.5*(ff + (n-l)*m*log(2*pi));

if saveinn, z1 = z1'; end
