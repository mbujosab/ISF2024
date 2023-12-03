function [Phix,Hx,Ex,Qx,A1x,A2x,A3x,A4x,A5x] = sidenti(y,u,i,n,s)
%
% SIDENTI-  General algorithm for subspace-based system identification.
%           It can be applied to models with or without inputs
%           The model identified is not in canonical form
%           The operation is the same as that of SIDENT but this function
%           determines the optimal values of i and j
%
% For models with inputs the syntax is:
%     [Phix,Hx,Ex,Qx,A1x,A2x,A3x,A4x,A5x] = sidenti(y,u,ij,n,s)
% where:
%     y  >  Matrix of endogenous variables data
%     u  >  Matrix of exogenous variables (inputs) data
%     i  >  Row vector that contains the range of values for the dimension of the
%           past (i) and future (j) information subspaces. Assumes that i=j.
%           The optimality criterion consists in minimizing the variance of
%           one-step-ahead forecast errors
%     n  >  System dimension
%     s  >  Seasonal period (if any)
%     Phix,Hx,Ex,Qx, ... < ???
%
%  For models without inputs the syntax is:
%     [Phix,Hx,Ex,Qx,A1x,A2x,A3x,A4x,A5x] = sidenti(y,[],ij,n,s)
%
% 23/12/96

i = i(:);
r = size(u,2);

if nargin < 5, s = 1; end

fmax = inf;
for k=1:size(i,1)
   [Phi,H,E,Q,A1,A2,A3,A4,A5] = sident(y, u, i(k), n, s);
   if r
      [thet, di] = ss2thd(Phi, A1, E, H, A2, [], Q);
      f = pem(thet,di,[y u]);
   else
      [thet, di] = ss2thd(Phi, [], E, H, [], [], Q);
      f = pem(thet,di,y);
   end
   if f < fmax
      fmax = f;
      Phix = Phi;
      Hx   = H;
      Ex   = E;
      Qx   = Q;
      A1x = A1;
      A2x = A2;
      A3x = A3;
      A4x = A4;
      A5x = A5;
   end
end
