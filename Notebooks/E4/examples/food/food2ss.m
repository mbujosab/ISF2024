function [Phi, Gam, E, H, D, C, Q, S, R] = food2ss(th, din)
% Returns the SS representation of Kmenta model 
%  th(1:7)  contains the structural parameters
%  th(8:10) contains the lower triangle of the noise covariance matrix
t = th(:,1);
F0 = [1 -th(2);
      1 -th(5)];
G0 = [th(1) th(3)   0      0;
      th(4)  0   th(6) th(7)];
V  = vech2m(th(8:10),2);
% THD formulation of structural model
iF0 = diag(1./diag(F0)); % Normalization
[theta,din] = str2thd(iF0*F0,[],[],[],iF0*V*iF0',1,iF0*G0,4);
% SS conversion
[Phi, Gam, E, H, D, C, Q, S, R] = thd2sss(theta,din);
