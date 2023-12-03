% AR(4) model with constant term.
% The data in gnpn is the implicit price deflator for US GNP
% See Greene, ch.19, pg. 576
e4init;
load gnpn.dat;
y = gnpn;
c = ones(size(y));

% Defines the AR(4) structure and computes preliminar estimates ...
[theta, din, lab] = arma2thd([0 0 0 0],[],[],[],0,4,0,1);
sete4opt('vcond','idej','econd','iu');
theta=e4preest(theta,din,[y c]);

% ... and then, computes ML estimates under homoskedasticity
[theta,it,lval,g,h]=e4min('lffast',theta,'', din, [y c]);
[std, corrm, varm, Im]= imod(theta, din, [y c]);
prtest(theta,din,lab,[y c],it,lval,g,h,std,corrm);

% Finally, computes the residuals and its squares, to detect
% ARCH effects
[ehat,vT,wT,vz1,vvT,vwT]=residual(theta,din,[y c]);
plotsers(ehat,0,'AR(4) residuals');
uidents(ehat,15,'AR(4) residuals');
uidents(ehat.^2,15,'AR(4) squared residuals');
