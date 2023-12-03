% Model: AR(4)+GARCH(1,1). Implicit price deflator for GNP data.
e4init;
load gnpn.dat;
y = gnpn;
c = ones(size(y));
% Model for the mean
[thetay, diny, lab1] = arma2thd([0 0 0 0],[],[],[],[0],4,[0],1);
% Model for the conditional variance
[thetae, dine, lab2] = arma2thd([0],[],[0],[],1,4);
% Composite model
[theta, din, lab3] = garc2thd(thetay,diny,thetae,dine,lab1,lab2);

% Computes preliminar estimates
sete4opt('vcond','idej','econd','iu');
theta=e4preest(theta,din,[y c]);
prtmod(theta, din, lab3);

% ... and ML estimates
[thopt,it,lval,g,h]=e4min('lfgarch',theta,'', din, [y c]);
[std,corrm,varm,Im]=igarch(thopt,din,[y c]);
prtest(thopt,din,lab3,[y c],it,lval,g,h,std,corrm);

% Validation. Note that in this case residual.m only returns
% the outputs 'ehat' and 'vz1'
[ehat,vT,wT,vz1]=residual(thopt,din,[y c]);
figure; % whitebg('w'); % MBB
plot(vz1);
title('Conditional variance')
plotsers(ehat,0,'original residuals');
stdres=ehat./sqrt(vz1);
plotsers(stdres,0,'standardized residuals');
uidents(stdres,15,'standardized residuals');
uidents(stdres.^2,15,'standardized squared residuals');
