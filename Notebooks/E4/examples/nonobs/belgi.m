% Non observable components model
e4init
load belgi.dat;
y = log(belgi)*1000;
[theta,din,lab]=ss2thd([1 1; 0 1],[],[0;1],[1 0],[],[1],[0],[0],[0]);
% All parameters except the variances are constrained
theta = [theta ones(12,1)];
theta(10,2)=0;
theta(12,2)=0;

% Compute preliminar estimates
sete4opt('vcond','idej','econd','zero');
theta=e4preest(theta,din,y);

% ... and optimize the likelihood function.
[thopt,it,lval,g,h]=e4min('lffast', theta,'', din, y);
[std,corrm,varm,Im] = imod(thopt,din,y);
prtest(thopt,din,lab,y,it,lval,g,h,std,corrm);
NVR=thopt(10,1)/thopt(12,1);
disp(sprintf('NVR = %4.2f ', NVR));

% Smoothing to estimate the trend
[xhat,phat,e]=fismod(thopt,din,y);
trend=xhat(:,1);
deltat=transdif(trend,1,1);

figure;
hold on
plot(trend/1000,'k-')
plot(y/1000,'ko');
grid
%whitebg('w');
hold off
title('plot of log(PIB) versus trend');

plotsers(deltat,0,'changes of the trend');

% Model for the trend derivative: 
y1=transdif(trend,1,2);
[theta3,din3,lab3] = arma2thd([0 0],[],[],[],[0],1);

% Computes preliminar and ML estimates
sete4opt('econd','zero','vcond','lyap');
theta3=e4preest(theta3,din3,y1);
[thetan,it,lval,g,h]=e4min('lffast',theta3,'',din3, y1);
[std,corrm,varm,Im] = imod(thetan,din3,y1);
prtest(thetan,din3,lab3,y1,it,lval,g,h,std,corrm);

% Residual diagnostics
[ehat,vT,wT,vz1,vvT,vwT]=residual(thetan,din3,y1);
descser(ehat,'residuals of the model for the trend');
uidents(ehat,10,'residuals of the model for the trend');
plotsers(ehat,0,'residuals of the model for the trend');

