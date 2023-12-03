% Wolfer sunspot data (corrected).
% Data input
e4init
load wolfercc.dat;
wolf10 = wolfercc/10;
wolf10=wolf10-mean(wolf10);

% Defines an ARMA(2,2) model and computes preliminar estimates
[theta1, din1, lab1] = arma2thd([0 0],[],[0 0],[],[0],1);
sete4opt('econd','zero','vcond','lyap','var','fac');
theta1=e4preest(theta1,din1,wolf10);

% ML estimation
[thopt1,it,lval,g,h]=e4min('lffast', theta1, '', din1, wolf10);
[std, corrm, varm, Im ] = imod(thopt1, din1, wolf10);
prtest(thopt1,din1,lab1,wolf10,it,lval,g,h,std,corrm);

% Computation of residuals and diagnosis
[ehat,vT,wT,vz1,vvT,vwT]=residual(thopt1,din1,wolf10);
descser(ehat,'sunspots data: residuals of ARMA(2,2)');
plotsers(ehat,0,'sunspots data: residuals of ARMA(2,2)');
uidents(ehat,25,'sunspots data: residuals of ARMA(2,2)');

% Defines an AR(2)+white noise and computes preliminar estimates
[th1, d1, l1] = arma2thd([0 0],[],[],[],[0],1);
[th2, d2, l2] = arma2thd([],[],[],[],[0],1);
[theta,din,lab] = stackthd(th1,d1,th2,d2,l1,l2);
[theta2,din2] = comp2thd(theta,din,lab);
lab2 = char(l1, 'Vu');
theta2=e4preest(theta2,din2,wolf10);

% Estimation
[thopt2,it,lval,g,h]=e4min('lffast', theta2, '', din2, wolf10);
[std, corrm, varm, Im ] = imod(thopt2, din2, wolf10);
prtest(thopt2,din2,lab2,wolf10,it,lval,g,h,std,corrm);
period = 2*pi/acos(-thopt2(1,1)/(2*sqrt(thopt2(2,1))));
disp(sprintf('period = %4.2f years', period));

% Diagnosis
[ehat,vT,wT,vz1,vvT,vwT]=residual(thopt2,din2,wolf10);
descser(ehat,'sunspots data: residuals of AR(2)+error');
plotsers(ehat,0,'sunspots data: residuals of AR(2)+error');
uidents(ehat,25,'sunspots data: residuals of AR(2)+error');

% Compute the 'clean' series
sunsp=wolf10-vT(:,2);
figure;clf;axis; % MBB
hold on
plot(sunsp,'k-')
plot(wolf10,'ko');
grid
%whitebg('w');
hold off
title('Plot of smoothed versus original sunspots (scaled deviations)');
