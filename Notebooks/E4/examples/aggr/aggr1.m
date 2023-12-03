% Disaggregation of value added in industry
% First, using nonstationary models
e4init;
load ipi.dat; load vai.dat;
x = ipi; y = vai;
[theta,din, lab] = tf2thd([],[-1.268 .268],[],[], [27.780], ...
                   12, [6.705], [0 0 0 0 0 0 0 0 0 0 0 -.268]);

yagr = NaN*zeros(252,1);
yagr(12:12:252) = y;

sete4opt('econ','ml','var','fac');
[yhat, vyhat] = aggrmod(theta, din, [yagr x],12);
figure; 
plot(vyhat);
title('Variance of the monthly VAI')
plotsers(yhat,0,'monthly VAI');

% Disaggregation with the alternative model
[theta2,din2]=tf2thd([-1],[-.003 -.07],[],[], [2.928],...
               12, [6.705], [0 0 0 0 0 0 0 0 0 0 0 -.268]);
[yhat2, vyhat2] = aggrmod(theta2, din2, [yagr x],12);
figure; 
plot(vyhat2);
title('Variance of the monthly VAI')
plotsers(yhat2,0,'monthly VAI');

% Disaggregation using a stationary model
dyagr=transdif(yagr,1,0,1,12)
dx=transdif(x,1,0,1,12);
[dtheta,ddin]=arma2thd([],[-.268],[],[], [27.780], 12,[6.705], 1);
[dyhat,dvyhat]=aggrmod(dtheta,ddin,[dyagr dx],12);
plotsers(dyhat,0,'annual increments of monthly VAI');
