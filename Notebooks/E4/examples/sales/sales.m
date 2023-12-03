% Interpolation example
load sales.dat
y = sales;
N = size(y,1);
yobj = y(N)*1.11; % End year target

[theta, din] = arma2thd([-1],[-1],[-.8],[-.6],[1],12); 
sete4opt('econ','ml','vcon','zero');
% Compute twelve month forecasts
[yfor, Bfor] = foremod(theta, din, y, 12);
Bfor=sqrt(Bfor);

% Create an augmented series with eleven NANs and the end year target,
% interpolate the growth path
yext = [y; NaN*ones(11,1); yobj];
[yhat Bhat] = fismiss(theta,din,yext);

% ... and plot the results
plot(1:12,yfor,'kx',1:12,yhat(N+1:N+12),'ko');
hold on
plot(1:12,yfor+1.96*Bfor,'k-.',1:12,yfor-1.96*Bfor,'k-.');
%whitebg('w');
