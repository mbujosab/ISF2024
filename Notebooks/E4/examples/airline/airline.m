e4init
load airline.dat
y=log(airline);

% Defines the nonstationary version of the airline model
% The parameters corresponding to the unit roots are constrained
[theta, din, lab] = arma2thd([-1],[-1],[0],[0],[0],12);
theta=[theta zeros(size(theta))];
theta(1,2)=1;
theta(2,2)=1;

% Computes preliminar estimates
sete4opt('econd','zero','vcond','idej','var','fac');
theta=e4preest(theta,din,y);

% ... and then ML estimates
[thopt,it,lval,g,h]=e4min('lffast', theta, '', din, y);
[std, corrm, varm, Im] = imod(thopt,din,y);
prtest(thopt,din,lab,y,it,lval,g,h,std,corrm);
N=size(airline,1);
[yfor,Bfor]=foremod(thopt,din,y,12);
Bfor=sqrt(Bfor);

% Computes the conditional forecasts and plots all the results
yobj = log(airline(N,1)*1.25); % End year target
yext = [y; NaN*ones(11,1); yobj];
[yhat Bhat] = fismiss(thopt,din,yext);

figure;
hold on
plot([exp(y((N-23):N,1));yfor*NaN],'k-')
plot([y((N-23):N,1)*NaN;exp(yfor)],'k--');
plot([y((N-23):N,1)*NaN;exp(yhat(N+1:N+12))],'ko');
plot([y((N-23):N,1)*NaN;exp(yfor+1.96*Bfor)],'k-.');
plot([y((N-23):N,1)*NaN;exp(yfor-1.96*Bfor)],'k-.');
grid on   % MBB
%whitebg('w'); % MBB
hold off
