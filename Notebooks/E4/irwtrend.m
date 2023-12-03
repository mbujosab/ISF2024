function [trend,dtrend,err,sigma2,pmat]=irwtrend(y,NVR)
% IRWTREND - extract the trend, derivative of trend and irregular component
%            using a Integrated Random Walk process 
% [trend,dtrend,err,sigma2]=irwtrend(y,nvr)
%     trend   <  Matrix of smoothed trend components (one column per series)
%     dtrend  <  Matrix of trend derivatives (one column per series)
%     err     <  Matrix of irregular components (one column per series)
%     sigma2  <  Matrix of estimated variances (one column per series)
%     pmat    <  Matrix of covariances of smoothed components (two columns per series)
%     y       >  Matrix of data to be smoothed. Each series must be in one column.
%     nvr     >  Noise-Variance-Ratio. If no NVR is specified the defaul HP value
%                1/1600 is assumed.
% The model is:
%   T(t+1)=T(t)+D(t)
%   D(t+1)=D(t)+nu(t)
%   y(t)=T(t)  +e(t)   
% Where T(t) is the trend in t, D(t) is the trend derivative in t.
% The variables nu(t) and e(t) are assumed to be zero-mean, uncorrelated
% and homoskedastic, with variances sig2nu and sig2e, respectively, and
% such that sig2nu/sig2e=nvr
%

global E4OPTION

% Saves old E4OPTION before choosing adequate initial conditions
oldeeopt=E4OPTION;
sete4opt('vcond','idej','econd','zero');

% If no NVR is chosen, select the default HP value for quarterly data
% For monthly data choose NVR=1/14400
if nargin == 1
  NVR=1/1600;
end

no_series=size(y,2);
trend=[]; dtrend=[]; err=[]; sigma2=[]; pmat=[];

for i=1:no_series
  sig2z=cov(transdif(y(:,1),1,2));
  sig2nu=sig2z/(1+6/NVR);
  sig2e=sig2nu/NVR;

  [theta,din,lab]=ss2thd([1 1; 0 1],[],[0;1],[1 0],[],[1],[sig2nu],[0],[sig2e]);

  [xhat,phat,e]=fismod(theta,din,y(:,i));
  sigma2=[sigma2 [sig2nu;sig2e]];
  trend=[trend xhat(:,1)];
  dtrend=[dtrend xhat(:,2)];
  err=[err e];
  pmat=[pmat phat];
end

% Restore old E4OPTION value and return
E4OPTION=oldeeopt;
