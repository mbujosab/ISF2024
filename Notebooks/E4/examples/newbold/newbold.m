clear
clc
e4init
% *** Series A. First read and transform the data
load seriesa.dat;
y=transdif(seriesa,0,1,1,12);

% Formulation and preestimation of the univariate model
[theta,din,lab]=arma2thd([],[],[0],[0],0,12);
sete4opt('vcond','lyap','econd','zero','var','fac');
theta=e4preest(theta,din,y);
prtmod(theta,din,lab);

% Optimize the likelihood function, compute the information matrix
% and print the results
[theta,it,lval,g,h]=e4min('lffast',theta,'',din,y);
[std,corrm,varm,Im]=imod(theta,din,y);
prtest(theta,din,lab,y,it,lval,g,h,std,corrm);

% *** Series B. Note that we do not define the THD model structure
% corresponding to series B and C, as it coincides with that of
% series A
load seriesb.dat;
y=transdif(seriesb,0,1,1,12);
theta=e4preest(theta,din,y);
prtmod(theta,din,lab);
[theta,it,lval,g,h]=e4min('lffast',theta,'',din,y);
[std,corrm,varm,Im]=imod(theta,din,y);
prtest(theta,din,lab,y,it,lval,g,h,std,corrm);

% *** Series C
load seriesc.dat;
y=transdif(seriesc,0,1,1,12);
theta=e4preest(theta,din,y);
prtmod(theta,din,lab);
[theta,it,lval,g,h]=e4min('lffast',theta,'',din,y);
[std,corrm,varm,Im]=imod(theta,din,y);
prtest(theta,din,lab,y,it,lval,g,h,std,corrm);

% Series D
load seriesd.dat;
seriesd=log10(seriesd);
y=transdif(seriesd,1,1,1,12);
[theta,din,lab]=arma2thd([0],[],[],[0],0,12);
theta=e4preest(theta,din,y);
prtmod(theta,din,lab);
[theta,it,lval,g,h]=e4min('lffast',theta,'',din,y);
[std,corrm,varm,Im]=imod(theta,din,y);
prtest(theta,din,lab,y,it,lval,g,h,std,corrm);

% Model diagnosis
[e,vt,wt,ve]=residual(theta,din,y);
titD='residuals from series D';
descser(e,titD);
plotsers(e,0,titD);
uidents(e,20,titD);

% Compute six months ahead forecasts
% Note that the nonstationary version of the model is used
phi=[-1+theta(1) -theta(1)];
sphi=-1;
sth=theta(2);
v=theta(3);
[thetaf,dinf,labf]=arma2thd([phi],[sphi],[],[sth],[v],12);
prtmod(thetaf,dinf,labf);
[yf,bf]=foremod(thetaf,dinf,seriesd,6);
% Forecasts of log sales
[(1:6)' yf bf]
% Forecasts of sales
10.^yf
