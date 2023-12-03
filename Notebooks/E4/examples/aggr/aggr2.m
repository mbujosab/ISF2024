e4init
load ipi.dat; load vai.dat;
x = ipi; y = vai;
yagr = NaN*zeros(252,1);
yagr(12:12:252) = y;

% Análisis del indicador
d112lx = transdif(x,0,1,1,12);
plotsers(d112lx,-1,'d112 log indicador');
uidents(d112lx,39,'d112 log indicador');
descser(d112lx,'d112 log indicador');

[theta, din, lab] = arma2thd([0],[0 0],[0],[0],[0],12);
prtmod(theta, din, lab)
sete4opt('econd','zero','vcond','lyap','var','fac');
theta=e4preest(theta,din,d112lx);

[thopt,it,lval,g,h]=e4min('lffast', theta, '', din, d112lx);
[std, corrm, varm, Im] = imod(thopt, din, d112lx);
prtest(thopt,din,lab,d112lx,it,lval,g,h,std,corrm);

[xhat, phat, e] = fismod(thopt, din, d112lx);
plotsers(e,-1,'residuals of indicator');
uidents(e,39,'residuals of indicator');
descser(e,'residuals of indicator');

% prueba multivariante

phi=[-1 NaN; NaN -1];
sphi=[-1 NaN; NaN -1];
the=[0 NaN; NaN 0];
sthe=[0 NaN; NaN 0];
v=[0 0;0 0];
[theta, din, lab] = arma2thd([phi],[sphi],[the],[sthe],[v],12);
theta=[theta zeros(size(theta,1),1)];
theta(1:4,2)=ones(4,1);
prtmod(theta, din, lab);
z=[yagr x];
z=log(z);
sete4opt('econd','zero','vcond','idej','var','fac');
theta=e4preest(theta,din,z);

[thopt,it,lval,g,h]=e4min('lfmiss', theta,'', din, z);
[std, corrm, varm, Im] = imiss(thopt, din, z);
prtest(thopt,din,lab,z,it,lval,g,h,std,corrm);
