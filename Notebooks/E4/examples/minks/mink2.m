% Bivariate modelling of the mink-muskrat series
% Model (5.9) of Jenkins and Alavi (1981) pag. 39
e4init
load mink.dat;
z1=transdif(mink(:,3),1,1);
z2=mink(:,2)-mean(mink(:,2));
z=[z1 z2(2:62)];

% Define the parameter matrices and generate
% the THD representation
phi1= [  0   0;   0  0];
phi2= [  0   0;   0  0];
phi3= [NaN NaN; NaN  0];
sigma=[  0   0;   0  0];

[theta,din,lab]=arma2thd([phi1 phi2 phi3],[],[],[],sigma,1);

% Compute preliminar estimates
sete4opt('econd','zero','vcon','lyap');
theta=e4preest(theta,din,z);
prtmod(theta,din,lab);

% Compute ML estimates, Information matrix and print the results
[thopt,it,l,g,H]=e4min('lffast',theta,'',din,z);
[std,corrm,varm,Im]=imod(thopt,din,z);
prtest(thopt,din,lab,z,it,l,g,H,std,corrm);
% Validation
[ehat,vT,wT,vz1,vvT,vwT]=residual(thopt,din,z);
tit=char('muskrat residuals','mink residuals');
descser(ehat,tit);
midents(ehat,10,tit);
plotsers(ehat,0,tit);
