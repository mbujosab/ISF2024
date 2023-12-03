% Univariate modelling of the mink-muskrat series
e4init
load mink.dat;
% The data is already in logs
z1=transdif(mink(:,3),1,1);
z2=mink(:,2)-mean(mink(:,2));
z=[z1 z2(2:62)];

% Define the Muskrat skins model
[theta, din, lab]=arma2thd([0 0 0 0 0 0],[],[0],[],[0],1);

% Compute preliminar estimates
sete4opt('econd','zero','vcon','lyap');
theta=e4preest(theta,din,z1);

% Compute ML estimates, Information matrix and print the results
[thopt,it,l,g,H]=e4min('lffast',theta,'',din,z1);
[std,corrm,varm,Im]=imod(thopt,din,z1);
prtest(thopt,din,lab,z1,it,l,g,H,std,corrm);

% Define the Mink skins model
[theta, din, lab]=arma2thd([0 0 0 0],[],[],[],[0],1);

% ... and repeat again the previous process
theta=e4preest(theta,din,z2);
prtmod(theta, din, lab);
[thopt,it,l,g,H]=e4min('lffast',theta,'',din,z2);
[std,corrm,varm,Im]=imod(thopt,din,z2);
prtest(thopt,din,lab,z2,it,l,g,H,std,corrm);
