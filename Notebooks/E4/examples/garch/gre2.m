% Model: AR(4)+constrained ARCH(8). Implicit price deflator for GNP data.
e4init;
load gnpn.dat;
y = gnpn;
c = ones(size(y));

% Model for the mean
[thetay, diny, lab1] = arma2thd([0 0 0 0],[],[],[],[0],4,[0],1);

% Model for the conditional variance
[thetae, dine, lab2] = arma2thd([0 0 0 0 0 0 0 0],[],[],[],1,4);

% Defines the composite model
[theta, din, lab3] = garc2thd(thetay, diny, thetae, dine, lab1, lab2);

% Computes preliminar estimates
sete4opt('vcond','idej','econd','iu');
theta=e4preest(theta,din,[y c]);
prtmod(theta, din, lab3);
thetan = theta(1:7,1);
labn   = lab3(1:7,:);
dinn = touser(din,'arch8');

% ... and ML estimates. Note the user function in the call
% to e4min
[thopt,it,lval,g,h]=e4min('lfgarch',thetan,'',dinn,[y c]);
prtest(thopt,dinn,labn,[y c],it,lval,g,h,[],[]);
