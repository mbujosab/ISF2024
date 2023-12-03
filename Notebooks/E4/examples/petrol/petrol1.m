% TF model of petrochemical consumption and industrial production
e4init
load petro.dat;
petro(:,1:2) = log(petro(:,1:2));
y=transdif(petro,1,1,1,4);

% Defines the structure of the transfer function
sar = [0 0]; ma = [0]; sma = [0]; v = [0];
w   = [0  NaN  NaN; 0 0 0];
[theta, din, lab] = tf2thd([],[sar],[ma],[sma],[v],4,[w],[]);

% Computes preliminar estimates
sete4opt('vcond','lyap','econd','ml','var','fac');
theta=e4preest(theta,din,y);
prtmod(theta,din,lab);

% ... and ML estimates
[thopt,it,lval,g,h]=e4min('lffast', theta, '', din, y);
[std,corrm,varm,Im]=imod(thopt,din,y);
prtest(thopt,din,lab,y,it,lval,g,h,std,corrm);

% Computes the period of the noise model
period = 2*pi/acos(-thopt(1,1)/(2*sqrt(thopt(2,1))));
disp(sprintf('period = %4.2f years', period));

% Validation
[ehat,vT,wT,vz1,vvT,vwT]=residual(thopt,din,y);
descser(ehat,'petro. consumption residuals');
uidents(ehat,10,'petro. consumption residuals');
plotsers(ehat,0,'petro. consumption residuals');

% Explicit period estimation
thetan = thopt;
thetan(2,1) = 3.7;
labn = char(lab(1,:),'Period',lab(3:size(lab,1),:) );
dinn = touser(din,'pcons1');

[theta,it,lval,g,h]=e4min('lffast',thetan,'',dinn,y);
prtest(theta,din,labn,y,it,lval,g,h);

% Constrained period estimation
thetan=theta(:,1);
thetan(2,1)=3;
thetan=[thetan, [0; 1; 0; 0; 0; 0; 0; 0; 0]];
dinn = touser(din,'pcons1');
[theta,it,lval,g,h]=e4min('lffast',thetan,'',dinn,y);
prtest(theta,din,labn,y,it,lval,g,h);
