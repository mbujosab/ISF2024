% Model for the supply and demand of food from Kmenta (1986)

e4init
load food.dat
Q = food(:,1); P = food(:,2);   D = food(:,3);
F = food(:,4); A = food(:,5); cte = food(:,6);
z = [Q P cte D F A];

% 2SLS estimates
t = [ 95; -0.24;  0.31; ...
      49;  0.24;  0.5;  0.25; 3.1; 1.7; 4.6];
% Model formulation
F0 = [1 -t(2);
      1 -t(5)];
G0 = [t(1), t(3),     0,     0; ...
      t(4),    0 ,  t(6), t(7)];
V  = vech2m(t(8:10),2);
lab = char(    'alpha1','alpha2','alpha3');
lab = char(lab,'beta1', 'beta2', 'beta3', 'beta4');
lab = char(lab, 'v1', 'c12', 'v2');

[tdum,din] = str2thd([F0],[],[],[],V,1,[G0],4);
din = touser(din,'food2ss');

sete4opt('econd', 'ml', 'vcond', 'idej');
[p,iter,lnew,g,h] = e4min('lffast',t,'',din,z);
prtest(p,din,lab,z,iter,lnew,g,h);

[e, vT, wT, Ve, VvT, VwT]=residual(p,din,z,1);
tit=['demand residuals';
     'supply residuals'];
midents(e,5,tit);
