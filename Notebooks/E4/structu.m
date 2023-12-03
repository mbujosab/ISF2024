e4init
cd C:\MATLAB\TOOLBOX\E4\expcode
%Stochastic trend
%Mt = Mt-1 + Bt-1 + Nt
%Bt = Bt-1 + St

Phi = [1 1;NaN 1];  H = [1 NaN]; Q = [.5;.25];
[theta1, din1, lab1] = ss2thd(Phi, [], [], H, [], [], Q, [], []);
theta1 = [theta1 [ones(size(theta1,1)-2,1);0;0]];

%Stochastic cycle
%|Ct|      |cos(lambda)  sin(lambda)| |Ct-1|   |Kt|
%|_ | = ro*|                        | |_   | + |_ |
%|Ct|      |-sin(lambda) cos(lambda)| |Ct-1|   |Kt|
%                    _
% with var(Kt) = var(Kt)
% Requires user function for reparametization

Phi = [.7 1;.25 1]; H = [1 NaN]; Q = [.35;.35];
[theta2, din2, lab2] = ss2thd(Phi, [], [], H, [], [], Q, [], []);
theta2 = [theta2 [0;0;1;1;1;0;1]];
din2 = touser(din2,'Scycle');
lab2(1,:) = ['RO' blanks(size(lab2,2)-2)];
lab2(2,:) = ['LAMBDA' blanks(size(lab2,2)-6)];
lab2(3,:) = ['IGNORE' blanks(size(lab2,2)-6)];
lab2(4,:) = lab2(3,:);
lab2(6,:) = ['SIGMA_K' blanks(size(lab2,2)-7)];
lab2(7,:) = lab2(3,:);

%Dummy variable seasonality
[theta3, din3, lab3] = arma2thd(ones(1,11),[],[],[],.4,12);
theta3 = [theta3 [ones(11,1);0]];

%Noise
[theta4,din4,lab4] = arma2thd([],[],[],[],.09,1);

%Structural model: trend + cycle + seasonality + noise
[tc,dc,labc] = stackthd(theta1,din1,theta2,din2,lab1,lab2);
[tc,dc,labc] = stackthd(tc,dc,theta3,din3,labc,lab3);
[tc,dc,labc] = stackthd(tc,dc,theta4,din4,labc,lab4);
[tc,dc,labc] = comp2thd(tc,dc,labc);

prtmod(tc,dc,labc);

%
y=simmod(tc,dc,200);
sete4opt('econd','zero','vcond','idej','var','var');
% tc=e4preest(tc,dc,y);
[thopt,it,fval,g,h]=e4min('lffast',tc,'',dc,y);
prtest(thopt,dc,labc,y,it,fval,g,h);

% Diagnosis
ehat=residual(thopt,dc,y);
tit='residuos';
descser(ehat,tit);
plotsers(ehat,-1,tit);
uidents(ehat,39,tit,12);

[trend,season,cycle,irreg]=e4trend(thopt,dc,y,0);

figure;
whitebg('w');
hold on
plot(trend,'k-')
%axis([0 550 4 8])
title('Tendencia');
xlabel('Obs #')
hold off

figure;
whitebg('w');
hold on
plot(cycle,'k-')
%axis([0 550 -2 2])
title('Ciclo');
xlabel('Obs #')
hold off

figure;
whitebg('w');
hold on
plot(season,'k-')
%axis([0 550 -2 2])
title('Estacionalidad');
xlabel('Obs #')
hold off

figure;
whitebg('w');
hold on
plot(irreg,'k-')
%axis([0 550 -2 2])
title('Componente irregular');
xlabel('Obs #')
hold off