clear; clc
e4init
sete4opt('econd','zero','vcond','idej','var','fac');

%Stochastic trend
[theta1,din1,lab1]=arma2thd([-2 1],[],[],[],[.2]);
theta1=[theta1 zeros(size(theta1,1),1)]; theta1(1,2)=1; theta1(2,2)=1;

%Stochastic cycle
[theta2,din2,lab2]=arma2thd([-.6 .4],[],[],[],[.3]);

%Dummy variable seasonality
[theta3,din3,lab3]=arma2thd(ones(1,11),[],[],[],[.4],12);
theta3=[theta3 ones(size(theta3,1),1)]; theta3(12,2)=0;

%Noise
[theta4,din4,lab4] = arma2thd([],[],[],[],.9,1);

%Structural model: trend + cycle + seasonality + noise
[theta,din,lab]=stackthd(theta1,din1,theta2,din2,lab1,lab2);
[theta,din,lab]=stackthd(theta,din,theta3,din3,lab,lab3);
[theta,din,lab]=stackthd(theta,din,theta4,din4,lab,lab4);
[theta,din,lab]=comp2thd(theta,din,lab);
prtmod(theta,din,lab);

%
y=simmod(theta,din,200);
[thopt,it,fval,g,h]=e4min('lffast',theta,'',din,y);
prtest(thopt,din,lab,y,it,fval,g,h);

% Diagnosis
ehat=residual(thopt,din,y);
tit='residuos';
descser(ehat,tit);
plotsers(ehat,-1,tit);
uidents(ehat,39,tit,12);

[trend,season,cycle,irreg]=e4trend(thopt,din,y,0);

figure;
whitebg('w');
hold on
plot(trend,'k-')
title('Trend');
xlabel('Obs #')

figure;
plot(cycle,'k-')
title('Cycle');
xlabel('Obs #')

figure;
plot(season,'k-')
title('Seasonal');
xlabel('Obs #')

figure;
plot(irreg,'k-')
title('Irregular');
xlabel('Obs #')
hold off