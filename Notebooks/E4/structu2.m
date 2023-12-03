clear; clc
e4init
sete4opt('econd','zero','vcond','idej','var','fac');

%Stochastic trend
[theta1,din1,lab1]=arma2thd([-2 1],[],[],[],[.2]);
theta1=[theta1 zeros(size(theta1,1),1)]; theta1(1,2)=1; theta1(2,2)=1;

%Stochastic cycle
%|Ct|      |cos(lambda)  sin(lambda)| |Ct-1|   |Kt|
%|_ | = ro*|                        | |_   | + |_ |
%|Ct|      |-sin(lambda) cos(lambda)| |Ct-1|   |Kt|
%                    _
% with var(Kt) = var(Kt)
% Requires a user function
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
[theta3,din3,lab3]=arma2thd(ones(1,11),[],[],[],[.4],12);
theta3=[theta3 ones(size(theta3,1),1)]; theta3(12,2)=0;

%Structural model: trend + cycle + seasonality + noise
[theta,din,lab]=stackthd(theta1,din1,theta2,din2,lab1,lab2);
[theta,din,lab]=stackthd(theta,din,theta3,din3,lab,lab3);
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