function [theta,din,lab] = varma2thd(F,sF,Th,sTh,Sigma)

% From subest4 to subest2 or EML

global e4option

m = size(F,1);
b = size(F,2);

fi=F(:,m+1:b);
th=Th(:,m+1:b);
fi(abs(F(:,m+1:b)./sF(:,m+1:b))<2)=NaN;
th(abs(Th(:,m+1:b)./sTh(:,m+1:b))<2)=NaN;

%th = reshape(th,m,b-m);

[theta, din, lab] = arma2thd([fi], [], [th], [], [Sigma]);
theta = [theta zeros(size(theta,1),1)];
theta(find(theta(:,1)==0),2)=1;