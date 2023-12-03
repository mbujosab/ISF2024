function [theta,din,lab,ikopt,std,innov] = subest4(z,u,n,CN);

% SUBEST4 - Computes a fast estimate of the parameters in a VARMA echelon form
% [nonseasonal multivariate series]
% [theta,din,lab,ikopt,std,innov] = subest4(z,u,n,CN) 
%
% z         > matrix of time series
% u         > matrix of input series [],
% n         > system order
% CN        > criterion used to discriminate Kronecker Indexes
%               CN = 0, (default) AIC
%               CN = 1, SBC
%               CN = 2, HQ
% theta     < vector of estimates
% din       < vector which stores a description of the model dynamics
% lab       < vector which contains the names for the parameters in theta
% ikopt     < optimal Kronecker Indexes
% std       < vector which contains the standard errors of the parameters in theta
% innov     < Matrix of one-step-ahead forecast errors

global E4OPTION

if nargin < 4, CN = 0; end
N = size(z,1);
m = size(z,2);
r = size(u,2);
i = max(round(log(N)),n+2);       
Crits = [];

if r 
    if m == 1 & n ~= 0
        ikopt = n;
        [F,sF,Th,sTh,Sigma,G,Gs,innov] = varmaech(z, u, [i], ikopt);
    elseif n == 0
        [F,sF,Th,sTh,Sigma,G,Gs,innov] = varmaech(z, u, [i], n);
    else
        FF = (n+1)*ones(m,1);
        comb = fullfact(FF)-1;      % Estamos restringiendo la suma de IK a n.
        ik = comb(find(sum(comb')==n),:);
        for j = 1:size(ik,1)
            [Phi,sPhi,H,sH,E,sE,Q,Gam,sGam,D,sD,innov] = sidechel(z, u, i, ik(j,:)'); % Forma canonica Luemberger: SS- innov
            Om = det((innov*innov')/N);
            L = -N*m/2*(1+log(2*pi))-N/2*log(Om);    
            LN2 = -2*L/N;
            [F1,Th1] = kron2str(ik(j,:)');
            k=sum(sum(F1==2))+sum(sum(Th1==2));         % ¡¡No tiene la forma tipica pq depende de la estructura de ik!! 
            if CN == 0                
                Crit = LN2+2*k/N;                        
            elseif CN == 1
                Crit = LN2+k*log(N)/N;
            else
                Crit = LN2+k*2*log(log(N))/N;
            end
            Crits = [Crits Crit];
        end
        Crits = Crits';
        [s1,opt] = min(Crits); ikopt = ik(opt,:)';
        [F,sF,Th,sTh,Sigma,G,sG,innov] = varmaech(z, u, [i], ikopt);
    end
else
    if m == 1 & n ~= 0
        ikopt = n;
        [F,sF,Th,sTh,Sigma,G] = varmaech(z,[], [i], ikopt);
    elseif n == 0
        [F,sF,Th,sTh,Sigma,G] = varmaech(z,[], [round(log(N))], n);
    else
        FF = (n+1)*ones(m,1);
        comb = fullfact(FF)-1;      % Estamos restringiendo la suma de IK a n.
        ik = comb(find(sum(comb')==n),:);
        for j = 1:size(ik,1)
            [Phi,sPhi,H,sH,E,sE,Q,innov] = sidechst(z, [i], ik(j,:)');
            Om = det((innov*innov')/N);
            L = -N*m/2*(1+log(2*pi))-N/2*log(Om);     
            LN2 = -2*L/N;
            [F1,Th1] = kron2str(ik(j,:)');
            k=sum(sum(F1==2))+sum(sum(Th1==2));         % ¡¡No tiene la forma tipica pq depende de la estructura de ik!! 
            if CN == 0                
                Crit = LN2+2*k/N;                    
            elseif CN == 1
                Crit = LN2+k*log(N)/N;
            else
                Crit = LN2+k*2*log(log(N))/N;
            end
            Crits = [Crits Crit];
        end
        Crits = Crits';
        [s1,opt] = min(Crits); ikopt = ik(opt,:)';       
        [F,sF,Th,sTh,Sigma,innov] = varmaech(z,[], [i], ikopt);
    end
end 
[phi1,theta1,ge1]=kron2str(ikopt,r);
F(find(phi1==0))=NaN;
sF(1:m,1:m)=eye(m)+sF(1:m,1:m); sF(find(phi1==0))=NaN;
Th=Th(:,m+1:size(Th,2)); sTh=sTh(:,m+1:size(sTh,2));
Th(find(theta1==0))=NaN; sTh(find(theta1==0))=NaN;
G(find(ge1==0))=NaN; sG(find(ge1==0))=NaN;
[theta, din, lab] = ech2thd(ikopt, F, Th, Sigma, 1, G, r);
std = ech2thd(ikopt, sF, sTh, zeros(m,m), 1, sG, r);
innov = innov';
