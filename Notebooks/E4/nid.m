function  [theta, din, lab, nopts, coinvec, urc, tab, innov] = nid(z, nvec, s, q, p); 

% NID - Computes a fast estimate of the McMillan index, the number of unit roots,
% the cointegrating rank and the cointegrating vector of a matrix of time series. 
%       [theta, din, lab, nopts, coinvec, urc, tab, innov] = nid(z, nvec, s, q, p);
% z         > matrix of time series
% nvec      > a vector of possible McMillan indices
% s         > seasonality (just for univariate models)
% q         > logical flag:
%               q=0,  ur results using urootm 
%               q=1, (default) ur results using urootm1 
% p         > logical flag:
%               p=0,  results are not printed
%               p=1, (default) they are 
% 
% theta-din < the theta-din estimates corresponding to the nopt-system
%             order.
% lab       < vector which contains the names for the parameters in theta
% nopts     < optimal McMillan index estimated for each criterion
% coinvec   < cointegrating vector estimates
% urc       < urc(1) is the number of unit roots and
%             urc(2) is the number of cointegrating rank
% tab       < results of the tests for each element of nvec
% innov     < Matrix of one step ahead forecast errors
% 30/04/04

%sidinit;

if nargin < 5, p = 1; end
if nargin < 4, q = 1; end
if nargin < 3, s = 1; end
if nargin < 2 | isempty(nvec), nvec = 0:1:2*m; end % La dimension de pasado es MUY importante para la estimacion 
                                                   % del orden del sistema en MUESTRAS CORTAS!!
[N m] = size(z);
nv = size(nvec,2); 
i = max(round(log(N)),max(nvec)+1); 

    %% Criterios basados en los VS: SVC1--Bauer con Omega_2, NIDC--New Criterion  
S1 = singval(z, i, 1, s);
SVC = nid1([N m], i, nvec, S1);     

    %% Criterios a posteriori (estiman modelos)  
ICs=[];
for k=1:nv;
    [Phi,H,E,Q,innov] = sident(z, [], i, nvec(k), s);
    Om = det((innov*innov')/N);
    L = -N*m/2*(1+log(2*pi))-N/2*log(Om);  
    k = 2*nvec(k)*m+m*(m+1)/2;
    LN2 = -2*L/N;
    AIC = LN2+2*k/N; 
    SBC = LN2+k*log(N)/N;
    HQ = LN2+k*2*log(log(N))/N; 
    ICs = [ICs AIC SBC HQ];
end
    
AkScHq=reshape(ICs,3,nv)';
if size(nvec)==[1 1]
   nopts = nvec*ones(1,5);
else
   [a1,a2] = min([AkScHq SVC]); nopts = nvec(a2);
end   
 
    %% Contraste de Tiao-Tsay (CHI^2)
tchi2 = sidang(z, [], i-1, s); 
pvchi2 = diag(tchi2); optchi = find(pvchi2>=.95);
if isempty(optchi), optchi2 = 0; else, optchi2 = max(optchi); end
if optchi2 > max(nvec), optchi2 = max(nvec);
elseif optchi2 < min(nvec), optchi2 = min(nvec);
end
if nvec(1,1)==0, pvchi2 = [NaN; pvchi2(1:nv-1)]; else, pvchi2 = pvchi2(nvec); end
    
tab = [nvec(1,1:nv)' AkScHq SVC pvchi2];
nopts = [nopts optchi2];
 
    %% Modelo en theta,din
nopt1 = tabulate(nopts+1); nopt = max(find(nopt1(:,2)==max(nopt1(:,2))))-1;  
[Phi,H,E,Q,inn] = sident(z, [], i, nopt, s);  innov = inn'; 
if isempty(Phi) 
   [theta, din, lab] = ss2thd(0, [], zeros(1,m), [1;zeros(m-1,1)], [], eye(m,m), Q); %% Seguro??
else
   C=eye(m,m);[theta, din, lab] = ss2thd(Phi, [], E, H, [], C, Q); 
end
Qeig = eig(Q);

    %% Uroots & Cointegration
if s == 1
    if m < 3
        S0 = singval(z, i, 0); % 0 = canonica, 1 = no canonica - V.S. No canonicos para el calculo de SVC's 
    else 
        S0 = singval(z, 2, 0);
    end
    if q
        ur = urootm1(N, i-1, S0);
    else
        ur = urootm(N, i-1, S0);
    end
    if ur > nopt, ur = nopt; end  %------%% ¡¡¡OJO CON ESTO!!! Seria mejor nopt = ur;
    [Phit, T]  = Bkdiag(Phi);
    Ht=H*T;
    c = m - ur;
    if c > 0 & ur > 0
        coinvec = Ht(1:c,1:ur)*inv(Ht(c+1:m,1:ur)); 
        coinvec = [eye(size(coinvec,1)) -coinvec];             
        urc = [ur c];
    else
        coinvec = 0;
        urc = [ur 0];
    end
else
    coinvec = [];
    urc = [];
end

        %% Output: Tablas
if p
   disp(' ');
   disp('**********Results from order estimation**********');
   disp(' ');
   disp('       n        AIC      SBC        HQ       SVC1      SVC2     PVCHI2');
   disp(' ');
   disp(tab);
   disp('*************************************************');
   disp('                nopts  ');
   disp(' ');
   disp('    AIC   SBC   HQ   SVC1  SVC2  CHI2(95%) ');
   disp(nopts);    
   disp('*************************************************');
   if s == 1
      disp(' Unit roots '),disp(ur);
  end
   if m ~= 1
      disp(' Cointegrating rank '),disp(urc(2));
      disp('*************************************************');
      disp(' ');
      disp(' Cointegrating vector'),disp(coinvec);
      disp('*************************************************');
   end
end
if any(Qeig)<0, disp('       Be careful, Q is not positive definite matrix!'); end
  
%% End function