e4init

close all % Close figure window(s)

n=250; % sample size
p=12;  % monthly series

TTVP  = [1 1;1 0]; 
TV    = [500, 1, 5] 
TP    = p./[0:p/2];

[y,C] = simdhr (n,p,TTVP,TV,TP);

%[V,P,TVP,oar,results,MCNN,NVR] = autodhr (y,p); % more output information
%[V,P,TVP,oar] = autodhr (y,p,[],[],[],TTVP)     % search models like TTVP with only one (complex pair) AR root for each harmonic)
[V,P,TVP,oar] = autodhr (y,p,[],[],[],TTVP,1)   % do not identify (force the model TTVP)
%[V,P,TVP,oar] = autodhr (y,p)                   % minimun input information (search AR, RW, SRW or IRW for all componentes)

%[trend,season,cycle,irreg]  = dhrfilt (y,P,TVP,V,p,0,1); % more output information
[trend,season,cycle,irreg]   = dhrfilt (y,P,TVP,V,p);

figure(11), plot([irreg,      C(:,1)])            % estimated and true irregular components
figure(12), plot([trend(:,1), C(:,2)])            % estimated and true trends
figure(13), plot([season,     sum(C(:,3:end),2)]) % estimated and true seasonals
