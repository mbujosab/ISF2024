function  [ur] = urootm1(N, i, S0); 

% Mayor potencia del criterio en las simulaciones (b): POR DEFECTO en NID
% [ur] = urootm1(N, i, S0);

j = size(S0,1);

% Matriz de estimaciones de la funcion de penalizacion
C = [exp(.60)*N^(-.50)*i^(-.10) exp(.60)*N^(-.50)*i^(-.10);
     exp(.50)*N^(-.43) exp(.50)*N^(-.43);
    -.3530+.03577*N-.00059*N^2+.0000030*N^3 exp(.18818)*N^(-.2854)*i^(-.17157);
    -.62117+.040609*N-.000556*N^2+.00000247*N^3 exp(1.5570)*N^(-.4694)*i^(-.4179);
    -.366001+.029814*N-.000301*N^2+.00000101*N^3-.063056*i exp(1.131294)*N^(-.360703)*i^(-.378167)];
 
Pos = zeros(j,1);  
if j>4
   for h = 1:4
       if N<88
       Pos(h,1) = 1-S0(h)^2-C(h,1);
       else
       Pos(h,1) = 1-S0(h)^2-C(h,2);    
       end
   end
   for h = 5:j
       if N<121
       Pos(h,1) = 1-S0(h)^2-C(5,1);
       else
       Pos(h,1) = 1-S0(h)^2-C(5,2);    
       end
   end
else 
   for h = 1:j
       if N<88
       Pos(h,1) = 1-S0(h)^2-C(h,1);
       else
       Pos(h,1) = 1-S0(h)^2-C(h,2);    
       end
   end
end
 
ur = sum(Pos<0); %Antes Pos(1:m)