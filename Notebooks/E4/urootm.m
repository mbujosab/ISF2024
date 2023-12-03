function  [ur] = urootm(N, i, S0); 

% Maximo tamaño del criterio en las simulaciones (a)
% [ur] = urootm(N, i, S0);

j = size(S0,1);

% Matriz de estimaciones de la funcion de penalizacion
C = [exp(.01)*N^(-.44)*i^(-.05) NaN;
     exp(.75)*N^(-.43) exp(.75)*N^(-.43);% PROBAR ADEMAS CON i!!
    -.30532+.03968*N-.00065*N^2+.0000033*N^3 exp(.78648)*N^(-.32804)*i^(-.2259); % exp(2.0482)*N^(-.6219)*i^(-.1439)
    -.634739+.044421*N-.000603*N^2+.00000267*N^3 exp(1.589335)*N^(-.436601)*i^(-.364576);%R4_90
    -.316565+.031558*N-.000315*N^2+.00000105*N^3-.075574*i exp(1.312803)*N^(-.38284)*i^(-.280176)];%-.326256+.030219*N-.0003*N^2+.000000989*N^3-.070474*i exp(1.122183)*N^(-.358629)*i^(-.300745)];

Pos = zeros(j,1);

if j>4
Pos(1,1) = 1-S0(1)-C(1,1);   
   for h = 2:4
       if N<88
       Pos(h,1) = 1-S0(h)^2-C(h,1);
       else
       Pos(h,1) = 1-S0(h)^2-C(h,2);    
       end
   end
   for h = 5:j
       if N<121 % OK!!
       Pos(h,1) = 1-S0(h)^2-C(5,1); 
       else
       Pos(h,1) = 1-S0(h)^2-C(5,2);    
       end
   end
else 
Pos(1,1) = 1-S0(1)-C(1,1);   
   for h = 2:j
       if N<88
       Pos(h,1) = 1-S0(h)^2-C(h,1);
       else
       Pos(h,1) = 1-S0(h)^2-C(h,2);    
       end
   end
end
    
   
ur = sum(Pos<0); %Antes Pos(1:m)