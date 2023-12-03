function P = lyapns(Phi1, Phi2, Q)
% LYAPNS  - Resuelve una ecuación algebraica de Lyapunov
%           P = Phi1*P*Phi2' + Q no simétrica
%           P = lyapunov(Phi1, Phi2, Q)
% 6/3/96
% (C@@)

[U1, Phi1] = schur(Phi1);
[U1, Phi1] = rsf2csf(U1,Phi1);
[U2, Phi2] = schur(Phi2);
[U2, Phi2] = rsf2csf(U2,Phi2);

n = size(Phi1,1);
P = zeros(n);

Q = U1' * Q * U2;

for j = n:-1:1
    for i = n:-1:1
        d = 1 - Phi1(i,i)*Phi2(j,j)';
        suma = Phi1(i,i:n) * P(i:n,j:n) * Phi2(j,j:n)' + Q(i,j);
        if (abs(suma) < eps)
           P(i,j) = 0;
        elseif (abs(d) < eps)
           error('Sistema sin solución');
        else
           P(i,j) = suma / d;
 	end
    end
end

P = real(U1 * P * U2');

