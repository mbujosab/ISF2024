function f = subes2(theta, din, R, ix, ij)
%
%   Algoritmo 
% 
%
% Marzo 2002
% Copyright (c) Jaime Terceiro, 1997

[Phi, Gam, E, H, D, C, Q, S, Ri] = thd2ss(theta, din);

n = size(Phi,1);
m = size(H,1);
r = max([size(Gam,2), size(D,2)]);

i = ij(1);
j = ij(2);
N = ij(3);
innov = ij(4);
garch = ij(5);

if ~innov
        [E,Q,U,iU,P0] = sstoinn(Phi, E, H, C, Q, S, Ri); 
else
	iC = pinv(C); E = E*iC;
    if garch    
        Q = eye(m,m);
    end
    P0 = lyapunov(Phi, E*Q*E', .9999);
end

O = [H; zeros((j-1)*m,n)];
Phi2 = Phi;

for k=1:j-1
    O(k*m+1:(k+1)*m,:) = H*Phi2;
    Phi2 = Phi2*Phi;
end

iO = pinv(O);
O1 = O(1:(j-1)*m,:);

for i=1:i   % variance equations of the KF in innovations form      
	B1  = H*P0*H' + Q;
    U = cholp(B1);
	K  = ((Phi*P0*H' + E*Q)/U)/U';
	Phib = Phi - K*H;	
	P0  = Phib*P0*Phib' + E*Q*E' + K*Q*K' - K*(E*Q)' - E*Q*K';
end

His = [[eye(m); O1*K] zeros(j*m,(j-1)*m)]; % E in stead of K
Phi2 = Phi;
for k=2:j
    His((k-1)*m+1:j*m,(k-1)*m+1:k*m) = His(1:(j-k+1)*m,1:m);
end

if r
        Hid = [[D;O1*Gam] zeros(j*m,(j-1)*r)];%
   for k=2:j
       Hid((k-1)*m+1:j*m,(k-1)*r+1:k*r) = Hid(1:(j-k+1)*m,1:r);
   end
          X = [iO*(R(ix(4,1):ix(5,2),ix(1,1):ix(3,2))-Hid*R(ix(1,1):ix(1,2),ix(1,1):ix(3,2))) zeros(n,m)];
        e = R(ix(4,1):ix(5,2),ix(1,1):ix(4,2))-O*X-Hid(1:j*m,1:j*r)*R(ix(1,1):ix(1,2),ix(1,1):ix(4,2));
        e = [e R(ix(4,1):ix(5,2),ix(5,1):ix(5,2))];
else
      X = [iO*R(ix(2,1):ix(3,2),ix(1,1):ix(1,2)) zeros(n,m)];
      e = R(ix(2,1):ix(3,2),ix(1,1):ix(2,2))-O*X;
      e = [e R(ix(2,1):ix(3,2),ix(3,1):ix(3,2))];
end
       Pi = O*P0*O';
       Qi = kron(eye(j),Q); %B1 in stead of Q
       BB = Pi + His(1:j*m,1:j*m)*Qi*His(1:j*m,1:j*m)';       
       f = log(det(BB))+trace(inv(BB)*e*e');