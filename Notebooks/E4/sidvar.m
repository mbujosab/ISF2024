function [Q,S,R] = sidvar(y, Phi, Gam, E, H, D, C, opt)
%
% SIDVAR -
% [Q,S,R] = sidvar(y, Phi, Gam, E, H, D, C, opt)
%
% opciones
%     0 : libre
%     1 : Q = R = S
%     2 : S = Q
%     3 : S = R
%     4 : S = []
%     5 : R = S = []
%     6 : Q = S = []
%  10+x : Q diagonal
%  20+x : R diagonal
%  30+x : Q y R diagonal 

if nargin < 8, opt = 0; end

diago = floor(opt/10);
opt  = rem(opt,10);

global SIDOPTION;

ext   = SIDOPTION(2);

if max(abs(eig(Phi))) >= 1 - eps & ext
   fprintf(1,'No debería usar matrices Block-Hankel extendidas en modelos no estacionarios');
end

n = size(Phi,1);
m = size(H,1);

O = zeros(m*(2*n+1),n);

O(1:m,:) = H;
HPhi = H*Phi;

for j=2:2*n+1
    O((j-1)*m+1:j*m,:) = HPhi;
    HPhi = HPhi*Phi;
end

[U S V] = svd(O(1:m*(n+1),:));
n2 = sum(diag(S) > eps);
if n2 < n
   n = n2;
   [U S V] = svd(O(1:m*(n+1),:));
end

Op0 = U(:,n+1:n+m)';
i = 2*n+1;
Op = zeros(m*(n+1),i*m);
for j=0:n
    Op(j*m+1:(j+1)*m,j*m+1:j*m+(n+1)*m) = Op0;
end

r = max(size(Gam,2), size(D,2));

Yi = blkhkel(y(:,1:m), i, 1, ext);
N = min(size(Yi,2), size(y,1));

if r

   Hu = zeros(i*m,i*r);

   for j=1:i
       Hu((j-1)*m+1:j*m,(j-1)*r+1:j*r) = D;
       Hu(j*m+1:i*m,(j-1)*r+1:j*r) = O(1:(i-j)*m,:)*Gam;
   end

   Ui = blkhkel(y(:,m+1:m+r), i, 1, ext);

   Yi = Yi - Hu*Ui;
end

Y = Op*Yi*Yi'*Op(1:m,:)'/N;
Y = Y(:);

sq = size(E,2);
sr = size(C,2);

im = i*m;
iq = i*sq;
ir = i*sr;
Hq = zeros(i*m,i*sq);

for j=1:i
    Hq(j*m+1:i*m,(j-1)*sq+1:j*sq) = O(1:(i-j)*m,:)*E;
end

Hq = Op*Hq;
Hq2 = Hq(1:m,:);
Hr = Op*kron(eye(i),C);
Hr2 = Hr(1:m,:);

Ii  = eye(i);
IIq = zeros(i*iq,sq);
Iq  = eye(sq);
for k=1:i
    IIq((k-1)*iq+1:k*iq,:) = kron(Iq,Ii(:,k));
end
IIq = kron(IIq,eye(sq));

IIr = zeros(i*ir,sr);
Ir  = eye(sr);
for k=1:i
    IIr((k-1)*ir+1:k*ir,:) = kron(Ir,Ii(:,k));
end
IIr = kron(IIr,eye(sr));

IIs = zeros(i*ir,sr);
for k=1:i
    IIs((k-1)*ir+1:k*ir,:) = kron(Ir,Ii(:,k));
end
IIs = kron(IIs,eye(sq));

Dp = reshape([1:sq*sr],sq,sr);
Dp = Dp';
Dp = Dp(:);

if diago == 1 | diago == 3
   Dq = [1:sq+1:sq*sq];
else
   Dq = [1:sq*sq]; 
end

if diago > 1
   Dr = [1:sr+1:sr*sr];
else
   Dr = [1:sr*sr]; 
end


if ~opt 
   X = [kron(Hq2,Hq)*IIq(:,Dq) kron(Hr2,Hr)*IIr(:,Dr) kron(Hr2,Hq)*IIs+kron(Hq2,Hr)*IIs(:,Dp)];
size(X)
elseif opt==1
   X = [kron(Hq2,Hq)*IIq(:,Dq)+kron(Hr2,Hr)*IIr(:,Dr)+kron(Hr2,Hq)*IIs(:,Dq)+kron(Hq2,Hr)*IIs(:,Dp(Dq))];
elseif opt==2
   X = [kron(Hq2,Hq)*IIq(:,Dq)+kron(Hr2,Hq)*IIs(:,Dq)+kron(Hq2,Hr)*IIs(:,Dp(Dq)) kron(Hr2,Hr)*IIr(:,Dr)];
elseif opt==3
   X = [kron(Hq2,Hq)*IIq(:,Dq) kron(Hr2,Hq)*IIs(:,Dr)+kron(Hq2,Hr)*IIs(:,Dp(Dr))+kron(Hr2,Hr)*IIr(:,Dr)];
elseif opt==4
   X = [kron(Hq2,Hq)*IIq(:,Dq) kron(Hr2,Hr)*IIr(:,Dr)];
elseif opt==5
   X = kron(Hq2,Hq)*IIq(:,Dq);
elseif opt==6
   X = kron(Hr2,Hr)*IIr(:,Dr);
end
svd(X)


Beta = X \ Y;

Q = [];
S = [];
R = [];

if opt < 6
   if diago == 1 | diago == 3
      Q = diag(Beta(1:sq));
      k = sq;
   else
      Q = reshape(Beta(1:sq*sq),sq,sq);
      k = sq*sq;
   end
   if opt == 1
      R = Q; S = Q;
   elseif opt < 5 
      if diago > 1
         R = diag(Beta(k+1:k+sr));
         k = k + sr;
      else
         R = reshape(Beta(k+1:k+sr*sr),sr,sr);
         k = k + sr*sr;
      end
      if ~opt
         S = reshape(Beta(k+1:size(Beta,1)),sq,sr);
      elseif opt==2
         S = Q;
      elseif opt==3
         S = R;
      end
   end
else
   if diago > 1
      R = diag(Beta);
   else
      R = reshape(Beta,sr,sr);
   end
end
