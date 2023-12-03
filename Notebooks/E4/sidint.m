function [Phi,sPhi,H,sH,E,sE,Qr,D,sD,innov] = sidint(y, u, i, ki, s)
%
% SIDINT -   Algoritmo para modelos de intervención
%
%     [Phi,sPhi,H,sH,E,sE,Q,D,sD,innov] = sidint(y, u, ij, ki, s)
%
% 02/01/97


global SIDOPTION

if nargin < 5, s = 1; end

exact = ~SIDOPTION(1);
ext   = SIDOPTION(2);
canon = SIDOPTION(3);
pond  = SIDOPTION(4);
oest  = SIDOPTION(5);

if size(i,1) == 1, j = i(1,1); else, j = i(2,1); end
i = i(1,1);

n = sum(ki);

if n > i | n > j, siderror(5); end

m = size(y,2);
mb = sum(ki > 0);
r = size(u,2);

Yi = blkhkel(y, i+j, s, ext);
Ui = blkhkel(u, i+j, s, ext);
N = min(size(Yi,2), size(y,1));

[Q, R] = qr([Ui;Yi]', 0);
R = R';
ix = zeros(6,2);
ix2= zeros(6,2);
ix(:,1) = [1; i*r+1; (i+1)*r+1; (j+i)*r+1; (j+i)*r+i*m+1; (j+i)*r+(i+1)*m+1];
ix(:,2) = [ix(2:6,1)-1; (i+j)*(m+r)];
ix2(:,1) = [0; 0; 0; 1; i*m+1; (i+1)*m+1];
ix2(:,2) = [ix2(2:6,1)-1; (i+j)*m];


% step 2
D = R(ix(5,1):ix(5,2),ix(1,1):ix(2,2))*R(ix(2,1):ix(2,2),ix(1,1):ix(2,2))'* ...
    inv(R(ix(2,1):ix(2,2),ix(1,1):ix(2,2))*R(ix(2,1):ix(2,2),ix(1,1):ix(2,2))');

tau = inf;
iter = 0;
tol = 0.0001;
Ij = eye(j);
II = zeros(j*j*r,r);
Ir = eye(r);
for k=1:j
    II((k-1)*j*r+1:k*j*r,:) = kron(Ir,Ij(:,k));
end
II = kron(II,eye(m));

R5514 = R(ix(5,1):ix(5,2),ix(1,1):ix(4,2));
R2214 = R(ix(2,1):ix(2,2),ix(1,1):ix(4,2));
R4413 = R(ix(4,1):ix(4,2),ix(1,1):ix(3,2));
R4414 = R(ix(4,1):ix(4,2),ix(1,1):ix(4,2));
R5613 = R(ix(5,1):ix(6,2),ix(1,1):ix(3,2));
R1113 = R(ix(1,1):ix(1,2),ix(1,1):ix(3,2));
R2313 = R(ix(2,1):ix(3,2),ix(1,1):ix(3,2));


while tau > tol
%
  Dold = D;
  iter = iter+1;

% step 3
     R2 = R(ix(4,1):ix(6,2),ix(1,1):ix(6,2)) - kron(eye(i+j),D)*R(ix(1,1):ix(3,2),ix(1,1):ix(6,2));
     R4414r = R2(ix2(4,1):ix2(4,2),ix(1,1):ix(4,2));
     R5614r = R2(ix2(5,1):ix2(6,2),ix(1,1):ix(4,2));

% step 4
     if pond == 0
         Om = eye(j*m);
     elseif pond == 1
         R5616r = R2(ix2(5,1):ix2(6,2),ix(1,1):ix(6,2));
         Om = (R5616r*R5616r' - R5614r*R4414r'*inv(R4414r*R4414r')*R4414r*R5614r')/N;
     else
         R5616r = R2(ix2(5,1):ix2(6,2),ix(1,1):ix(6,2));
         Om = R5616r*R5616r';
     end

     sOm = sqrtm(Om);
     isOm = inv(sOm);

% step 5
     iRR = (R4414r*R4414r')^(-.5);
     [U S V] = svd(isOm*R5614r*R4414r'*iRR);
     AB = U(:,1:n)*S(1:n,1:n)*V(:,1:n)'*iRR;
% step 7

% step 8
     Y  = isOm*R5613-AB*R4413;
     X  = (kron(R2313',isOm)-kron(R1113',AB))*II;
     D = reshape(X \ Y(:), m,r);
     tau = norm(D-Dold)/norm(Dold);
%
end  % of while

O = sOm * U(:,1:n);
iO = U(:,1:n)'*isOm;
X = iO*R5614r*R4414r'*inv(R4414r*R4414r')*R4414r;
V = R5514 / [X;R2214];
H = V(:,1:n);
D = V(:,n+1:r+n);

% step 9

R4515r = R2(ix2(4,1):ix2(5,2),ix(1,1):ix(5,2));
R6615r = R2(ix2(6,1):ix2(6,2),ix(1,1):ix(5,2));
R5515 = R(ix(5,1):ix(5,2),ix(1,1):ix(5,2));
R5555 = R(ix(5,1):ix(5,2),ix(5,1):ix(5,2));
R6666r = R2(ix2(6,1):ix2(6,2),ix(6,1):ix(6,2));
R6615 = R(ix(6,1):ix(6,2),ix(1,1):ix(5,2));
R2215 = R(ix(2,1):ix(2,2),ix(1,1):ix(5,2));

%  La estimación de D no depende del sistema de coordenadas del vector de estado
ZR = [R5514-V*[X;R2214] R5555];
VV = diag(inv([X;R2214]*[X;R2214]'))';
sD = sqrt(diag(ZR*ZR'/N)*VV(n*m+1:(n+r)*m));

if oest == 0
   if pond == 0
      Om1 = eye((j-1)*m);
   elseif pond == 1
      R6616r = R2(ix2(6,1):ix2(6,2),ix(1,1):ix(6,2));
      Om1 = (R6616r*R6616r' - R6615r*R4515r'*inv(R4515r*R4515r')*R4515r*R6615r')/N;
   else
      R6616r = R2(ix2(6,1):ix2(6,2),ix(1,1):ix(6,2));
      Om1 = R6616r*R6616r';
   end
   sOm1 = sqrtm(Om1);
   isOm1 = inv(sOm1);
end

% step 6
if oest == 0
   [U1 S1 V1] = svd(isOm1*R6615r*R4515r'*(R4515r*R4515r')^(-.5));
   O1 = sOm1 * U1(:,1:n);
   iO1 = pinv(O(1:(j-1)*m,:))*O1*U1(:,1:n)'*isOm1;
else
   iO1  = pinv(O(1:(j-1)*m,:));
end

X1 = iO1*R6615r*R4515r'*inv(R4515r*R4515r')*R4515r;
X = [X zeros(n,m)];
Phi = X1 / X;
% solo aprox
%ZR = R5515-V*[X;R2215]
%Qr= ZR*ZR'/N
%E = X1*ZR'*inv(Qr)/N;

[Phi, H, T] = echelon(ki, Phi, H);
iT = inv(T);
[strH, strPhi] = echelstr(ki);
strH = strH(:);
strPhi = strPhi(:);

P = inv(chol(R5555*R5555')');
X  = iT*X;
Y = R5515;
ik = find(ki > 0);
YY = Y;
YY(ik,:) = Y(ik,:) - X(1:mb,:);
YY = P*YY;
X2 = kron(X(1:mb,:)',P);
idx = find(strH > 0);
sidx = size(idx,1);
X2 = [X2(:,idx) kron(R2215',P)];
V = X2 \ YY(:);

H1 = zeros(m*mb,1);
H1(idx) = V(1:sidx);
H1 = reshape(H1,m,mb);
D  = reshape(V(sidx+1:size(V,1)),m,r);
H = eye(m);
H = H1 + H(:,ik);
err = [YY(:)-X2*V];
VV  = sqrt(diag((err'*err + m)/(m*N)*inv(X2'*X2)));
H1(idx) = VV(1:sidx);
sH = reshape(H1,m,mb);

Zt = Y-H*X(1:mb,:)-D*R2215;

if n > mb
   Y = iT*iO1*R6615r - Phi(:,mb+1:n)*X(mb+1:n,:);
else
   Y = iT*iO1*R6615r;
end

P = iT*iO1*R6666r;
P = inv(chol(P*P')');
X2 = kron((H(ik,:)*X(1:mb,:))',P);
YY = P*Y;
idx = find(strPhi > 0);
sidx = size(idx,1);

V = X2 \ YY(:);

Phi1 = zeros(n*mb,1);
Phi1(idx) = V;
Phi(:,1:mb) = reshape(Phi1,n,mb);
err = [YY(:)-X2*V];
sPhi = reshape(sqrt(diag((err'*err + n)/(n*N)*inv(X2'*X2))),n,mb);

X2 = [H(ik,:)*X(1:mb,:); R5515-D*R2215];

V = Y / X2;

E = V(:,mb+1:mb+m);
E(:,ik) = -V(:,1:mb);
err = [Y-V*X2 iT*iO1*R6666r];
VV = reshape(sqrt(diag(kron(err*err'/N,inv(X2*X2')))),mb+m,n)';
sE = VV(:,mb+1:mb+m);
sE(:,ik) = VV(:,1:mb);

Qr =  Zt*Zt'/N;

if nargout > 4
   if ~ext
           innov = Zt*Q(:,ix(1,1):ix(5,2))';
        else
           innov = Zt*Q((j-1)*s+1:N+(j-1)*s,ix(1,1):ix(5,2))';
        end
     end
   end
%
% end function
