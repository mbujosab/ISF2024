function [f, zh, QT] = fvbekk(theta, z, nostat)
% function [f, hominnov, QT] = fvbekk(theta, z, nostat)
%
% theta    > parameter vector.
%            theta(1:m*(m+1)/2): cholesky factor of Lambda (lower triangle)
%            theta(m*(m+1)/2+1:m*(m+1)/2+m*m): A
%            theta(m*(m+1)/2+m*m+1:m*(m+1)/2+2*m*m): B
%            where
%            Q(t+1) = Lambda + A*z(t,:)'*Z(t,:)*A' + B*Q(t)*B'
% z        > matrix of observable variables.
% nostat   > 1 for non stationary models (Q(1) = Lambda)
%            0 for stationary models (Q(1) = Lambda + A*Q(1)*A' + B*Q(1)*B')
% f        < value of the likelihood function.
% hominnov < (optional) stores the sequence of standarized innovations.
% QT       < (optional) stores the sequence of conditional variances.
%

global EEOPTION

if nargin < 2,  e4error(3); end
if nargin < 3, nostat = 0; end

saveinn = 0; if nargout > 1, saveinn = 1; end
scaleb = EEOPTION(2);
zeps   = EEOPTION(15);

m = round((sqrt(40*size(theta,1) + 1) - 1)/10);
if m ~= size(z,2), e4error(11); end
m2 = m*(m+1)/2;
R = vech2m(theta(1:m2,1),m,1);
A = reshape(theta(m2+1:m2+m*m,1),m,m);
B = reshape(theta(m*m+m2+1:m2+2*m*m,1),m,m);
L = R*R';

if nostat
   Q = L;
else
   Q = reshape(pinv(eye(m*m) - kron(A,A) - kron(B,B))*L(:),m,m);
   if sum(eig(Q) < 0), disp('Modelo no estacionario'); end
end

n  = size(z,1);

if saveinn, QT = zeros(n*m,m); end
zh = zeros(size(z));
ff = 0;

for t = 1:n

    if scaleb, U = cholp(Q, abs(Q));
    else       U = cholp(Q); end
    if saveinn, QT((t-1)*m+1:t*m,:) = Q; end
    zh(t,:) = z(t,:)/U; 
    Q = L + A*z(t,:)'*z(t,:)*A' + B*Q*B';
    ff  = ff + 2*sum(log(diag(U))) + zh(t,:)*zh(t,:)';

end  % t

f = 0.5*(ff + n*m*log(2*pi));
