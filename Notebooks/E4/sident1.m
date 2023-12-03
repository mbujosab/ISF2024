function [S,Phi,H,E,Qr,A1] = sident1(y, u, i, n, s)
%
% SIDENT -    Algoritmo general de identificación basado en subespacios, tanto para
%             modelos con variables exógenas como modelos puramente estocásticos
%
% Modelos con variables exógenas
%     [Phi,H,E,Q,Gam,D,innov] = sident(y, u, ij, n, s)
%
% Modelos puramente estocásticos
%     [Phi,H,E,Q,innov] = sident(y, [], ij, n, s)
%
% 23/12/96

   global SIDOPTION

   if nargin < 5, s = 1; end

   exact = ~SIDOPTION(1);
   ext   = 0;
   canon = SIDOPTION(3);
   %pond  = SIDOPTION(4);
   oest  = SIDOPTION(5);

   if nargout > 1, t = 1; else, t = 0; end
   
   if size(i,1) == 1, j = i(1,1); else, j = i(2,1); end
   i = i(1,1);

   if n > i | n > j, siderror(5); end

   m = size(y,2);
   r = size(u,2);

   %
     Yi = blkhkel(y, i+j, s, 1);
     N = min(size(Yi,2), size(y,1));

% step 1
     [Q, R] = qr(Yi', 0);
     R = R';
     ix = zeros(3,2);
     ix(:,1) = [1; i*m+1; (i+1)*m+1];
     ix(:,2) = [ix(2:3,1)-1; (i+j)*m];

% step 2
         Om = R(ix(2,1):ix(3,2),ix(2,1):ix(3,2))*R(ix(2,1):ix(3,2),ix(2,1):ix(3,2))';
    
     sOm = sqrtm(Om); % Produce lo mismo que sOm=chol(Om)';
     isOm = pinv(sOm);

     
% step 3 & 4
      [U,S,V] = svd(R(ix(2,1):ix(3,2),ix(2,1):ix(3,2)) \ R(ix(2,1):ix(3,2),ix(1,1):ix(1,2)));    
     %[U S V] = svd(isOm*R(ix(2,1):ix(3,2),ix(1,1):ix(1,2)));

 if t
         
     O = sOm * U(:,1:n);
     iO = U(:,1:n)'*isOm;
     C = S(1:n,1:n) * V(:,1:n)'* R(ix(1,1):ix(1,2),ix(1,1):ix(1,2))'/ N;

     if oest == 0
        [U1 S1 V1] = svd(isOm1*R(ix(3,1):ix(3,2),ix(1,1):ix(2,2)));
        O1 = sOm1 * U1(:,1:n);
        iO1 = pinv(O(1:(j-1)*m,:))*O1*U1(:,1:n)'*isOm1;
     else
        iO1  = pinv(O(1:(j-1)*m,:));
     end

% step 5
     Phi = iO1*O(m+1:j*m,:);
     H   = O(1:m,:);
     G = C(:,m*(i-1)+1:m*i);
     A2 = G;
     L0 = R(ix(2,1):ix(2,2),ix(1,1):ix(2,2))*R(ix(2,1):ix(2,2),ix(1,1):ix(2,2))'/N;
     A3 = L0;

     if exact

% step 6 & 7
        [P, E, Qr] = Ricnit(Phi, H, G, L0);

     else
     %   aprox
     % step 6'

        Qr = (R(ix(2,1):ix(2,2),ix(2,1):ix(2,2))*R(ix(2,1):ix(2,2),ix(2,1):ix(2,2))' + sOm(1:m,:)*U(:,n+1:j*m)*S(n+1:j*m,n+1:j*m)^2*U(:,n+1:j*m)'*sOm(1:m,:)')/N;
        E  = iO1*(R(ix(3,1):ix(3,2),ix(2,1):ix(2,2))*R(ix(2,1):ix(2,2),ix(2,1):ix(2,2))' + sOm(m+1:j*m,:)*U(:,n+1:j*m)*S(n+1:j*m,n+1:j*m)^2*U(:,n+1:j*m)'*sOm(1:m,:)')/N * pinv(Qr);
     %
     end
     if nargout > 5
        if ~ext
           X = S(1:n,1:n)*V(:,1:n)'*Q(:,ix(1,1):ix(1,2))';
           A1 = Yi(ix(2,1):ix(2,2),:) - H*X;
        else
           X = S(1:n,1:n)*V(:,1:n)'*Q((j-1)*s+1:N+(j-1)*s,ix(1,1):ix(1,2))';
           A1 = Yi(ix(2,1):ix(2,2),(j-1)*s+1:N+(j-1)*s) - H*X;
        end
     end    
 end
% end function
