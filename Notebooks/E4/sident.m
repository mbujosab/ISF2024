function [Phi,H,E,Qr,A1,A2,A3,A4,A5] = sident(y, u, i, n, s)
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
   ext   = SIDOPTION(2);
   canon = SIDOPTION(3);
   pond  = SIDOPTION(4);
   oest  = SIDOPTION(5);

   if size(i,1) == 1, j = i(1,1); else, j = i(2,1); end
   i = i(1,1);

   if n > i | n > j, siderror(5); end

   m = size(y,2);
   r = size(u,2);

   if r
   %
     Yi = blkhkel(y, i+j, s, ext);
     Ui = blkhkel(u, i+j, s, ext);
     N = min(size(Yi,2), size(y,1));

% step 1
     [Q, R] = qr([Ui(r*(i+1)+1:r*(i+j),:);Ui(r*i+1:r*(i+1),:);Ui(1:r*i,:);Yi]', 0);
     R = R';
     ix = zeros(7,2);
     ix(:,1) = [1; (j-1)*r+1; j*r+1; (j+i)*r+1; (j+i)*r+m+1; (j+i)*r+i*m+1; (j+i)*r+(i+1)*m+1];
     ix(:,2) = [ix(2:7,1)-1; (i+j)*(m+r)];

% step 2
     if pond == 0
         Om = eye(j*m);
     elseif pond == 1
         Om = R(ix(6,1):ix(7,2),ix(6,1):ix(7,2))*R(ix(6,1):ix(7,2),ix(6,1):ix(7,2))'/N;
     else
         Om = R(ix(6,1):ix(7,2),ix(3,1):ix(7,2))*R(ix(6,1):ix(7,2),ix(3,1):ix(7,2))';
     end

     sOm = sqrtm(Om);
     isOm = pinv(sOm);

     if oest == 0
        if pond == 0
           Om1 = eye((j-1)*m);
        elseif pond == 1
           Om1 = R(ix(7,1):ix(7,2),ix(7,1):ix(7,2))*R(ix(7,1):ix(7,2),ix(7,1):ix(7,2))'/N;
        else
           Om1 = R(ix(7,1):ix(7,2),ix(2,1):ix(7,2))*R(ix(7,1):ix(7,2),ix(2,1):ix(7,2))';
        end
        sOm1 = sqrtm(Om1);
        isOm1 = pinv(sOm1);
     end

% step 3
     [U S V] = svd(isOm*R(ix(6,1):ix(7,2),ix(3,1):ix(5,2)));

     O = sOm * U(:,1:n);
     iO = U(:,1:n)'*isOm;

     if oest == 0
        [U1 S1 V1] = svd(isOm1*R(ix(7,1):ix(7,2),ix(2,1):ix(6,2)));
        O1 = sOm1 * U1(:,1:n);
        iO1 = pinv(O(1:(j-1)*m,:))*O1*U1(:,1:n)'*isOm1;
     else
        iO1  = pinv(O(1:(j-1)*m,:));
     end

% step 4
     Phi = iO1*O(m+1:j*m,:);
     H   = O(1:m,:);

% step 5
     OHd  = U(:,n+1:j*m)'*[R(ix(6,1):ix(7,2),ix(1,1):ix(2,2)) R(ix(4,1):ix(5,2),ix(1,1):ix(3,2))]*...
             pinv([[R(ix(2,1):ix(2,2),ix(1,1):ix(2,2)); R(ix(1,1):ix(1,2),ix(1,1):ix(2,2))] R(ix(3,1):ix(3,2),ix(1,1):ix(3,2))]);

% step 6
     [D Gam] = sidverg(OHd, U(:,n+1:j*m)', O, j, m, r);


     if exact

   % step 7

        P0d = lyapunov(Phi, iO*(R(ix(4,1):ix(5,2),ix(4,1):ix(5,2))*R(ix(4,1):ix(5,2),ix(4,1):ix(5,2))' - R(ix(5,1):ix(6,2),ix(4,1):ix(6,2))*R(ix(5,1):ix(6,2),ix(4,1):ix(6,2))')*iO'/N);

%        if find(eig(P0d) < 0), disp('P0d no presenta el signo correcto, probar con método aproximado'); end
        if isempty(P0d) | ~isempty(find(eig(P0d) < 0)), P0d = zeros(n); end

        L = R(ix(4,1):ix(5,2),ix(4,1):ix(5,2))*R(ix(4,1):ix(5,2),ix(4,1):ix(5,2))'/N - O*P0d*O';
        C = iO*R(ix(6,1):ix(7,2),ix(4,1):ix(5,2))*R(ix(4,1):ix(5,2),ix(4,1):ix(5,2))'/N - Phi^i*P0d*O';

        G = C(:,m*(i-1)+1:m*i);
        A4 = G;
        L0 = zeros(m);
        for k=1:i
            L0 = L0 + L(m*(k-1)+1:k*m,m*(k-1)+1:k*m);
        end
        L0 = L0/i;
        A5 = L0;

   % step 8 & 9
        [P, E, Qr] = Ricnit(Phi, H, G, L0);

     else
     %   aprox
     % step 7'

        Qr = (R(ix(6,1):ix(6,2),ix(6,1):ix(6,2))*R(ix(6,1):ix(6,2),ix(6,1):ix(6,2))' + sOm(1:m,:)*U(:,n+1:j*m)*S(n+1:j*m,n+1:j*m)^2*U(:,n+1:j*m)'*sOm(1:m,:)')/N;
        E  = iO1*(R(ix(7,1):ix(7,2),ix(6,1):ix(6,2))*R(ix(6,1):ix(6,2),ix(6,1):ix(6,2))' + sOm(m+1:j*m,:)*U(:,n+1:j*m)*S(n+1:j*m,n+1:j*m)^2*U(:,n+1:j*m)'*sOm(1:m,:)')/N * pinv(Qr);
     %
     end

     A1 = Gam;
     A2 = D;

     if nargout > 6
        if ~ext
           A3 = R(ix(6,1):ix(6,2),ix(6,1):ix(6,2))*Q(:,ix(6,1):ix(6,2))' + sOm(1:m,:)*U(:,n+1:j*m)*S(n+1:j*m,n+1:(i+r)*m)*V(:,n+1:(i+r)*m)'*Q(:,ix(3,1):ix(5,2))';
        else
           A3 = R(ix(6,1):ix(6,2),ix(6,1):ix(6,2))*Q((j-1)*s+1:N+(j-1)*s,ix(6,1):ix(6,2))' + sOm(1:m,:)*U(:,n+1:j*m)*S(n+1:j*m,n+1:(i+r)*m)*V(:,n+1:(i+r)*m)'*Q((j-1)*s+1:N+(j-1)*s,ix(3,1):ix(5,2))';
        end
     end
%
   else % r = 0, puramente estocásticos
   %
     Yi = blkhkel(y, i+j, s, ext);
     N = min(size(Yi,2), size(y,1));

% step 1
     [Q, R] = qr(Yi', 0);
     R = R';
     ix = zeros(3,2);
     ix(:,1) = [1; i*m+1; (i+1)*m+1];
     ix(:,2) = [ix(2:3,1)-1; (i+j)*m];

% step 2
     if pond == 0
         Om = eye(j*m);
     elseif pond == 1
         Om = R(ix(2,1):ix(3,2),ix(2,1):ix(3,2))*R(ix(2,1):ix(3,2),ix(2,1):ix(3,2))'/N;
     else
         Om = R(ix(2,1):ix(3,2),ix(1,1):ix(3,2))*R(ix(2,1):ix(3,2),ix(1,1):ix(3,2))';
     end

     sOm = sqrtm(Om);
     isOm = pinv(sOm);

     if oest == 0
        if pond == 0
           Om1 = eye((j-1)*m);
        elseif pond == 1
           Om1 = R(ix(3,1):ix(3,2),ix(3,1):ix(3,2))*R(ix(3,1):ix(3,2),ix(3,1):ix(3,2))'/N;
        else
           Om1 = R(ix(3,1):ix(3,2),ix(1,1):ix(3,2))*R(ix(3,1):ix(3,2),ix(1,1):ix(3,2))';
        end
        sOm1 = sqrtm(Om1);
        isOm1 = pinv(sOm1);
     end

% step 3 & 4
     [U S V] = svd(isOm*R(ix(2,1):ix(3,2),ix(1,1):ix(1,2)));

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

     if nargout > 4
        if ~ext
           X = S(1:n,1:n)*V(:,1:n)'*Q(:,ix(1,1):ix(1,2))';
           A1 = Yi(ix(2,1):ix(2,2),:) - H*X;
        else
           X = S(1:n,1:n)*V(:,1:n)'*Q((j-1)*s+1:N+(j-1)*s,ix(1,1):ix(1,2))';
           A1 = Yi(ix(2,1):ix(2,2),(j-1)*s+1:N+(j-1)*s) - H*X;
        end
     end
   end
%
% end function
