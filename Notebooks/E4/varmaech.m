function [F,sF,Th,sTh,Sigma,A1,A2,A3] = varmaech(y, u, i, ki, s)
%
% VARMAECH -  Algoritmo de identificación basado en subespacios, que permite
%             obtener una forma VARMAX echelon forzando unos índices de observabilidad
%
% Modelos con variables exógenas
%     [F,sF,Th,sTh,Sigma,G,sG,innov] = varmaech(y, u, ij, ki, s)
%
% Modelos puramente estocásticos
%     [F,sF,Th,sTh,Sigma,innov] = varmaech(y, [], ij, ki, s)
%
% Para el modelo:
%    F(B)y(t) = G(B)u(t) + Th(B)e(t)
% identifica las matrices:
% F = [F0 F1 ... Fk], Th = [F0 Th1 ... Thk], G = [G0 G1 ... Gk], V = V(e(t))
% donde k es el máximo de los índices de observabilidad.
% Estas matrices incorporan las restricciones dadas por la estructura
% de índices de Kronecker. Se devuelven además las desviaciones
% de los parámetros estimados en sF, sTh y sG.
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
     isOm = inv(sOm);

     if oest == 0
        if pond == 0
           Om1 = eye((j-1)*m);
        elseif pond == 1
           Om1 = R(ix(7,1):ix(7,2),ix(7,1):ix(7,2))*R(ix(7,1):ix(7,2),ix(7,1):ix(7,2))'/N;
        else
           Om1 = R(ix(7,1):ix(7,2),ix(2,1):ix(7,2))*R(ix(7,1):ix(7,2),ix(2,1):ix(7,2))';
        end
        sOm1 = sqrtm(Om1);
        isOm1 = inv(sOm1);
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
     [Phi, H, T] = echelon(ki, Phi, H);
     iT = inv(T);
     [strH, strPhi] = echelstr(ki);
     strH = strH(:);
     strPhi = strPhi(:);
     P = inv(R(ix(6,1):ix(6,2),ix(6,1):ix(6,2)));

     X = [iT*iO*R(ix(6,1):ix(7,2),ix(1,1):ix(5,2)) zeros(n,m)];
     Y = R(ix(6,1):ix(6,2),ix(1,1):ix(6,2));

     ik = find(ki > 0);
     YY = Y;
     YY(ik,:) = Y(ik,:) - X(1:mb,:);
     YY = P*YY;

     X2 = kron(X(1:mb,:)',P);
     idx = find(strH > 0);
     sidx = size(idx,1);
     X2 = [X2(:,idx) kron(R(ix(1,1):ix(2,2),ix(1,1):ix(6,2))',P)];

     V = X2 \ YY(:);

     H1 = zeros(m*mb,1);
     H1(idx) = V(1:sidx);
     H1 = reshape(H1,m,mb);
     Di = reshape(V(sidx+1:size(V,1)),m,r*j);
     H = eye(m);
     H = H1 + H(:,ik);
     err = [YY(:)-X2*V];
     VV  = sqrt(diag((err'*err + m)/(m*N)*inv(X2'*X2)));
     H1(idx) = VV(1:sidx);
     sH = reshape(H1,m,mb);

     Zt = Y-H*X(1:mb,:)-Di*R(ix(1,1):ix(2,2),ix(1,1):ix(6,2));

     if n > mb
        Y = iT*iO1*R(ix(7,1):ix(7,2),ix(1,1):ix(6,2)) - Phi(:,mb+1:n)*X(mb+1:n,:);
     else
        Y = iT*iO1*R(ix(7,1):ix(7,2),ix(1,1):ix(6,2));
     end

     X2 = [X(1:mb,:); R(ix(1,1):ix(2,2),ix(1,1):ix(6,2))];
     V = Y / X2;

     P = [Y-V*X2 iT*iO1*R(ix(7,1):ix(7,2),ix(7,1):ix(7,2))];
     P = inv(chol(P*P')');
     X2 = kron((H(ik,:)*X(1:mb,:))',P);

     YY = P*Y;
     idx = find(strPhi > 0);
     sidx = size(idx,1);

     X2 = [X2(:,idx) kron(R(ix(1,1):ix(2,2),ix(1,1):ix(6,2))',P)];

     V = X2 \ YY(:);

     Phi1 = zeros(n*mb,1);
     Phi1(idx) = V(1:sidx);
     Phi(:,1:mb) = reshape(Phi1,n,mb);
     Gami = reshape(V(sidx+1:sidx+n*r*j), n, r*j);
     VV = inv(N*X2'*X2);
     sPhi = zeros(n*mb,1);
     sPhi(idx) = sqrt(diag(VV(1:sidx,1:sidx)));
     sPhi = reshape(sPhi,n,mb);

     X2 = [H(ik,:)*X(1:mb,:);Zt+H*X(1:mb,:);R(ix(1,1):ix(2,2),ix(1,1):ix(6,2))];

     V = Y / X2;
     E = V(:,mb+1:mb+m);
     E(:,ik) = -V(:,1:mb);
     err = [Y-V*X2 iT*iO1*R(ix(7,1):ix(7,2),ix(7,1):ix(7,2))];
     VV = reshape(sqrt(diag(kron(err*err'/N,inv(X2*X2')))),mb+m+r*j,n)';
     sE = VV(:,mb+1:mb+m);
     sE(:,ik) = VV(:,1:mb);

     Sigma =  Zt*Zt'/N;

     [D Gam] = sidverg([[Di(:,(j-1)*r+1:j*r) Di(:,1:(j-1)*r)]; [Gami(:,(j-1)*r+1:j*r) Gami(:,1:(j-1)*r)]], [eye(m) zeros(m,(j-1)*m); zeros(n,m) iT*iO1]*[eye(j*m)-O*iO], O*T, j, m, r);
     F0 = eye(m);
     F0(:,ik) = H;
     F0 = inv(F0);
     sF0 = zeros(m);
     sF0(:,ik) = sH;
     G0 = F0*D;
     G1 = Gam - Phi(:,1:mb)*D(ik,:);

     X = [zeros(n,j*r) iT*iO*R(ix(6,1):ix(7,2),ix(3,1):ix(5,2)) zeros(n,m)];
     if n > mb
        Y = iT*iO1*R(ix(7,1):ix(7,2),ix(1,1):ix(6,2)) - Phi(:,mb+1:n)*X(mb+1:n,:);
     end

     X2 = [H(ik,:)*X(1:mb,:)+D(ik,:)*R(ix(2,1):ix(2,2),ix(1,1):ix(6,2));R(ix(2,1):ix(2,2),ix(1,1):ix(6,2))];
     V = Y / X2;
     err = [Y-V*X2 iT*iO1*R(ix(7,1):ix(7,2),ix(7,1):ix(7,2))];
     VV = reshape(sqrt(diag(kron(err*err'/N,inv(X2*X2')))),mb+r,n)';
     sG1 = VV(:,mb+1:mb+r);

     X2 = [X;R(ix(2,1):ix(2,2),ix(1,1):ix(6,2))];
     V = F0*R(ix(6,1):ix(6,2),ix(1,1):ix(6,2)) / X2;
     err = F0*R(ix(6,1):ix(6,2),ix(1,1):ix(6,2))-V*X2;
     VV = reshape(sqrt(diag(kron(err*err'/N,inv(X2*X2')))),n+r,m)';
     sG0 = VV(:,n+1:n+r);

     if nargout > 6
        if ~ext
           A3 = Zt*Q(:,ix(1,1):ix(6,2))';
        else
           A3 = Zt*Q((j-1)*s+1:N+(j-1)*s,ix(1,1):ix(6,2))';
        end
     end

     [F,Th,A1] = sstoech(ki,F0,-Phi(:,1:mb),F0,E,G0,G1);
     [sF,sTh,A2] = sstoech(ki,sF0,sPhi,zeros(m),sE,sG0,sG1);
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
     isOm = inv(sOm);

     if oest == 0
        if pond == 0
           Om1 = eye((j-1)*m);
        elseif pond == 1
           Om1 = R(ix(3,1):ix(3,2),ix(3,1):ix(3,2))*R(ix(3,1):ix(3,2),ix(3,1):ix(3,2))'/N;
        else
           Om1 = R(ix(3,1):ix(3,2),ix(1,1):ix(3,2))*R(ix(3,1):ix(3,2),ix(1,1):ix(3,2))';
        end
        sOm1 = sqrtm(Om1);
        isOm1 = inv(sOm1);
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

     [Phi, H, T] = echelon(ki, Phi, H);
     iT = inv(T);

     [strH, strPhi] = echelstr(ki);
     strH = strH(:);
     strPhi = strPhi(:);
     P = inv(R(ix(2,1):ix(2,2),ix(2,1):ix(2,2)));

     X = [iT*S(1:n,1:n)*V(:,1:n)' zeros(n,m)];
     Y = R(ix(2,1):ix(2,2),ix(1,1):ix(2,2));

     ik  = find(ki > 0);
     YY = Y;
     YY(ik,:) = Y(ik,:) - X(1:mb,:);
     YY = P*YY;

     X2 = kron(X(1:mb,:)',P);
     idx = find(strH > 0);
     X2 = X2(:,idx);

     V = X2 \ YY(:);

     H1 = zeros(m*mb,1);
     H1(idx) = V;
     H1 = reshape(H1,m,mb);
     H = eye(m);
     H = H1 + H(:,ik);
     err = [YY(:)-X2*V];
     VV  = sqrt(diag((err'*err + m)/(m*N)*inv(X2'*X2)));
     H1(idx) = VV;
     sH = reshape(H1,m,mb);
     Zt = Y-H*X(1:mb,:);

     if n > mb
        Y = iT*iO1*R(ix(3,1):ix(3,2),ix(1,1):ix(2,2)) - Phi(:,mb+1:n)*X(mb+1:n,:);
     else
        Y = iT*iO1*R(ix(3,1):ix(3,2),ix(1,1):ix(2,2));
     end

     P = [Y-Phi(:,1:mb)*X(1:mb,:) iT*iO1*R(ix(3,1):ix(3,2),ix(3,1):ix(3,2))];
     P = inv(chol(P*P')');
     X2 = kron((H(ik,:)*X(1:mb,:))',P);

     YY = P*Y;
     idx = find(strPhi > 0);
     sidx = size(idx,1);

     V = X2(:,idx) \ YY(:);

     Phi1 = zeros(n*mb,1);
     Phi1(idx) = V(1:sidx);
     Phi(:,1:mb) = reshape(Phi1,n,mb);
     sPhi = zeros(n*mb,1);
     sPhi(idx) = sqrt(diag(inv(N*X2(:,idx)'*X2(:,idx))));
     sPhi = reshape(sPhi,n,mb);

     X2 = [H(ik,:)*X(1:mb,:);R(ix(2,1):ix(2,2),ix(1,1):ix(2,2))];

     V = Y / X2;
     E = V(:,mb+1:mb+m);
     E(:,ik) = -V(:,1:mb);
     err = [Y-V*X2 iT*iO1*R(ix(3,1):ix(3,2),ix(3,1):ix(3,2))];
     VV = reshape(sqrt(diag(kron(err*err'/N,inv(X2*X2')))),mb+m,n)';
     sE = VV(:,mb+1:mb+m);
     sE(:,ik) = VV(:,1:mb);

     Sigma =  Zt*Zt'/N;

     if nargout > 4
        if ~ext
           A1 = Zt*Q(:,ix(1,1):ix(2,2))';
        else
           A1 = Zt*Q((j-1)*s+1:N+(j-1)*s,ix(1,1):ix(2,2))';
        end
     end

     F0 = eye(m);
     F0(:,ik) = H;
     F0 = inv(F0);
     sF0 = zeros(m);
     sF0(:,ik) = sH;

     [F,Th] = sstoech(ki,F0,-Phi(:,1:mb),F0,E);
     [sF,sTh] = sstoech(ki,sF0,sPhi,zeros(m),sE);

   end
%
% end function




