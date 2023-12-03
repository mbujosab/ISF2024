function [tchi2, desv, angle, rvchi, dg] = sicoin(y, u, i, s)
%
%   Contraste de Tsay y Tiao para determinar el orden del sistema
%
%   [tchi2, desv, angle, rvchi, dg] = sidang(y, u, i, s)
%

    global SIDOPTION
    if nargin < 4, s = 1; end

%    ext   = SIDOPTION(2);
    ext   = 0;   % No utilizar este test con matrices extendidas

    dy = transdif(y,1,1);

    m = size(y,2);

    if isempty(u)
       r  = 0;
       Yi = blkhkel(y, 2*i+1, s, ext);
       Yj = Yi((i+1)*m+1:(2*i+1)*m,:) - Yi(i*m+1:(2*i)*m,:);
       Yi = Yi(1:i*m,:);
    else
       r = size(u,2);
       Yi = blkhkel(y, 2*i, s, ext);
       Ui = blkhkel(u, 2*i, s, ext);
       [Q, R] = qr([Ui(r*i+1:2*r*i,:);Ui(1:r*i,:);Yi]', 0);
       R = R';
       Yi = R(i*r+1:i*(2*r+m),i*r+1:i*(2*r+m))*Q(:,i*r+1:i*(2*r+m))';
       Yj = R(i*(2*r+m)+1:2*i*(m+r),i*r+1:2*i*(m+r))*Q(:,i*r+1:2*i*(m+r))';
    end


    [Qi, R] = qr(Yi', 0);
    [Qj, R] = qr(Yj', 0);

    N = size(Qi,1);

    desv  = zeros(m*i);
    tchi2 = zeros(m*i);
    dg    = zeros(m*i);
    angle = zeros(m*i);
    rvchi = zeros(m*i);

    for j = 1:i
    for k = 1:m
    %
        s = (j-1)*m+k;
        [U S V] = svd(Qi'*Qj(:,1:s));
        ang = diag(S(j:s,j:s));
        U = Qi * U(:,j:s);
        V = Qj(:,1:s) * V(:,j:s);

        sz = s-j+1;
        ds = zeros(sz,1);

        for h = 1:j-1
        %
            ds = ds + diag((U(1:N-h,:)'*U(1+h:N,:)).*(V(1:N-h,:)'*V(1+h:N,:)));
        %
        end

        dv = 1 + 2*ds;
        desv(s,j:j+sz-1) = dv';
        vchi2 =  -N*flipud(cumsum(log(1-(ang(sz:-1:1).^2)./dv)));
        tchi2(s,j:j+sz-1) = chi2cdf(vchi2, (s-[1:sz]'-j+2).*(i*(m+r)-[1:sz]'-j+2))';
        rvchi(s,j:j+sz-1) = vchi2';
        dg(s,j:j+sz-1) = ((s-[1:sz]'-j+2).*(i*(m+r)-[1:sz]'-j+2))';
        angle(s,j:j+sz-1) = ang(1:sz)';
    %
    end
    end

    dv = sqrt(diag(desv))/sqrt(N);

