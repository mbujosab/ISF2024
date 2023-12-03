function S = corrcan(y, u, i);
%
%
%

    ext = 0;
    s = 1;
    m = size(y,2);

    if isempty(u)
       r  = 0;
       Yi = blkhkel(y, 2*i, s, ext);
       Yj = Yi(i*m+1:(2*i)*m,:);
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

    [U S V] = svd(Qi'*Qj);
    S = diag(S);