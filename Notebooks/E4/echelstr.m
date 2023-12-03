function [strH, strF] = echelstr(ki)
%
  m = size(ki,1);
  n = sum(ki);

  ik = find(ki > 0);
  mb = size(ik,1);
  ki2 = ki(ik);

  strH = tril(ones(m),-1);
  strH = strH(:,ik);
  strF = zeros(n,mb);

  for k=1:m
      idx = find(ki(ik(ik < k)) <= ki(k));
      strH(k,idx) = zeros(1,size(idx,1));
  end

  din = 1;
  k = 1;
  while k <= n
        ik2 = find(ki2 >= din);
        for h=1:size(ik2,1)
            idx = find(ki2 + din > ki2(ik2(h)));
            strF(k,idx) = ones(1, size(idx,1));
            k = k + 1;
        end
        din = din + 1;
  end

