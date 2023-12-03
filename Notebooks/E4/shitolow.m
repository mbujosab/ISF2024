function [J] = shitolow(th,dh, idxagg, typeagg, s)

[H_D, type, m, r, s0, n, np, userflag, userf, innov, szpriv] = e4gthead(dh);
if nargin < 5
   s = s0;
   if nargin < 4
      typeagg = 0;
      if nargin < 3
         idxagg = ones(m,1);           
      end
   end   
end

dif = 1e-8;
[tl,dl,ll] = hitolow(th,dh, idxagg, typeagg, s);

if size(th,2) < 2
   sth = size(th,1);    
   k = (1:sth)';
else
   k = find(~th(:,2));
   sth = size(k,1);
end   

J = zeros(size(tl,1),sth);

for i = 1:sth
    th2 = th;
    th2(k(i),1) = th2(k(i),1) + dif;
    J(:,i) = (hitolow(th2,dh, idxagg, typeagg, s) - tl)/dif;
end