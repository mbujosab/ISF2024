function [F,Th,G] = sstoech(ki,F0,F1,Th0,Th1,G0,G1)
%
   k  = max(ki);
   ik = find(ki > 0);
   mb = sum(ki > 0);
   m  = size(F0,1);
   r  = 0;

   F  = [F0 zeros(m,k*m)];
   Th = [Th0 zeros(m,k*m)];
   if nargin > 5
      r  = size(G0,2);
      G  = [G0 zeros(m,k*r)];      
   end

   n = 0;
   for i=1:k
       ik2 = find(ki >= i);
       mb2 = sum(ki >= i);
       F(ik2,m*i+ik)         = F1(n+1:n+mb2,:);
       Th(ik2,m*i+1:m*(i+1)) = Th1(n+1:n+mb2,:);
       if r, G(ik2,r*i+1:r*(i+1))  = G1(n+1:n+mb2,:); end
       n = n+mb2;
   end