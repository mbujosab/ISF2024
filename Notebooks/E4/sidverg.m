function [D, Gam] = sidverg(OHd, Op, O, i, m, r)
%
% SIDVERG -   Resuelve el sistema de ecuaciones propuesto por Viberg (1994)
%            para la identificación de las matrices D y Gam
%
%    [D, Gam] = sidverg(OHd, Op, O, i, m, r)
%
% 2/1/97


   n = size(O,2);
   l = size(OHd,1);

   YY = zeros(l*i,r);
   XX = zeros(l*i,n+m);

   for k=1:i-1
       YY((k-1)*l+1:l*k,:) = OHd(:,(k-1)*r+1:r*k);
       XX((k-1)*l+1:l*k,:) = [Op(:,(k-1)*m+1:m*k) Op(:,m*k+1:m*i)*O(1:m*(i-k),:)];
   end

   YY((i-1)*l+1:l*i,:) = OHd(:,(i-1)*r+1:r*i);
   XX((i-1)*l+1:l*i,:) = [Op(:,(i-1)*m+1:m*i) zeros(l,n)];

   DGam = XX \ YY;
   D = DGam(1:m,:);
   Gam = DGam(m+1:m+n,:);