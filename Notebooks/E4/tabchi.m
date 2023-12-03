function [tchi2,desv,angle,rvchi,dg] = tabchi(y,u,i,s)
%
% TABCHI - Computes an display statistics to determine
%          the dimension of the system.
%

[tchi2,desv,angle,rvchi,dg] = sidang(y,u,i,s);
disp('  n  ccorr     test df  p-val')
for row=1:i
  fprintf(' %2.0f %6.3f %8.3f %2.0f %6.3f\n', ...
     row,angle(row,row),rvchi(row,row),dg(row,row),1-tchi2(row,row));
end
disp('')
