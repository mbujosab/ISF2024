function [tchi1,tchi2,desv1,desv2,angle1,angle2,rvchi1,rvchi2,dg1,dg2] = uroot(y,u,i,s)
% Contraste de raíces unitarias utilizando el estdístico de Tsay and Tiao (19xx)
% para determinar la dimensión del sistema

[tchi1,desv1,angle1,rvchi1,dg1]=sidang(y,u,i,s);
[tchi2,desv2,angle2,rvchi2,dg2]=sidcoin(y,u,i,s);
disp('   n   CanCorr    test  df  p-val   CanCorr    test  df  p-val')
for row=1:i
  fprintf(' %3.0f %8.3f %8.3f %3.0f %6.3f %8.3f %8.3f %3.0f %6.3f\n', ...
     row,angle1(row,row),rvchi1(row,row),dg1(row,row),1-tchi1(row,row),...
         angle2(row,row),rvchi2(row,row),dg2(row,row),1-tchi2(row,row));
end
