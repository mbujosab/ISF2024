function chanhead(file, data)
% CHANHEAD - IRÁ FUERA

% Copyright (C) 1997 Jaime Terceiro
% 
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2, or (at your option)
% any later version.
% 
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details. 
% 
% You should have received a copy of the GNU General Public License
% along with this file.  If not, write to the Free Software Foundation,
% 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.

fid = fopen(file,'r');
F = fread(fid);
fclose(fid);

s = char(F');

k = find(s == 10);

sout = [];
i = 1;
j = 1;
while 1   % quitamos las primeras líneas en blanco
%
   cad1 = s(j:k(i));
   k2 = find(cad1 ~= 32 & cad1 ~= 13 & cad1 ~= 10);
   cad2 = cad1(k2);
   if cad2(1)=='%', break; end
   sout = [sout cad1];
%   if ~isempty(k2),  break; end
   j = k(i)+1;
   i = i + 1;
%
end
   
cad2 = cad1(k2);

insert = 1;

if cad2(1) ~= '%'   % No aparece primero el nombre de la función
   if size(cad2,2) >= 8
      if cad1(1:8) == 'function'
         sout = cad1;
         j = k(i)+1;
         i = i + 1;
      else
         sout = s;
         insert = 0;
      end
   end
end

if insert

   while i <= size(k,2)
   %
      cad1 = s(j:k(i));
      k2 = find(cad1 ~= 32 & cad1 ~= 13 & cad1 ~= 10);
      if isempty(k2), break; end
      if cad1(k2(1)) ~= '%', break; end
      j = k(i)+1;
      i = i + 1;
   %
   end

   %sout = [sout data 13 10 s(j:size(s,2))];
   sout = [sout data s(j:size(s,2))];

end 

fid = fopen(file,'w');
fwrite(fid, sout');
fclose(fid);