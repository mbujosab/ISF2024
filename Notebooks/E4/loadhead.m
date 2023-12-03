function loadhead
% loadhead - IRÁ FUERA

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


[file, direc] = uigetfile('*.TXT', 'Archivo de Ayudas')
eval(['cd ' direc ]);

fid = fopen(file,'r');
F = fread(fid);
fclose(fid);

s = char(F');

if s(1,size(s,2))~=10, s = [s 13 10]; end

kl = find(s == 10);
k  = find(s == '@');
k2 = find(s(k+1) == '$');   % identificamos fichero por @$
k  = k(k2)+2;

for h=1:size(k2,2)
%
    smin = min(k(h)+32,size(s,2));
    pf = find(s(k(h):smin) == 32 | s(k(h):smin) == 10 | s(k(h):smin) == 13);
    file = s(k(h):k(h)+pf-2);

    kp = find(kl > k(h));

    i = kp(1)+1;
    j = kl(kp(1))+1;

    datos = [];

    while j <= size(s,2)
    %
       cad1 = s(j:kl(i));
       k2 = find(cad1 ~= 32 & cad1 ~= 13 & cad1 ~= 10);
       if isempty(k2),  break; end
       if sum(sum(k2(1)+j+1 == k)), break; end

       %datos = [datos '% ' s(j:kl(i))];
       datos = [datos s(j:kl(i))];

       j = kl(i)+1;
       i = i + 1;
    %
    end

    file, datos
    chanhead(file, datos);
end





