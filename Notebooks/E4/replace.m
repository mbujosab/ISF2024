function N = replace(filein, fileout, c, c2, ci)
%
% ci : case insensitive

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

if nargin < 5, ci = 0; end

fid = fopen(filein,'r');
F = fread(fid);
fclose(fid);

s = char(F');
l = size(c,2);

if ci, k = findstr(upper(s),upper(c)); else, k = findstr(s,c); end
N = size(k,2);

if ~isempty(k)
   s2 = [];
   ptr = 1;
   for i=1:N
       s2 = [s2 s(ptr:k(i)-1) c2];
       ptr = k(i) + l;
   end
   s2 = [s2 s(ptr:size(s,2))];
else
   s2 = s;
end

fid = fopen(fileout,'w');
fwrite(fid, s2');
fclose(fid);
