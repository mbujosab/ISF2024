function [opt] = sidopt(o1,v1,o2,v2,o3,v3,o4,v4,o5,v5,o6,v6, ...
                   o7,v7,o8,v8,o9,v9,o10,v10)
%
% SIDOPT -   Modifies the default optiions stored in SIDOPTION.
%   The general syntax is:
%     [opt] = sidopt('option', 'value', ...)
% There can be up to 10 'option','value' pairs
% The possible combinations of options and values are the following:
%
%     option      values
%     ----------  ----------------------------------------
%     Metodo      exacto (d), aproximado
%     extendidas  si (d), no
%     canonica    si, no (d)
%     ponderación no, varianza residual (d), ángulos
%     Obs(i-1)+   reestimada (d), Moore-Penrose Obs(i-1)
%     verbose     si (d), n0
%
% Some special calls are:  
%   SIDOPT         restores default options
%   SIDOPT('show') displays current options
%
% 23/12/96
% (C@@)

global SIDOPTION

if ~exist('SIDOPTION'), siderror(1); end
if ((nargin > 1) & (rem(nargin,2) ~= 0)) | (nargin > 20)
   siderror(2);
end

if nargin == 0 % Valores por defecto
   SIDOPTION(1:6) = [0 1 0 1 0 0];
   opt = SIDOPTION;
   disp('    Fijados valores por defecto');
   nargin = 1; o1 = 'show'; % Fuerza que se enseñen
end
if nargin == 1
   if size(o1,2) < 3, siderror(3); end
   if (o1(1:3) == 'sho')
   % Muestra los valores
    disp('************* Opciones fijadas por el usuario ************');
    if SIDOPTION(1) == 1, s1 = 'Aproximado'; else s1 = 'Exacto'; end
    if SIDOPTION(2) == 1, s2 = 'Sí'; else s2 = 'No'; end
    if SIDOPTION(3) == 1, s3 = 'Sí'; else s3 = 'No'; end
    if SIDOPTION(4) == 0, s4 = 'No'; elseif SIDOPTION(4) == 1, s4 = 'Varianza residual'; else, s4 = 'Angulos'; end
    if SIDOPTION(5) == 0, s5 = 'reestimada'; elseif SIDOPTION(5) == 1, s5 = 'Moore-Penrose O(i-1)'; end    
    if SIDOPTION(6) == 0, s6 = 'Sí'; else s6 = 'No'; end

    disp(['Método ' s1 ]);
    disp(['Matrices Block Hankel extendidas ' s2]);
    disp(['Forma canónica ' s3]);
    disp(['Ponderación ' s4]);
    disp(['Obs(i-1)+ ' s5]);
    disp(['verbose ' s6]);
    disp('**********************************************************');
    disp(' '); disp(' ');
   else
    siderror(4, o1);
   end
   return
end

siddisp('************* Modificadas las siguientes opciones **********');

for i=1:(nargin/2)
   optstr = lower(eval(['o' int2str(i)]));
   optval = eval(['v' int2str(i)]);
   if size(optstr,2) < 3, siderror(3, optstr); end
   if size(optval,2) < 1, siderror(4, optsrt); end
   if optstr(1:3) == 'met'
      optval = lower(optval);
      if optval(1) == 'e', SIDOPTION(1) = 0;
      elseif optval(1) == 'a'
             SIDOPTION(1) = 1;
      else siderror(4, optval); end
      siddisp(['Método ' upper(optval)]);

   elseif optstr(1:3) == 'ext'
      optval = lower(optval);
      if optval(1) == 's', SIDOPTION(2) = 1;
      elseif optval(1) == 'n', SIDOPTION(2) = 0;
      else siderror(4, optval); end
      siddisp(['B-H extendidas ' upper(optval)]);
      
   elseif optstr(1:3) == 'can'
      optval = lower(optval);
      if optval(1) == 's', SIDOPTION(3) = 1;
      elseif optval(1) == 'n', SIDOPTION(3) = 0;
      else siderror(4,optval); end
      siddisp(['Forma canónica ' upper(optval)]);

   elseif optstr(1:3) == 'pon'
      optval = lower(optval);
      if optval(1) == 'n', SIDOPTION(4) = 0;
      elseif optval(1) == 'v', SIDOPTION(4) = 1;
      elseif optval(1) == 'a', SIDOPTION(4) = 2;
      else siderror(4,optval); end
      siddisp(['Ponderación ' upper(optval)]);

   elseif optstr(1:3) == 'obs'
      optval = lower(optval);
      if optval(1) == 'm', SIDOPTION(5) = 1;
      elseif optval(1) == 'r', SIDOPTION(5) = 0;
      else siderror(4,optval); end
      siddisp(['Obs(i-1)+ ' upper(optval)]);

   elseif optstr(1:3) == 'ver'
      optval = lower(optval);
      if optval(1) == 's', SIDOPTION(6) = 0;
      elseif optval(1) == 'n', SIDOPTION(6) = 1;
      else siderror(4,optval); end
      disp(['Verbose ' upper(optval)]);
   
   else
      disp(['Opción no reconocida ' optstr ]);
   end
end
opt = SIDOPTION;
siddisp('************************************************************');
siddisp(' '); siddisp(' ');

