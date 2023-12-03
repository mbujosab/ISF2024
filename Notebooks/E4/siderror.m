function ercod = siderror(code, P1, P2, P3, P4, P5)
% SIDERROR  - Devuelve una cadena de error y detiene la ejecución del paquete
%               ercod = e4error(code, P1, P2, P3, P4, P5)
% 23/12/96
% (C@@)

global 	SIDERROR

str1  = ['ERROR ' SIDERROR(code,:)];
if nargin > 1
    evstr = getparms(nargin-1);
    eval(['str1 = sprintf(str1 ' evstr ');']);
end
error(str1);

