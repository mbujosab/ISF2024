%
% SIDINIT  - Initializes the default options stored in the global vector
%            SIDOPTION
%
% This function should be executed at least once at the beginning of a session
%
% 23/12/96
% (C@@)

global SIDOPTION SIDERROR SIDWARN SIDDISP
if isempty(SIDOPTION)
   SIDOPTION = zeros(1,10);
end
sidopt;

% Inicialización de los mensajes de error
SIDERROR = ['1. Ejecute SIDINIT antes de usar la librería                         ';
            '2. Número de argumentos incorrecto                                   ';
            '3. SIDOPT: Opción no reconocida: %s                                  ';
            '4. SIDOPT: Valor no reconocido: %s                                   ';
            '5. n debe ser menor o igual a i y a j                                '];
            
% Inicialización de los mensajes de aviso
SIDWARN = [];

% Inicialización de mensajes de información
SIDDISP = [];

