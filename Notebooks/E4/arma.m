function [Phi, H, T, E, G] = arma(A, C, E, G, D)
%
%   Modelos con variables exógenas
%        [Phic, Hc, T, Ec, Gc] = arma(Phi, H, E, Gam, D)
%   Modelos puramente estocásticos
%        [Phic, Hc, T, Ec] = arma(Phi, H, E)
%
	n = size(A,1);
        m = size(C,1);

        if rem(n,m) > 0
           error('n debe ser múltiplo entero de m');
        end

        k = n / m;

        Phi = zeros(n);
        H   = zeros(m,n);
        H(:, 1:m) = eye(m);

        T = zeros(n);
	O = zeros(n);

        for i=1:k
	%
            O(n-i*m+1:n-(i-1)*m,:) = C * A^(i-1);
        %
        end

	rk = rank(O);

        if rk < m
           error('rango de la matriz de observabilidad inferior a la dimensión del vector de observación');
        end

	if rk < n
	%

	O
	   disp('Matriz de observabilidad deficiente en rango');
	   T(1:rk, n-m+1:n) = inv(O(1:rk,1:rk)) * H(:,1:rk)';

	T

	O * T(:, n-m+1:n)
	%
	else
	%
           T(:, n-m+1:n) = inv(O) * H';
	%
	end

        for i=2:k
       	%
            T(:, n-i*m+1:n-(i-1)*m) = A * T(:, n-(i-1)*m+1:n-(i-2)*m);
        %
        end

	if (rk < n)
	   Phi(1:rk,:) = inv(T(1:rk,1:rk)) * A(1:rk,:) * T;
	else
	   Phi = inv(T) * A * T;
	end

	E = inv(T) * E - Phi(:,1:m);

        if nargin > 3, G = inv(T) * G - Phi(:,1:m)*D; end

	H = C * T;

%
        