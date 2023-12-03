function [irf] = plotirf(Phi, Gam, E, H, D, C, Q, S, R, hor)
% PLOTIRF - Displays standarized graphs of an impulse response function of a serie ****
%	irf = plotirf(Phi, Gam, E, H, D, C, Q, S, R, hor)
% hor	> horizon
% irf	< matrix of standarized impulse response function
% **argumentos
% **argumentos
% **argumentos
% **argumentos
% **argumentos
% Copyright (c) Jaime Terceiro, 1998


% Covariance matrix

V = [Q S; S' R];

if Q == S & Q == R,
V = Q;
end

% Cholesky Decomposition

K = chol(V);
P = K';

% Eigenvalues and eigenvectors

[W,L]=eig(V);

% Z matrix

Z = inv(W)*P;

% Impulse Response Function

[m_states,l_exog] = size(Gam);
[m_states,k_noises] = size(E);
[n_endog,l_exog] = size(D);
[n_endog,m_states]  = size(H);

Time_axis = 0 : hor;
Resp_mat = [];
for shock_counter = 1 : k_noises,
	Response = zeros(n_endog,hor);
	Response(n_endog,1) = H*E*P;
	for time_counter = 2 : hor,
		Response(:,time_counter) = H*Phi^(time_counter-1)*E*P;
	end;
	Resp_mat =	[ Resp_mat
			  Response ];
end;

irf = Resp_mat;

% Plotting

name = ['Impulse Response Function'];
%%% figure; whitebg('w'); close; %MBB
figure;
plot(Resp_mat);
grid;
title(name);
hold on;
