function [Phi, Gam, E, H, D, C, Q, S, R] = Scycle(theta, din)
% Reparametrization of theta for building stochastic cycle
% theta(1,1) => ro
% theta(2,1) => lambda
% theta(5,1) => sigma_k

ro = theta(1,1);
lambda = theta(2,1);
sigma_k = theta(6,1);

theta(1,1) = ro*cos(lambda);
theta(4,1) = theta(1,1);
theta(3,1) = ro*sin(lambda);
theta(2,1) = -theta(3,1);
theta(7,1) = sigma_k;

din = tomod(din);

[Phi, Gam, E, H, D, C, Q, S, R] = thd2ss(theta, din);
