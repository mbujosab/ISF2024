function [Phi, Gam, E, H, D, C, Q, S, R] = pcons1(thetan, dinn);
% SS formulation of the tf model with period as explicit parameter
% stored in theta(2,1).
theta = thetan;
theta(2,1) = (-theta(1,1)/(2*cos(2*pi/thetan(2,1))))^2;
din   = tomod(dinn);
[Phi, Gam, E, H, D, C,Q, S, R]= thd2ss(theta,din);
