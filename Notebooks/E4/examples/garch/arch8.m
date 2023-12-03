function[Phi, Gam, E, H, D, C, Q, Phig, Gamg, Eg, Hg, Dg] = arch8(thetan,dinn);
% User function for constrained ARCH(8) model

theta = zeros(14,1);

theta(1:6) = thetan(1:6,1);
alpha1     = thetan(7,1);

% Now constrained parameters
for i=1:8
    theta(i+6) = alpha1*(9-i)/36;
end

din= tomod(dinn);
[Phi Gam E H D C Q Phig Gamg Eg Hg Dg]= garch2ss(theta,din);
