function [fri,gain] = irf(thopt,din,maxlag,titles)

[Phi, Gam, E, H, D, C, Q, S, R] = thd2ee(thopt,din)

% Impulse-respose
n_endog=eye(size(z,2))

fri=[eye(n_endog) H*E];
for i=1:(maxlag-2)
  fri=[fri H*Phi^i*E];
end

figure; whitebg('w'); close;

figure
hold on
plot(0:(maxlag-1),fri(1,1:2:80),'k-')
grid
xlabel('lag')
ylabel('impulse response')
title('Response of mink skins to a unit pulse in mink skins')
axis([0 maxlag -1 1])

gain1=sum(fri(1,1:2:80))

