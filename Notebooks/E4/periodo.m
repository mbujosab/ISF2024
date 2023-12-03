function [maxper,freq,power] = periodo(y,tit)
% Computes the periodogram of a time series
Y = fft(y);
n = length(Y);
Y(1) = [];
n=length(Y);

power = abs(Y(1:n/2)).^2;
nyquist = 1/2;
freq = ((1:n/2)/(n/2)*nyquist)';

figure;
%whitebg('w');
hold on
plot(freq,power,'k-')
xlabel('cycles/unit time')
if nargin == 2
  title(['Periodogram of ' tit])
end
grid
hold off

period=(1./freq)';
index=find(power==max(power));
maxper=(period(index));

return

Pyy=Y.*conj(Y)/n;
fre=1000/n*(0:(n/2-1));

figure;
%whitebg('w');
hold on
plot(fre,Pyy(1:n/2),'k-')
xlabel('Hz')
if nargin == 2
  title(['Spectral density of series ' tit])
end
grid
hold off
