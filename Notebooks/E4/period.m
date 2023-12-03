function [powr,freq,seasfreq,tdfreq] = period(z,s,tit);
% arguments: z data, s seasonal period

n=size(z,1); 		          	% Sample size
powr=zeros(61,1); 		      % Power
a=zeros(61,1); b=zeros(61,1);	% Coefficients of the harmonic regression
freq=(0:60)'./120;	      	% Frequencies from 0 to 1/2 cycles/year
seq=1:n;

if s ~= 1
% Computes a seasonal indicator: 1 seasonal frequency, 0 otherwise
  seasfreq=NaN*zeros(size(freq));
  seasf=[1:fix(s/2)]'/s;
  for i=1:size(seasf,1)
    seasfreq(find(freq==seasf(i)))=1;
  end

  if s == 12
% Computes an indicator of trading day effects: 1 trading day frequency, 0 otherwise
    tdfreq=NaN*zeros(size(freq));
    tdf=[0.304 0.348 0.432];
    j=1;
    for k=2:61			%Inserts tdf into the range of frequencies
      if freq(k-1) < tdf(j) & freq(k) >= tdf(j)
        if abs(freq(k) - tdf(j)) < abs(freq(k-1)-tdf(j))
          freq(k)=tdf(j);
          tdfreq(k)=1;
        else
          freq(k-1)=tdf(j);
          tdfreq(k-1)=1;
        end
        j=j+1;
      end
      if j > 3; break; end
    end
  end
else
  seasfreq=[];
end

a(1)=mean(z);
k=0;
for k=2:60
  a(k) = (2/n)*cos(2*pi*freq(k)*seq)*z;
  b(k) = (2/n)*sin(2*pi*freq(k)*seq)*z;
end

a(61)=(1/n)*(cos(2*pi*freq(61)*seq)*z);
powr=(n/2*(a.^2+b.^2));
powr(1)=2*powr(1);
powr(61)=2*powr(61);

% Plot the results
figure
%whitebg('w');
hold on
plot(freq,powr,'k-')

if s ~= 1
  plot(freq,powr.*seasfreq,'ko')
  if s == 12
    plot(freq,powr.*tdfreq,'kx')
  end
end

grid
xlabel('Frequency (cycles/unit time)')
ylabel('Power')
ax=axis;
axis([0 .5 0 ax(4)])
if nargin == 3
  title(['Periodogram of ' tit])
end

index=find(powr==max(powr));
MaxFreqStr=num2str(freq(index));
% plot(freq(index),powr(index),'k.','MarkerSize',26,'EraseMode','none')
% text(freq(index)+2,powr(index),['Max. Freq = ' MaxFreqStr],'EraseMode','none')

hold off