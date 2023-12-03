function [theta,din,lab,AR,MA,s2a] =  dhrgraph (y,TVP,P,VAR,p,oar)
  
  nP=size(P,2);

  [theta,din,lab,AR,MA,s2a] = dhr2rf (TVP,P,VAR,[],p);

  xx=(filter(AR,1,y)); 
  xx=xx(size(AR,2):size(y,1));
  rx=size(xx,1);  
  nfrec=[1:rx-1]/rx*2*pi;
  Pxx = (abs (fft (xx - mean (xx)))) .^ 2 / rx; % Periodograma
  Pxx(1)=[]; 
    
  X=nfrec(ones(length(P),1),:)';
  frecMOD=2*pi./P;
  XX=frecMOD(ones(length(nfrec),1),:);
  [f,c]=find((abs(X-XX)<2*eps));
  if ~isempty(f), nfrec(f)=[];Pxx(f)=[]; end
  f=nfrec>pi;Pxx(f)=[];nfrec(f)=[];
  [ares,s2e]=aresp(y,oar,[],[],nfrec);

  %xlab=sprintf('frequencies [0\\ \\pi]');
  xlab=sprintf('\\omega');
  ylab=sprintf('f(\\omega) (Log scale)');
  leyen=['Periodogram        ';
  	 'DHR pseudo-spectrum'
  	 'AR-spectrum        '];
  titulo = sprintf('Periodogram, DHR pseudo-spectrum and AR(%g)-spectrum',oar);
  
  %% periodograma - modelo DHR - ARespectro
  figure; 
  semilogy(nfrec,[Pxx.*esptarma(1,AR,nfrec,1),esptarma(MA,AR,nfrec,s2a),ares*s2e])
  %legend(leyen,1);
  set (gcf(), 'DefaultTextFontName', 'arial') 
  set (gca,'XTick',0:pi/2:pi)
  set (gca,'XTickLabel',{'0','pi/2','pi'})
  title(titulo); xlabel(xlab); ylabel(ylab)
    
  %% periodograma (MA) - modelo DHR (MA)
  figure; semilogy(nfrec,[Pxx,esptarma(MA,1,nfrec,s2a)])
  set (gcf(), 'DefaultTextFontName', 'arial') 
  set(gca,'XTick',0:pi/2:pi)
  set(gca,'XTickLabel',{'0','pi/2','pi'})
  title('Filtered (MA stationary part) periodogram and DHR spectrum'); xlabel(xlab); ylabel(ylab)
  
  
				%     %% ganancia filtro de la tendencia
				%     ylab=sprintf('f(w)');
				%     if any(P==Inf)
				%       nfrec=[1:512-1]/512*pi;
				%       [thetat,dint,labt,ARt,MAt,s2at] = dhr2thd (TVP(:,1),P(1),VAR([1,2]),[],p); 
				%       [gaint,phaset]=gainandphase(ARt,MAt,nfrec);
				%       figure; plot(nfrec,gaint)    
				%       set(gca,'XTick',0:pi/2:pi)
				%       set(gca,'XTickLabel',{'0','pi/2','pi'})
				%       title('Gain function of trend-extraction filter'); xlabel(xlab); ylabel(ylab)
				%       flecha=sprintf('\\leftarrow %5.0f periodos',2*pi/min(nfrec(gaint>=.5)));
				%       text(min(nfrec(gaint>=.5)+.03),.5,flecha,'HorizontalAlignment','left')
				%     end
  
  

 % figure; plot(nfrec,[Pxx,esptarma(MA,1,nfrec,s2a)])

  
  %% ganancia filtro desestacionalizador
  ylab=sprintf('f(\\omega)');
  if p>1 & any(P(2:nP))
    nfrec=[1:512-1]/512*pi;
    [thetaus,dinus,labus,ARus,MAus,s2aus] = dhr2rf (TVP(:,2:nP),P(2:nP),VAR([1,3:nP+1]),[],p); 
				%      figure; plot(nfrec,esptarma(ARus,MAus,nfrec,1)./max(esptarma(ARus,MAus,nfrec,1)))
    %figure; plot(nfrec,gainandphase(ARus,MAus,nfrec,1))
    figure; plot(nfrec,gainandphase(ARus,MAus,nfrec))
    set (gcf(), 'DefaultTextFontName', 'arial') 
    set(gca,'XTick',0:pi/2:pi)
    set(gca,'XTickLabel',{'0','pi/2','pi'})
    title('Gain function of unseasonal series-extraction filter'); xlabel(xlab); ylabel(ylab)
  end
  
  
  %% ganancia filtro DHR completo
  [gain,phase]=gainandphase(AR,MA,nfrec);
  figure; plot(nfrec,gain)    
  set (gcf(), 'DefaultTextFontName', 'arial') 
  set(gca,'XTick',0:pi/2:pi)
  set(gca,'XTickLabel',{'0','pi/2','pi'})
  title('Gain function of DHR filter'); xlabel(xlab); ylabel(ylab)
  flecha=sprintf('\\leftarrow %5.2f periodos',2*pi/min(nfrec(gain>=.5)));
  text(min(nfrec(gain>=.5)+.03),.5,flecha,'HorizontalAlignment','left')
  
				%      for i=1:size(P,2)
				%        [theta,din,lab,AR,MA,s2a] = dhr2thd (TVP(:,i),P(i),VAR([1,1+i]),[],p);
				%        figure(i+3);
				%        plot(nfrec,gainandphase(AR,MA,nfrec,1))
				%      end
  
  drawnow ()
 
