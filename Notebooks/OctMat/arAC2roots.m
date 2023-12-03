function  [PAC,RAC,rAC,nAC,wAC,pAC]       =  arAC2roots (aR)

  if length(aR)==1; PAC=[]; RAC=[]; return; end 

  rAC        = roots(aR);
  nAC        = abs(rAC);
  precision  = 10^15; 
  wAC        = fix(precision*angle(rAC))/precision;  
  pAC        = inf*ones(size(wAC));
  pAC(wAC~=0)= round(precision*2*pi./wAC(wAC~=0))/precision;
  wAC        = 2*pi./pAC;

  x   = fix((cos (wAC) .* nAC)*10^16)/10^16;
  y   = fix((sin (wAC) .* nAC)*10^16)/10^16;
  rAC = x+y*i;

  [kk, ind] = sort (abs(pAC),'descend');
  pAC = pAC(ind);
  rAC = rAC(ind);
  nAC = nAC(ind);
  wAC = wAC(ind);

  %% calculo de P 
  pp = pAC(pAC>0);		% solo periodos positivos

  %% cuales hay (descontando repeticiones)
  f = find([pp~=circshift(pp,1)]); f(1)=1;
  TU= triu(single(pp(f,ones(length(f),1))==(pp(f,ones(length(f),1)))'),1);
  c = find(TU*f);              f(c)=[];
  
  PAC = pp(f)';

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% calculo de RAC (raices AR) %t TVP
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  RAC = zeros(4,size(PAC,2));

  p  = [pAC;NaN*ones(20,1)];
  ii = zeros(4,size(PAC,2));
  tt = zeros(2,size(PAC,2));

%%%%%%%%%%%%%%%%%%%%%%
  ff = find(abs(p)~=circshift(abs(p),-1)&abs(p)~=circshift(abs(p),1)&~isnan(p))'; % solo 1 raiz (necesariamente real)

  [kk,col]=size(ff);

  fr1=zeros(4,col);
  fr1(1:kk,:)=ff;

  tr1=zeros(2,col);
  tr1(1:kk,:)=ff;

  pr1=p(fr1(1,:))';
  [kk,col]=find(PAC(ones(size(fr1,2),1),:)==pr1(ones(size(PAC,2),1),:)');

  ii(:,col)=fr1;
  tt(:,col)=tr1;

  %%%%%%%%%%%%%%%%%%%%%%
 
  ff=[find(p==circshift(p,-1))';
      find(p==circshift(p, 1))'];	% 2 raices reales
  [kk,col]=size(ff);

  fr2=zeros(4,col);
  fr2(1:kk,:)=ff;

  tr2=zeros(2,col);
  tr2(1:kk,:)=ff;

  pr2=p(fr2(1,:))';
  [kk,col]=find(PAC(ones(size(fr2,2),1),:)==pr2(ones(size(PAC,2),1),:)');

  ii(:,col)=fr2;
  tt(:,col)=tr2;

%%%%%%%%%%%%%%%%%%%%%%

  ff=[find(p==circshift(-p,-1)&p~=circshift(-p, 1)&p~=circshift(-p,-3))';
      find(p==circshift(-p, 1)&p~=circshift(-p,-1)&p~=circshift(-p, 3))']; % un par de raices conjugadas

  [kk,col]=size(ff);

  fc1=zeros(4,col);
  fc1(1:kk,:)=ff;

  f=[find(p==circshift(-p,-1)&p~=circshift(-p, 1)&p~=circshift(-p,-3))'];
  [kk,col]=size(f);

  tc1=zeros(2,col);
  tc1(1:kk,:)=f;

  pc1=p(fc1(1,:))';
  [kk,col]=find(PAC(ones(size(fc1,2),1),:)==pc1(ones(size(PAC,2),1),:)');

  ii(:,col)=fc1;
  tt(:,col)=tc1;

%%%%%%%%%%%%%%%%%%%%%%

  ff=[find(p==circshift(-p,-1)&p==circshift(-p,-3))';
      find(p==circshift(-p, 1)&p==circshift( p,-2))';
      find(p==circshift(-p,-1)&p==circshift( p, 2))';
      find(p==circshift(-p, 1)&p==circshift(-p, 3))']; % 2 pares de raices conjugadas

  [kk,col]=size(ff);

  fc2=zeros(4,col);
  fc2(1:kk,:)=ff;


  f=[find(p==circshift(-p,-1)&p==circshift(-p,-3))';
     find(p==circshift(-p,-1)&p==circshift( p, 2))'];

  [kk,col]=size(f);
  
  tc2=zeros(2,col);
  tc2(1:kk,:)=f;
  
  pc2=p(fc2(1,:))';
  [kk,col]=find(PAC(ones(size(fc2,2),1),:)==pc2(ones(size(PAC,2),1),:)');
  
  ii(:,col)=fc2;
  tt(:,col)=tc2;

%%%%%%%%%%%%%%%%%%%%%%

  RAC(find(ii))=rAC(ii(ii>0));

 
