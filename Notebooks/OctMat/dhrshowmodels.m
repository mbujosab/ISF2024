function dhrshowmodels (PaP,mostrados,model,mcnn,R2,NVR,oar)

  %% mostrando los modelos alternativos
  cabecera= '\n%s  %s %s     ';
  if any(PaP==Inf), cabecera=strcat (cabecera,'%s         '); end
  
  if any(PaP(1,:)<inf&PaP(1,:)>=100)
    strtab='%3.1f     ';
    cabecera=strcat (cabecera,sprintf('%c',vec(strtab' *ones(size(find(PaP(1,:)<inf&PaP(1,:)>=100))))));
  end
  if any(PaP(1,:)<inf&PaP(1,:)>=10)
    strtab='%3.1f      ';
    cabecera=strcat (cabecera,sprintf('%c',vec(strtab' *ones(size(find(PaP(1,:)<inf&PaP(1,:)>=10))))));
  end
  if any(PaP(1,:)<10)
    strtab='%3.1f       ';
    cabecera=strcat (cabecera,sprintf('%c',vec(strtab' *ones(size(find(PaP(1,:)<10))))));
  end
  
  if any(PaP==Inf)
    fprintf(cabecera,'AR', 'NN', 'R2', 'T', PaP(PaP(1,:)<inf))
  else
    fprintf(cabecera,'AR', 'NN', 'R2', PaP(PaP(1,:)<inf))
  end
  
  tabla='\n%2.0f  %1.0f % 2.3f  ';
  fff='%0.4f-%c  ';
  tabla=strcat (tabla,sprintf('%c',vec(fff' *ones(size(find(PaP(1,:)))))));
  
  for i=1:size(mostrados,1)
    Modelo=model(mostrados(i)).TVP;
    L=(' '-0)*ones(size(Modelo(1,:)));
    L(all(Modelo==1))=('I'-0);
    L((Modelo(1,:)==1)&(Modelo(2,:)<1))=('S'-0);
    L((Modelo(1,:)==1)&(Modelo(2,:)==0))=('R'-0);
    L((Modelo(1,:)<1))=('A'-0);
    L(all(Modelo==0))=('-'-0);
    s2i=model(mostrados(i)).VAR(1);
    if mostrados(i) == oar&&i>1;fprintf('\n');end
    fprintf(tabla,[mostrados(i),mcnn(mostrados(i)),R2(mostrados(i)),reshape([model(mostrados(i)).VAR(2:end)./s2i;L],1,size(NVR,2)*2)])
    if mostrados(i) == oar;fprintf('\n');end
  end
  fprintf('\n\n');
  
 
