function [periodos,raices,normas,Wj,AR] = rootsana (AR)

  ## usage:  [periodos,raices,normas] = rootsana (AR)
  ##
  ## Ordenadas por normas de mayor a menor
  ## 
  ## Copyright (C) 2003, 2004, 2006 Marcos Bujosa  
  
  ## Author:  MB <marcos.bujosa@ccee.ucm.es> 
  ## Description:  ANÁLISIS DE LAS RAÍCES DEL POLINOMIO AR 

  if size(AR)==[1 1]
    periodos=[];raices=[];normas=[];Wj=[];
  else
    ## reducimos precisión para evitar errores por redondeo
    ##     precision=100000;		# suficiente para órdenes < 50
    ##     AR=real(round(AR*precision))/precision;

    raices=roots(AR); 		# raíces del polinomio AR 
    normas=real(abs(raices));	# norma de las raíces 
    Wj=angle(raices);

    ## Bug en octave 2.1.35 y 2.1.50
    Wj(abs(Wj)        <10e-8)= 0; raices(Wj== 0)=real(raices(Wj== 0));
    Wj(abs(abs(Wj)-pi)<10e-8)=pi; raices(Wj==pi)=real(raices(Wj==pi));

    periodos=zeros(size(Wj));	# periodos asociados 
    periodos(Wj~=0)=real((2*pi)./Wj(Wj~=0));
    periodos(Wj==0)=Inf;  	# para no dividir por cero

    ## ordenamos por normas 
    [kk,orden]=sort(abs(1-normas)); 
    periodos=periodos(orden);
    raices=raices(orden);
    normas=normas(orden);
    Wj=Wj(orden);
    
    ## ignoramos los periodos negativos 
    pos=find(periodos>0); 
    periodos=periodos(pos);
    raices=raices(pos);
    normas=normas(pos);
    Wj=Wj(pos);
  end

