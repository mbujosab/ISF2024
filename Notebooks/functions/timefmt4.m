function Y = timefmt4 (P,Y,M,N)

  ## usage:  Y = timefmt4 (P,Y,M,N)
  ##
  ## P > Periodicity (1 yearly, 4 quaterly, 12 monthly)
  ## Y > First year, ('1929')
  ## M > First month, ('3' for march)
  ## N > Number of dates 
  ##
  ## Y < time axis

  mfy=(P+1-M);			# months of the first year
  middle_years=floor((N-mfy)/P); # number of middle years
  last_year=N-mfy-(middle_years*P);
  mly=Y+middle_years+1;		# months of the last year

  timeaxis= repmat (blanks(12), N, 1);
  td=M;  
  for i=1:mfy
      a=strjust ([num2str(td);"AA"]);  
      timeaxis(i,:)=cat(2,num2str(Y),"-",strrep (a(1,:)," ","0"),"-01  ") ; 
      td=td+1;
  end
  
  j=mfy+1;
  for y=[Y+1:mly-1]
      for td=1:P
          a=strjust ([num2str(td);"AA"]);  
          timeaxis(j,:)=cat(2,num2str(y),"-",strrep (a(1,:)," ","0"),"-01  "); 
          j=j+1;
      end
  end
  
  for td=1:last_year
      a=strjust ([num2str(td);"AA"]);  
      timeaxis(j-1+td,:)=cat(2,num2str(mly),"-",strrep (a(1,:)," ","0"),"-01  ");
      td=td+1;
  end
  
  Y = cat(1,"obs         ",timeaxis);
  