function timeaxis = timefmt2 (P,Y,M,N)

  ## usage:  timeaxis = timefmt (P,Y,M,N)
  ##
  ## P > Periodicity (1 yearly, 4 quaterly, 12 monthly)
  ## Y > First year, ('1929')
  ## M > First month, ('3' for march)
  ## N > Number of dates 
  ##
  ## timeaxis < time axis

  mfy=(P+1-M);			# months of the first year
  middle_years=floor((N-mfy)/P); # number of middle years
  last_year=N-mfy-(middle_years*P);
  mly=Y+middle_years+1;		# months of the last year
  
  timeaxis=zeros(N,1);

  if P==4
    td=(M-1)*.25;
    for i=1:mfy
      timeaxis(i)=Y+td;
      td=td+.25;
    end
    
    j=mfy+1;
    for y=[Y+1:mly-1]
      for td=0:0.25:.75
	timeaxis(j)=y+td;
	j=j+1;
      end      
    end
    
    for td=1:last_year
      timeaxis(j-1+td)=mly+(td-1)*.25;
      td=td+1;
    end

  else
    td=M;
    for i=1:mfy
      timeaxis(i)=Y*100+td;
      td=td+1;
    end
    
    j=mfy+1;
    for y=[Y+1:mly-1]
      for td=1:P
	timeaxis(j)=y*100+td;
	j=j+1;
      end
    end
    
    for td=1:last_year
      timeaxis(j-1+td)=mly*100+td;
      td=td+1;
    end
  end
