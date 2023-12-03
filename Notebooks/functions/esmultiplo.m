function t = esmultiplo (m,f)

  ## usage:  t = esmultipo (f,m)
  ##
  ## 1 si m=f*(algun entero)

  t=any(factor(m)==f);

endfunction