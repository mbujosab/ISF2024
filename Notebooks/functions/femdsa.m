## usage:  [out] = femdsa (act, fcast)
##
## Forecast error measures for different steps ahead
##
## act   * >  Actual values               (Nx1)
## fcast * >  Forecasts (1:p steps ahead) (Nxp)
##
## out     <  is the matrix [RMSE, MAE, RMSPE, MAPE]
##  
##   where columns are:
##
## RMSE    :  Root Mean Squared Errors
## MAE     :  Mean Absulte Error
## RMSPE   :  Root Mean Squared Percentage Errors
## MAPE    :  Mean Absulte Percentage Error
##
##   and the j-th row corresponds to j-steps ahead forecasts

## Copyright (C) 2008  Marcos Bujosa

## Author: MB <marcos.bujosa@ccee.ucm.es>
## Description: Forecast error measures.

function out = femdsa (act,fcast)
  
 
  if size(act,1)~=size(fcast,1)
    error('The number of actual and forecast values must be the same')
  endif

  p=size(fcast,2);

  out=zeros(s,4);

  for s=0:p-1
    out(s+1,:)=fem(act(s+1:end),diag(fcast,-s));
  endfor

endfunction
