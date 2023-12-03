## usage:  [out] = fem (act, fcast)
##
## Forecast error measures
##
## act   * >  Actual values (Nx1)
## fcast * >  Forecasts     (Nx1)
##
## out     <  [RMSE, MAE, RMSPE, MAPE]
##  
##   where:
##
## RMSE    :  Root Mean Squared Errors
## MAE     :  Mean Absulte Error
## RMSPE   :  Root Mean Squared Percentage Errors
## MAPE    :  Mean Absulte Percentage Error
##
## Copyright (C) 2008  Marcos Bujosa

## Author: MB <marcos.bujosa@ccee.ucm.es>
## Description: Forecast error measures.

function [out] = fem (act,fcast)
 
  if size(act,1)~=size(fcast,1)
    error('The number of actual and forecast values must be the same')
  endif

  err=act-fcast;

  RMSE  = sqrt(mean(err.^2));
  MAE   = mean(abs(err));
  RMSPE = sqrt(mean(err./act.^2));
  MAPE  = mean(abs(err./act));

  out=[RMSE,MAE,RMSPE,MAPE];

endfunction
