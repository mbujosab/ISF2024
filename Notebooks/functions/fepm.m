## usage:  [RMSE,MAE] = fepm (actval, fcast1, fcast2, ...)
##
## Forecast error percentage measures
##
## actval  * >  Actual values of variable to make comparisons (nx1).
## fcast1  * >  Forecasts 1 (nx1)
## fcast2    >  Forecasts 2 (nx1)
## ...       >  ...
##
## RMSE      <  Root Mean Squared Errors
## MAE       >  Mean Absulte Error
##
## Copyright (C) 2002  Marcos Bujosa

## Author: MB <marcos.bujosa@ccee.ucm.es>
## Description: Forecast error percentage measures.

#function [RMSE,MAE] = fepm (actval,fcast1,...)
function [RMSE,MAE] = fepm (actval,fcast)
  
 
  if size(actval,1)~=size(fcast,1)
    error('The number of actual and forecast values must be the same')
  endif

  if is_scalar(actval)
    RMSE = sqrt((actval-fcast).^2);
    MAE  = abs(actval-fcast);
  else
    actval=kron(actval,ones(1,size(fcast,2)));
    RMSE = sqrt(mean((actval-fcast).^2));
    MAE  = mean(abs(actval-fcast));
  endif

endfunction
