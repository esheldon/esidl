;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; Function to compute Mallow's C_p, an estimate of the chi-square
; between the true function and the fit. This function assumes that
; both y and yfit have been normalized by the covariance matrix of y.
;
; Author : Brandon C. Kelly, Steward Obs., Feb 2005
;
; Inputs :
;    Y - The data, normalized by the covariance matrix so that it is
;        normally distributed with covariance equal to the identity
;        matrix.
;    YFIT - The fit to Y, normalized by the covariance matrix of Y.
;    M - The number of parameters used in the linear model :
;                  YFIT = S ## Y ,
;        For S some matrix that converts Y to YFIT.
;
; Output :
;    The mean-squared error of YFIT with the true function that
;    generated, i.e., the mean of the distribution that generated Y.
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function cp_mse, y, yfit, m

  if n_params() ne 3 then begin
      print, 'Syntax- C_p = cp_mse( y, yfit, m )'
      return, -1
  endif

  n = n_elements(y)
  cp = mean( (y - yfit)^2, /nan ) - 1 + 2d * m / n
  
  return, cp
end

