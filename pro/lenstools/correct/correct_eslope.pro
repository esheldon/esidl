PRO correct_eslope, fitStruct, $
                    r, e1, e2, psfe1, psfe2, $
                    e1c, e2c, $
                    original=original, $
                    nsmooth=nsmooth

  ;; nsmooth = 5 works well
  IF n_params() LT 8 THEN BEGIN 
      print,'-Syntax: correct_eslope, fitStruct, r, e1, e2, psfe1, psfe2, nsmooth, e1out, e2out, /original, nsmooth='
      return
  ENDIF 

  IF n_elements(nsmooth) EQ 0 THEN nsmooth=5

  ;; Interpolate to the r-value
  int_ge1pe1 = $
    interpol( smooth(fitStruct.intercept_ge1pe1, nsmooth), fitStruct.mean_r, r)
  slope_ge1pe1 = $
    interpol( smooth(fitStruct.slope_ge1pe1, nsmooth), fitStruct.mean_r, r)
  
  int_ge2pe2 = $
    interpol( smooth(fitStruct.intercept_ge2pe2, nsmooth), fitStruct.mean_r, r)
  slope_ge2pe2 = $
    interpol( smooth(fitStruct.slope_ge2pe2, nsmooth), fitStruct.mean_r, r)

  ;; Correct the shapes

  e1c = e1 - (int_ge1pe1 + slope_ge1pe1*psfe1)
  e2c = e2 - (int_ge2pe2 + slope_ge2pe2*psfe2)

END 
