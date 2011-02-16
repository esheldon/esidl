FUNCTION masssis_trunc, sigma, cutoff, radius, core=core

  IF n_params() LT 3 THEN BEGIN 
      print,'-Syntax: mass=masssis_trunc(sigma_v, cutoff, radius, core=core)'
      print,'cutoff/radius in kpc'
      print,'if /core then 1 kpc core is used'
      return,-1.
  ENDIF 

  IF keyword_set(core) THEN BEGIN 
      term1 = cutoff/sqrt(1. + (radius)^2 )
      term2 = sqrt( 1.0 + cutoff^2/(1.+radius^2) )
  ENDIF ELSE BEGIN 
      term1 = cutoff/radius
      term2 = sqrt(1.0+(cutoff/radius)^2)
  ENDELSE 

  return, 7.307e14*(sigma/1000.)^2*(radius/1000.)*(1.0 + term1 - term2 )
                                ;solar masses


END 
