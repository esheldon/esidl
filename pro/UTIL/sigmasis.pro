FUNCTION sigmasis, sigma, radius, core=core

  IF n_params() LT 2 THEN BEGIN
      print,'Syntax:  s = sigdiffsis_trunc(sigma, radius, core=core)'
      print,'if /core, then add 1kpc core'
      print,'Sigma in km/s.   Radius in kpc.'
      return,-1
  ENDIF 

  IF keyword_set(core) THEN docore = 1 ELSE docore = 0

  IF docore THEN term = 1./sqrt(1. + radius^2) ELSE term = 1./radius

  return, 3.36e3*(sigma/170.)^2*term ;Msolar/pc^2

END 
