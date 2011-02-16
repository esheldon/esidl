FUNCTION sigdiffsis_trunc, sigma, cutoff, radius, core=core

  IF n_params() LT 3 THEN BEGIN
      print,'Syntax:  s = sigdiffsis_trunc(sigma, cutoff, radius, core=core)'
      print,'if /core, then add 1kpc core'
      print,'Sigma in km/s.   Radius/cutoff in kpc.'
      return,-1
  ENDIF 

  IF keyword_set(core) THEN docore = 1 ELSE docore = 0

  IF docore THEN BEGIN
      term1 =  1./sqrt(1. + radius^2)
      term2 =  2.*cutoff*term1^2
      term3 = -1./sqrt( radius^2 + cutoff^2 )
      term4 =  2.*cutoff^2*term1^2*term3
  ENDIF ELSE BEGIN
      term1 =  1./radius
      term2 =  2.*cutoff*term1^2
      term3 = -1./sqrt( radius^2 + cutoff^2 )
      term4 =  2.*cutoff^2*term1^2*term3
  ENDELSE 

  return, 3.3e3*(sigma/170.)^2*( term1 + term2 + term3 + term4 ) ;Msolar/pc^2

END 
      
