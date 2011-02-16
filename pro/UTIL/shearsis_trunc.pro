FUNCTION shearsis_trunc, sigma, cutoff, zsource, zlens, angle, core=core

  IF n_params() LT 4 THEN BEGIN
      print,'Syntax:  s = shearsis_trunc(sigma, cutoff, zsource, zlens, angle, core=core)'
      print,'if /core, then add 1" core'
      print,'Sigma in km/s.   Angle in arcsec.   Assumes matter only universe'
      return,-1
  ENDIF 

  IF keyword_set(core) THEN docore = 1 ELSE docore = 0

  ratio = angdist(zsource,zlens,h=1.)/angdist(zsource,h=1.)

  IF docore THEN BEGIN
      term1 =  1./sqrt(1. + angle^2)
      term2 =  2.*cutoff*term1^2
      term3 = -1./sqrt( angle^2 + cutoff^2 )
      term4 =  2.*cutoff^2*term1^2*term3
  ENDIF ELSE BEGIN
      term1 =  1./angle
      term2 =  2.*cutoff*term1^2
      term3 = -1./sqrt( angle^2 + cutoff^2 )
      term4 =  2.*cutoff^2*term1^2*term3
  ENDELSE 

  return, 1.4*(sigma/220.)^2*ratio/2.*( term1 + term2 + term3 + term4 )

END 
      
