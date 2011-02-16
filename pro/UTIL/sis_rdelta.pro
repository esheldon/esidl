FUNCTION sis_rdelta, delta, z, sigma

  IF n_params() LT 3 THEN BEGIN 
      print,'Syntax: result=sis_rdelta(delta, z, sigma)'
      print,'  sigma in km/s.  returns r in Mpc'
      return,-1.
  ENDIF 


  return, 20.0/sqrt(delta)/(1. + z)^(3./2.)*(sigma/1000.)

END 
