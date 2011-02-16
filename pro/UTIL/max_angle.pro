FUNCTION max_angle, z, rmax, h=h, omegamat=omegamat

  IF n_params() LT 2 THEN BEGIN 
      print,'-Syntax: result = max_angle( z, rmax, h=h, omegamat=omegamat)'
      return,-1.
  ENDIF 

  return, rmax/angdist_lambda(z, h=h, omegamat=omegamat)*180.0/!pi

END 
