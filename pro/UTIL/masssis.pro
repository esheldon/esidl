FUNCTION masssis, sigma, radius

  IF n_params() LT 2 THEN BEGIN 
      print,'-Syntax: mass=masssis(sigma_v, radius)'
      return,-1.
  ENDIF 

  return, 7.307e14*(sigma/1000.)^2*(radius/1000.) ;solar masses


END 
