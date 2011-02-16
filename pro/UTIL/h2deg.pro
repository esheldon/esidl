FUNCTION h2deg, hr, min, sec

  IF n_params() LE 1 THEN BEGIN 
      print,'Syntax: h2deg, hr, min, sec, deg'
      return, 0.d
  ENDIF 

  np = n_params()

  IF np EQ 1 THEN return, hr*15.d
  IF np EQ 2 THEN return, hr*15.d + min/60.*15.d
  return, hr*15.d + min/60.*15.d + sec/3600.*15.d

END 
