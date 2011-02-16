PRO fix_lambda, lambda

  IF n_params() LT 1 THEN BEGIN 
      print,'-Syntax: fix_lambda, lambda'
      return
  ENDIF 

  w=where(lambda GT 90.0d, nw)
  IF nw NE 0 THEN lambda[w] = 180.0d - lambda[w]

  w=where(lambda LT -90.0d, nw)
  IF nw NE 0 THEN lambda[w] = -180.0d - lambda[w]


END 
