PRO rotate_lambda, lambda

  IF n_params() LT 1 THEN BEGIN 
      print,'-Syntax: rotate_lambda, lambda'
      return
  ENDIF 

  wL = where(lambda LT 0., nL)
  wG = where(lambda GE 0., nG)

  IF nL NE 0 THEN BEGIN 
      lambda[wL] = lambda[wL] + 180d
  ENDIF 
  IF nG NE 0 THEN BEGIN 
      lambda[wG] = lambda[wG] - 180d
  ENDIF 

  return
END
