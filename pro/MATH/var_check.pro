
FUNCTION var_check, X, Double = Double, NaN = NaN

  ON_ERROR, 2

  IF n_elements(x) LT 2 THEN return,0. ELSE $
  RETURN, (moment( X, Double=Double, Maxmoment=2, NaN = NaN ))[1]
END
