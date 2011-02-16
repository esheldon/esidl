
FUNCTION mean_check, X, Double = Double, NaN = nan

  ON_ERROR, 2

  IF n_elements(X) LT 2 THEN return,X[0] ELSE $
  RETURN, (moment( X, Double=Double, Maxmoment=1, NaN=nan ))[0]
END
