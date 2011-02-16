PRO range2error, low, mean, high, errhigh, errlow

  IF n_params() LT 3 THEN BEGIN 
      print,'-Syntax: range2error, low, mean, high, errhigh, errlow'
      return
  ENDIF 

  errlow = mean-low
  errhigh = high-mean

END 
