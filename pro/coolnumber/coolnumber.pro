PRO coolnumber, nn

  IF n_params() LT 1 THEN BEGIN 
      print,'-Syntax: coolnumber, nprint'
      return
  ENDIF 

  cn = 999999L/7

  FOR i=1, nn DO print,i,i*cn


END 
