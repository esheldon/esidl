;; old procedure version.  See cov2corr()
PRO calc_correlation_matrix, cov, corr, status

  IF n_params() LT 1 THEN BEGIN 
      print,'-Syntax: calc_correlation_matrix, cov, corr, status'
      print,' status = 0 is good'
      return
  ENDIF 

  corr = cov2corr(cov, status=status)
  return

END 
