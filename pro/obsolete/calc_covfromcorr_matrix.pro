PRO calc_covfromcorr_matrix, corr, diagerr, cov

  IF n_params() LT 1 THEN BEGIN 
      print,'-Syntax: calc_covfromcorr_matrix, corr, diagerr, cov'
      return
  ENDIF 

  cov = corr2cov(corr, diagerr)
  return
END 
