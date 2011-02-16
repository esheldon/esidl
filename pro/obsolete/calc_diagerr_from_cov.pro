PRO calc_diagerr_from_cov, covariance, diagerr

  IF n_params() LT 2 THEN BEGIN 
      print,'-Syntax: calc_diagerr_from_cov, covariance, diagerr'
      return
  ENDIF 

  sz = size(covariance)

  nbin = sz[1]

  diagerr = dblarr(nbin)

  FOR i=0L, nbin-1 DO diagerr[i] = sqrt(covariance[i, i])

END 
