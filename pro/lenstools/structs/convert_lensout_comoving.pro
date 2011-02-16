PRO convert_lensout_comoving, struct

  ;; This converts from the mean redshift to redshift zero (comoving)

  IF n_params() LT 1 THEN BEGIN 
      print,'-Syntax: convert_lensout_comoving, struct'
      return
  ENDIF 

  factor = 1 + struct.zmean
  
  struct.rmin = struct.rmin*factor
  struct.rmax = struct.rmax*factor
  
  struct.rmin_act = struct.rmin_act*factor
  struct.rmax_act = struct.rmax_act*factor
  
  struct.rsum = struct.rsum*factor
  struct.meanr = struct.meanr*factor
  
  struct.sigma = struct.sigma/factor^2
  struct.sigmaerr = struct.sigmaerr/factor^2
  struct.sigmaerr2 = struct.sigmaerr2/factor^2
  
  struct.orthosig = struct.orthosig/factor^2
  struct.orthosigerr = struct.orthosigerr/factor^2
  struct.orthosigerr2 = struct.orthosigerr2/factor^2
  
  struct.tsigma = struct.tsigma/factor^2
  struct.tsigmaerr = struct.tsigmaerr/factor^2
  struct.tsigmaerr2 = struct.tsigmaerr2/factor^2
  
  struct.sigerrsum = struct.sigerrsum/factor^2
  struct.orthosigerrsum = struct.orthosigerrsum/factor^2
  
  IF tag_exist(struct, 'covariance') THEN $
    struct.covariance = struct.covariance/factor^4

  return

END 
