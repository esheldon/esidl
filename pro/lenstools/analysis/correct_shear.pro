FUNCTION correct_shear, instruct, rstruct

  IF n_params() EQ 0 THEN BEGIN 
      print,'-Syntax: st = correct_shear(struct, rstruct)'
      return,-1
  ENDIF 

  narr = n_elements(instruct.rsum)
  arrval = fltarr(narr)
  struct = create_struct(instruct, $
                         'corr', arrval+1, $
                         'corr_err', arrval, $
                         'frac_clust', arrval, $
                         'frac_clust_err', arrval)

  ssh = struct.ssh

  corr = struct.wsum_mean/rstruct.wsum_mean; > 1.0
  corr_err = corr*sqrt( (struct.wsum_err/struct.wsum_mean)^2 + $
                        (rstruct.wsum_err/rstruct.wsum_mean)^2 )
  
  frac_clust = 1.0 - rstruct.wsum_mean/struct.wsum_mean 
  frac_clust_err = frac_clust/corr*corr_err


  struct.sigma = struct.sigma/ssh*(corr > 1.0)
  struct.sigmaerr = struct.sigmaerr/ssh*(corr > 1.0)

  struct.corr = corr
  struct.corr_err = corr_err
  struct.frac_clust = frac_clust
  struct.frac_clust_err = frac_clust_err

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; total within radius
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;
  
  wts = 1./struct.sigmaerr^2
  tsigma_sum = total(struct.sigma*wts, /cumulative)
  twsum = total(wts, /cumulative)

  struct.tsigma = tsigma_sum/twsum
  struct.tsigmaerr = sqrt(1./twsum)

  return, struct
      
END 
