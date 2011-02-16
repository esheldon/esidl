PRO reconstruct_resamp_cov, boot_covmean, boot_coverr, Nvar, cov, coverr

  IF n_params() LT 2 THEN BEGIN 
      print,'-Syntax: reconstruct_resamp_cov, boot_covmean, boot_coverr, cov, coverr'
      return
  ENDIF 

  ;; The point here is to take the diagonal elements from
  ;; bootstrap, which are actually the variance and covariance,
  ;; and write them into a covariance matrix on the terms.
  ;; Same with the means

  cov = replicate(boot_covmean[0], nvar, nvar)
  coverr = cov

  FOR j=0L, nvar-1 DO BEGIN 
      FOR i=j, nvar-1 DO BEGIN 

          cov[i,j] = boot_covmean[j*Nvar + i]
          coverr[i,j] = boot_coverr[j*Nvar + i]

          IF i NE j THEN BEGIN 
              cov[j,i] = cov[i,j]
              coverr[j,i] = coverr[i,j]
          ENDIF 

      ENDFOR 
  ENDFOR 

END 
