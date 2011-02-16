FUNCTION _fitsis_sismodel, radius, sigmav
  return, 3.36e3*(sigmav/170.)^2/radius ;Msolar/pc^2
END 

PRO fitsis, x, y, error, sigmav_guess, yfit, sigmav, sigmav_err

  itmax = 300
  FUNCNAME = '_fitsis_sismodel'
  sigmav = MPFITFUN(FUNCNAME, x, y, error, sigmav_guess, $
               AUTODERIVATIVE=1, MAXITER=itmax, niter=iter, $
               yfit=yfit, perror=sigmav_err, covar=covar,/quiet)

  IF NOT keyword_set(silent) THEN BEGIN 

      chi2 = total( ( (y-yfit)/error )^2, /double)
      dof = n_elements(y)-2
      chi2per = chi2/dof
      print
      print,'--------------------------------------------------------'
      print,'iterations = ',iter
      print,'sigmav     = ',ntostr(sigmav)+' +/- '+ntostr(sigmav_err)
      print,'chi^2/dof  = ',ntostr(chi2)+'/'+ntostr(dof)+' = '+ntostr(chi2per)
  ENDIF 


END 
