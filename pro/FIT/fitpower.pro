PRO fitpower, x, y, error, aguess, yfit, a, sigma, $
              silent=silent

  IF N_params() EQ 0 THEN BEGIN 
     print,'-Syntax: fitpower, x, y, error, aguess [,yfit, a, sigma, silent=silent]'
     print,''
     print,'aguess = [norm, power]'
     print,'Use doc_library,""  for more help.'  
     return
  ENDIF 
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Parameters
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  a = aguess                    ;because it will change.
  na = n_elements(a)
  IF na NE 2 THEN BEGIN
      print,'a guess array must have 2 elements'
      return
  ENDIF 
  nx = n_elements(x)
  ny = n_elements(y)
  IF ny NE nx THEN BEGIN
      print,'X and Y must be same size'
      return
  ENDIF 

  w=where(error EQ 0., nw)
  IF nw NE 0 THEN BEGIN
      print,'Errors must be non-zero'
      return
  ENDIF 

  itmax = 300
  MYFUNCT = 'mpfit_power_law'
  a = MPFITFUN(MYFUNCT, x, y, error, aguess, $
               AUTODERIVATIVE=1, MAXITER=itmax, niter=iter, $
               yfit=yfit, perror=sigma, covar=covar,/quiet)
  
  IF NOT keyword_set(silent) THEN BEGIN 

      chi2 = total( ( (y-yfit)/error )^2, /double)
      dof = n_elements(y)-2
      chi2per = chi2/dof
      print
      print,'--------------------------------------------------------'
      print,'iterations = ',iter
      print,'Norm  (a0) = ',ntostr(a[0])+' +/- '+ntostr(sigma[0])
      print,'Index (a1) = ',ntostr(a[1])+' +/- '+ntostr(sigma[1])
      print,'chi^2/dof  = ',ntostr(chi2)+'/'+ntostr(dof)+' = '+ntostr(chi2per)
  ENDIF 



  return 
END 
