function fit_r200power, rmpc, rhobyrhocrit, error, aguess, $
                        silent=silent

  IF N_params() EQ 0 THEN BEGIN 
      on_error, 2
      print,'-Syntax: fitst = fit_r200power(r, rhobyrhocrit, error, aguess, /silent)'
      print,''
      print,'aguess = [r200, power]'
      print,'rmpc, r200 in Mpc'
      print,'Use doc_library,""  for more help.'  
      message,'Halting'
  ENDIF 
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Parameters
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  a = aguess                    ;because it will change.
  na = n_elements(a)
  IF na NE 2 THEN message,'a guess array must have 2 elements'
  nrmpc = n_elements(rmpc)
  nrhobyrhocrit = n_elements(rhobyrhocrit)
  IF nrhobyrhocrit NE nrmpc THEN message,'RMPC and RHOBYRHOCRIT must be same size'

  w=where(error EQ 0., nw)
  IF nw NE 0 THEN print,'Errors must be non-zero'

  itmax = 300
  MYFUNCT = 'mpfit_rhobyrhocrit'
  a = MPFITFUN(MYFUNCT, rmpc, rhobyrhocrit, error, aguess, $
               AUTODERIVATIVE=1, MAXITER=itmax, niter=iter, $
               yfit=yfit, perror=sigma, covar=covar,/quiet)

  chi2 = total( ( (rhobyrhocrit-yfit)/error )^2, /double)
  dof = n_elements(rhobyrhocrit)-2
  chi2per = chi2/dof  

  IF NOT keyword_set(silent) THEN BEGIN 
      print
      print,'--------------------------------------------------------'
      print,'iterations = ',iter
      print,'r200  (a0) = ',ntostr(a[0])+' +/- '+ntostr(sigma[0])
      print,'Index (a1) = ',ntostr(a[1])+' +/- '+ntostr(sigma[1])
      print,'chi^2/dof  = ',ntostr(chi2)+'/'+ntostr(dof)+' = '+ntostr(chi2per)
  ENDIF 

  struct = {$
             rhobyrhocrit_fit: yfit, $
             r200: a[0], $
             r200_err: sigma[0], $
             r200_index: a[1], $
             r200_index_err: sigma[1], $
             chi2: chi2, $
             dof: dof, $
             chi2per: chi2per $
           }

  return, struct
END 
