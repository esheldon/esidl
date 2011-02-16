; This is specifically designed for the area fraction stuff
; for maxbcg galaxy counts

FUNCTION fitareapoly, x, y, error, degree=degree, $
                      yfit=yfit, sigma=sigma, $
                      silent=silent


  IF N_params() LT 3 THEN BEGIN 
      on_error, 2
      print,'-Syntax: parms = FitAreaPoly(x, y, error, yfit=, sigma=, silent=)]'
      print,''
      message,'Halting'
  ENDIF 
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Parameters
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF n_elements(degree) EQ 0 THEN degree = 5

  ;; First we do a poly fit on the data - 1.0 and use these
  ;; as our guess for the fit with constant = 1.0

  aguess = poly_fit(x, y, degree, measure_error=error,yfit=tyfit)
  aguess = reform(aguess)
  aguess = aguess[1:degree]

  nx = n_elements(x)
  ny = n_elements(y)
  IF ny NE nx THEN BEGIN
      print,'X and Y must be same size'
      on_error, 2
      message,'Halting'
  ENDIF 

  w=where(error EQ 0., nw)
  IF nw NE 0 THEN BEGIN
      print,'Errors must be non-zero'
      on_error, 2
      message,'Halting'
  ENDIF 

  itmax = 300
  MYFUNCT = 'areapoly'
  a = MPFITFUN(MYFUNCT, x, y, error, aguess, $
               AUTODERIVATIVE=1, MAXITER=itmax, niter=iter, $
               yfit=yfit, perror=sigma, covar=covar,/quiet,/double)

  IF NOT keyword_set(silent) THEN BEGIN 

      chi2 = total( ( (y-yfit)/error )^2, /double)
      dof = n_elements(y)-(degree-1)
      chi2per = chi2/dof
      print
      print,'--------------------------------------------------------'
      print,'best fit: '
      FOR i=0,degree-1 DO BEGIN 
          print,'  a['+ntostr(i)+'] = '+ntostr(a[i])+' '+!plusminus+' '+ntostr(sigma[i])
      ENDFOR 
      print,'iterations = ',iter
      print,'chi^2/dof  = ',ntostr(chi2)+'/'+ntostr(dof)+' = '+ntostr(chi2per)
  ENDIF 



  return, a
END 
