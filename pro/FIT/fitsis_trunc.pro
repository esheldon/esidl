PRO fitsis_trunc, x, y, error, aguess, yfit, a, asigma

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
;
; NAME:
;    
;       
; PURPOSE:
;    
;
; CALLING SEQUENCE:
;    
;
; INPUTS: 
;    
;
; OPTIONAL INPUTS:
;    
;
; KEYWORD PARAMETERS:
;    
;       
; OUTPUTS: 
;    
;
; OPTIONAL OUTPUTS:
;    
;
; CALLED ROUTINES:
;    
; 
; PROCEDURE: 
;    
;	
;
; REVISION HISTORY:
;    
;       
;                                      
;-                                       
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF N_params() EQ 0 THEN BEGIN 
     print,'-Syntax: fitpower, x, y, error, aguess [,yfit, a, asigma]'
     print,''
     print,'Use doc_library,""  for more help.'  
     return
  ENDIF 
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Parameters
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  asend = aguess                    ;because it will change.
  na = n_elements(asend)
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

  ;; a[0] = sigma_v km/s  convert to units of 1000km/s
  asend[0] = asend[0]/1000.
  ;; a[1] is the outer scale in units of kpc, convert to Mpc
  asend[1] = asend[1]/1000.
  ;; x is distance in Mpc.

  ;; y is in units of Msolar/pc^2  Convert to units of Msolar/Mpc^2
  ;; then convert to units of 1.15e14 Msolar/Mpc^2

  fac1 = 1.e12
  fac = fac1/1.15e14
  ysend = y*fac
  errorsend = error*fac

  weights = 1./errorsend^2
  
  print
  print,'asend = ',aguess
;  print,max(ysend),max(errorsend)

  itmax=1000
  tol = 1.e-5
  a = comfit2(x, ysend, asend, itmax=itmax, tol=tol, $
              /sigdiff, weights=weights, sigma=asigma, chisq=chisq, yfit=yfit, iter=iter)

  a[0] = a[0]*1000.
  a[1] = a[1]*1000.
  asigma[0] = asigma[0]*1000.
  asigma[1] = asigma[1]*1000.
  yfit = yfit/fac

  print,'      afit       err   '
  colprint,a,asigma
  print
  print,'chisq ',chisq



  return 
END 
