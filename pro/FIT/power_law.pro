FUNCTION power_law, x, a

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
;
; NAME:
;    power_law
;       
; PURPOSE:
;    Fitted function to be called by lmfit.
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
;    returns the value of the function  F = a0*x^a1 and dF/da0, dF/da1
;	
; Use with lmfit:
;
; IDL> fn = 'power_law'
; IDL> a= [guess(a0), guess(a1)]
; IDL> weights = 1/err(y)^2
; IDL> yfit = lmfit(x, y, a, FUNCTION_NAME=FUNCTION_NAME, weights=weights, $
;                   sigma=sigma, chisq=chisq, convergence=convergence)
; REVISION HISTORY:
;    
;       
;                                      
;-                                       
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


  IF N_params() EQ 0 THEN BEGIN 
     print,'-Syntax: result = power_law(x, a)'
     print,''
     print,'Use doc_library,""  for more help.'  
     return,-1.
  ENDIF 

  ind1 = a[1]
  ind2 = a[1]-1.
  return, [ [ a[0]*x^ind1 ], $
            [      x^ind1 ], $
            [ a[0]*ind1*x^ind2] ]

END 
