FUNCTION gammfac, beta, inverse=inverse

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


  IF N_params() LT 1 THEN BEGIN 
     print,'-Syntax: result = gammfac(beta, inverse=inverse)'
     print,''
     print,'Use doc_library,"gammfac"  for more help.'  
     return,0.
  ENDIF 
  
  
  
  IF NOT keyword_set(inverse) THEN BEGIN
      return, 1./sqrt(1.-beta^2)
  ENDIF ELSE BEGIN
      return, sqrt(1. - 1./beta^2)
  ENDELSE 

END 
