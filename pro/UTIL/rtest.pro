PRO rtest, v1, v2, r, s

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
     print,'-Syntax: rtest, v1, v2, r, s'
     print,' Returns r and s (significance)'
     print,''
     print,'Use doc_library,"rtest"  for more help.'  
     return
  ENDIF 

  n1 = n_elements(v1)
  n2 = n_elements(v2)
  IF n1 NE n2 THEN BEGIN
      print,'v1 and v2 must be same size'
      return
  ENDIF 

  m1 = mean(v1)
  m2 = mean(v2)

  d1 = sqrt( total( (v1 - m1)^2 ) )
  d2 = sqrt( total( (v2 - m2)^2 ) )
  
  r = total( (v1 - m1)*(v2 - m2) )/d1/d2

  s = 1. - errorf( abs(r)*sqrt(n1/2.) )

  return 
END 
