PRO nhist, array, binsize, xhist, nyhist, fac=fac, overplot=overplot, psym=psym, _extra=extra, noplot=noplot

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


  IF N_params() LT 2 THEN BEGIN 
     print,'-Syntax: nhist, array, binsize, xhist, nyhist, overplot=overplot, _psym=psym, extra=extra'
     print,''
     print,'Use doc_library,""  for more help.'  
     return
  ENDIF 



  IF n_elements(psym) EQ 0 THEN psym = 10

  plothist, array, xhist, yhist, bin=binsize, /noplot

  nyhist = float(yhist)/n_elements(array)

  IF NOT keyword_set(noplot) THEN BEGIN 
      IF keyword_set(overplot) THEN BEGIN 
          oplot, xhist, nyhist, psym=psym, _extra=extra
      ENDIF ELSE BEGIN 
          plot, xhist, nyhist, psym=psym, _extra=extra
      ENDELSE 
  ENDIF 
  
  return 
END 
