PRO mextract, k, ksig, kcat, n, nsig, ncat, step=step

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
;
; NAME:
;    MEXTRACT
;       
; PURPOSE:
;    Extract peaks (both negative and positive) from a 2-d image.  Only 
;    works well for smoothed images.
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


  IF N_params() LT 3 THEN BEGIN 
     print,'-Syntax: mextract, k, ksig, kcat, n, nsig, ncat, step=step'
     print,''
     print,'Use doc_library,"mextract"  for more help.'  
     return
  ENDIF 

  IF n_elements(n) NE 0 THEN DO_n = 1 ELSE DO_n = 0
  IF DO_n THEN print,'Will do noise'
  t=systime(1)

  IF keyword_set(step) THEN BEGIN 

      sz = size(k)
      sx = sz[1]
      sy = sz[2]

      ly = lindgen(sy)

      stepfac = step*sx
      nstep = sy/stepfac
      print,'Using ',ntostr(nstep),' steps'

      FOR i=0, nstep-1 DO BEGIN
          
          yr = ly[i*stepfac:(i+1)*stepfac - 1]

          find_peak, k[*, yr], tkcat
          tkcat.s2n = ksig[ tkcat.x, yr[ tkcat.y ] ]
          tkcat.wild = i
          IF i EQ 0  THEN kcat = tkcat ELSE kcat = [kcat, tkcat]

          IF DO_n THEN BEGIN 
              find_peak, n[*, yr], tncat
              tncat.s2n = nsig[ tncat.x, yr[ tncat.y ] ]
              tncat.wild = i
              IF i EQ 0 THEN ncat = tncat ELSE ncat = [ncat, tncat]
          ENDIF 

      ENDFOR 

  ENDIF ELSE BEGIN 

      print,'Finding Peaks'
      find_peak, k, kcat
      print,'Finding Peak S/N'
      kcat.s2n = ksig[ kcat.x, kcat.y ]
      kcat.wild = 0

      IF DO_n THEN BEGIN 
          print,'Finding Noise Peaks'
          find_peak, n, ncat
          print,'Finding Noise Peak S/N'
          ncat.s2n = nsig[ ncat.x, ncat.y ]
          ncat.wild = 0
      ENDIF 

  ENDELSE 

  ptime,systime(1) - t

return
END 

