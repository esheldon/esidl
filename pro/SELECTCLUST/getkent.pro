PRO getkent, run, camcol, field, ra, dec, names=names, clusters=clusters

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
;
; NAME:
;    GETKENT
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
     print,'-Syntax: getkent, run, camcol, field, [, ra, dec, names=names]'
     print,' IMPORTANT: you need to start netscape before running this procedure'
     print,''
     print,'Use doc_library,""  for more help.'  
     return
  ENDIF 

  IF NOT keyword_set(clusters) THEN clusters = 0

  n=n_elements(run)

  FOR i=0, n-1 DO BEGIN 

      url =  getkenturl(run[i], camcol[i], field[i], clusters=clusters)

      IF n_elements(names) NE 0 THEN BEGIN
          print
          print,'Object Name:   ',+names[i]
      ENDIF 
      
      IF n_elements(ra) NE 0 THEN BEGIN
          IF n_elements(names) EQ 0 THEN print
          radecstr, ra[i], dec[i], rastr, decstr
          print,'RA: ',rastr,'  DEC: ',decstr
      ENDIF 

      spawn,"netscape -remote 'openURL("+url+")'",result

      IF (n GT 1) AND (i NE n-1) THEN BEGIN
          print,'Next object (y/n)?'
          key=get_kbrd(1)
          IF (key EQ 'n') OR (key EQ 'N') THEN return
      ENDIF 
  ENDFOR 

  return 
END 
