PRO rmclose_lameta_seeing, struct, indices, out, tol=tol, $
                           seeing_index=seeing_index

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
;
; NAME: 
;    rmclose_lameta
;       
; PURPOSE: 
;    Remove duplicate objects from struct; that is, objects that have
;    the same position within some tolerance.  Must be sorted by
;    lambda.
;	
;
; CALLING SEQUENCE: 
;    rmclose_radec, struct, indices, out           
;
; INPUTS: struct   Structure that must contain ra and dec as tags.
;
; OUTPUTS: indices.  The indices of non-degenerate entries.
;          out.      The stuff we threw out.
; 
;
; REVISION HISTORY:
;       
;                                      
;-                                       
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


  IF N_params() EQ 0 THEN BEGIN 
     print,'-Syntax: rmclose_radec_seeing, struct, indices, out, tol=tol, seeing_index=seeing_index'
     print,''
     print,'Use doc_library,"rmclose_radec"  for more help.'  
     return
  ENDIF 

  IF n_elements(seeing_index) EQ 0 THEN sid = 0 ELSE sid = seeing_index

  ns = n_elements(struct)
  print
  print,'Throwing out duplicates (keeping best-seeing)'
  print,'Beginning with '+ntostr(ns)+' galaxies'
  print


  IF NOT tag_exist(struct, 'clambda', index=lind) THEN BEGIN
      IF NOT tag_exist(struct, 'lambda', index=lind) THEN BEGIN
          message,'struct must contain "clambda" or "lambda" tags'
      ENDIF 
  ENDIF 

  IF NOT tag_exist(struct, 'ceta', index=eind) THEN BEGIN
      IF NOT tag_exist(struct, 'eta', index=eind) THEN BEGIN
          message,'struct must contain "ceta" or "eta" tags'
      ENDIF 
  ENDIF 

  d2r=!dpi/180.
                                ;1 arcsecond default tolerance
  IF n_elements(tol) EQ 0 THEN tol1 = double(1.)/3600. ELSE tol1=tol
  tol2=tol1*d2r

  ss = intarr(ns)
  FOR i=0L,ns-2L DO BEGIN 
      nn=i+1
      
      dlam = abs( struct[nn].(lind) - struct[i].(lind) )
      WHILE (dlam lt tol1) DO BEGIN
          gcirc, 0, struct[i].(eind)*d2r, struct[i].(lind)*d2r, $
            struct[nn].(eind)*d2r, struct[nn].(lind)*d2r, dist
          IF (dist LE tol2) THEN BEGIN

              ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
              ;; now keep the one that had the smallest seeing
              ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

              CASE 1 OF
                  (ss[i] NE -1) AND (ss[nn] NE -1): BEGIN
                      IF struct[i].seeing[sid] LT struct[nn].seeing[sid] THEN ss[nn] = -1 $
                      ELSE ss[i] = -1
                  END 
                  ELSE: IF (ss[i] NE -1) THEN ss[nn]=-1
              ENDCASE 

          ENDIF 
          nn=nn+1
          IF (nn LT ns) THEN BEGIN 
              dlam = abs( struct[nn].(lind) - struct[i].(lind) )
          ENDIF ELSE dlam=1.e6
      ENDWHILE 
  ENDFOR 

  indices=where(ss NE -1)
  out = where(ss EQ -1)

return
END
