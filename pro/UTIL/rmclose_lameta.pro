PRO rmclose_lameta, lambda, eta, indices, out, tol=tol

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
     print,'-Syntax: rmclose_radec, lambda, eta, indices, out, tol=tol'
     print,''
     print,'WARNING: struct must be sorted by lambda'
     print,'Use doc_library,"rmclose_lameta"  for more help.'  
     return
  ENDIF 

  ns = n_elements(lambda)
  print
  print,'Throwing out duplicates'
  print,'Beginning with '+ntostr(ns)+' galaxies'
  print


  d2r=!dpi/180.
                                ;1 arcsecond default tolerance
  IF n_elements(tol) EQ 0 THEN tol1 = double(1.)/3600. ELSE tol1=tol
  tol2=tol1*d2r

  ss = intarr(ns)
  FOR i=0L,ns-2L DO BEGIN 
      nn=i+1
      
      dlam = abs( lambda[nn] - lambda[i] )
;      dlam = abs( struct[nn].(lind) - struct[i].(lind) )
      WHILE (dlam lt tol1) DO BEGIN
;          gcirc, 0, struct[i].(eind)*d2r, struct[i].(lind)*d2r, $
;            struct[nn].(eind)*d2r, struct[nn].(lind)*d2r, dist
          gcirc, 0, eta[i]*d2r, lambda[i]*d2r, $
                 eta[nn]*d2r, lambda[nn]*d2r, dist
          IF (dist LE tol2) THEN BEGIN
              IF (ss[i] NE -1) THEN ss[nn]=-1
          ENDIF 
          nn=nn+1
          IF (nn LT ns) THEN BEGIN 
;              dlam = abs( struct[nn].(lind) - struct[i].(lind) )
              dlam = abs( lambda[nn] - lambda[i] )
          ENDIF ELSE dlam=1.e6
      ENDWHILE 
  ENDFOR 

  indices=where(ss NE -1)
  out = where(ss EQ -1)

return
END
