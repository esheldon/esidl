PRO rmclose_radec, ra, dec, indices, out, tol=tol

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
;
; NAME: 
;    rmclose_radec
;       
; PURPOSE: 
;    Remove duplicate objects from struct; that is, objects that have
;    the same position within some tolerance. Must be sorted by
;    ra and not cross RA=0 degree mark!!!!
;	
;
; CALLING SEQUENCE: 
;    rmclose_radec, ra, dec, indices, out           
;
; INPUTS: ra, dec
;
; OUTPUTS: indices.  The indices of non-degenerate entries.
;          out.      The stuff we threw out.
; 
;
; REVISION HISTORY:
;	Erin Scott Sheldon  UofMichigan  8/25/99
;       Modified version of Phil Fischer's code.
;       
;                                      
;-                                       
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  on_error,2

  IF N_params() LT 2 THEN BEGIN 
     print,'-Syntax: rmclose_radec, ra, dec, indices, out, tol=tol'
     print,''
     print,'Use doc_library,"rmclose_radec"  for more help.'  
     return
  ENDIF 

  ns = n_elements(ra)
  print
  print,'Throwing out duplicates'
  print,'Beginning with '+ntostr(ns)+' galaxies'
  print

  d2r=!dpi/180.
                                ;1 arcsecond default tolerance
  IF n_elements(tol) EQ 0 THEN tol1 = double(1.)/3600. ELSE tol1=tol
  tol2 = tol1*d2r

  ss = intarr(ns)
  FOR i=0L,ns-2L DO BEGIN 
      nn=i+1
      dra=abs( (ra[nn] - ra[i])/cos(dec[nn]) ) 

      WHILE (dra lt tol1) DO BEGIN

          gcirc, 0, ra[i]*d2r, dec[i]*d2r,  ra[nn]*d2r, dec[nn]*d2r, dist

          IF (dist LE tol2) THEN BEGIN
              IF (ss[i] NE -1) THEN ss[nn]=-1
          ENDIF 
          nn=nn+1
          IF (nn LT ns) THEN BEGIN 
              dra=abs((ra[nn]-ra[i])/cos(dec[nn]))

          ENDIF ELSE dra=1.e6
      ENDWHILE 

      IF ( (i MOD 10000) EQ 0 ) THEN print,'.',format='($,a)'

  ENDFOR 

  indices=where(ss NE -1)
  out = where(ss EQ -1)

return
END
