PRO rmclose_radec_seeing_test, struct, indices, out, tol=tol

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
;    rmclose_radec, struct, indices, out           
;
; INPUTS: struct   Structure that must contain ra and dec as tags.
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


  IF N_params() EQ 0 THEN BEGIN 
     print,'-Syntax: rmclose_radec, struct, indices, out, tol=tol'
     print,''
     print,'Use doc_library,"rmclose_radec"  for more help.'  
     return
  ENDIF 

  d2r=!dpi/180.
                                ;1 arcsecond default tolerance
  IF n_elements(tol) EQ 0 THEN tol1 = double(1.)/3600. ELSE tol1=tol
  tol2 = tol1*d2r

  ns = n_elements(struct)
  ss = intarr(ns)
  FOR i=0L,ns-2L DO BEGIN 
      nn=i+1
      dra=abs( (struct[nn].ra-struct[i].ra)/cos(struct[nn].dec) )
      WHILE (dra lt tol1) DO BEGIN
          mygcirc, struct[i].ra, struct[i].dec, struct[nn].ra, struct[nn].dec, $
            dist, /radians_out
          IF (dist LE tol2) THEN BEGIN
              CASE 1 OF
                  (ss[i] NE -1) AND (ss[nn] NE -1): BEGIN
                      IF struct[i].seeing LT struct[nn].seeing THEN ss[nn] = -1 $
                      ELSE ss[i] = -1
                  END 
                  ELSE: IF (ss[i] NE -1) THEN ss[nn]=-1
              ENDCASE 
          ENDIF 
          nn=nn+1
          IF (nn LT ns) THEN BEGIN 
              dra=abs((struct[nn].ra-struct[i].ra)/cos(struct[nn].dec))
          ENDIF ELSE dra=1.e6
      ENDWHILE 
  ENDFOR 

  indices=where(ss NE -1)
  out = where(ss EQ -1)

return
END
