;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
;
; NAME:
;    htmIntersectRadec
;       
; PURPOSE:
;    THIS PROGRAM IS OBSOLETE.  For IDL v5.5 or later, use lf=htm_intersect()
;
;    Find htm triangles that are within angle of input (ra,dec) position.
;    Some triangles may only partially intersect.
;
; CALLING SEQUENCE:
;    htmIntersectRadec, ra, dec, angle, depth, leaflist
;
; INPUTS: 
;    ra,dec: double precision position, must be scalar.
;
; OPTIONAL INPUTS:
;    NONE
;
;       
; OUTPUTS: 
;    leaflist: the leaf id's for the (ra,dec) pairs.
;
; OPTIONAL OUTPUTS:
;    NONE
;
; CALLED ROUTINES:
;    DATATYPE
;    ISARRAY
;    htmIntersectRadec.so (shared object library)
;
; PROCEDURE: 
;    
;	
;
; REVISION HISTORY:
;    ??-NOV-2000 Erin Scott Sheldon UofMich
;       
;                                      
;-                                       
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


PRO htmIntersectRadec, ra, dec, angle, depth, leaflist, numdefault=numdefault

  IF n_params() LT 4 THEN BEGIN 
      print,'THIS PROGRAM IS OBSOLETE.  For IDL v5.5 or later, use result=htm_intersect()'
      print,'-Syntax: htmIntersectRadec, ra, dec, angle, depth, leaflist, numdefault=numdefault'
      print
      print,' ra,dec,angle must be double precision'
      print,' ra,dec in degrees'
      print,' angle in radians'
      print,'Use doc_library,"htmIntersectRadec"'
      return
  ENDIF 

  time=systime(1)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; do some type checking
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF (datatype(ra) NE 'DOU') OR $
     (datatype(dec) NE 'DOU') OR $
     (datatype(angle) NE 'DOU')THEN $
    message,'ra,dec,angle must be of type double'

  IF isarray(ra) OR isarray(dec) OR isarray(angle) THEN $
    message,'ra,dec,angle must be scalars'

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; check for bad values
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF (ra LT 0.) THEN message,'ra < 0'
  IF (ra GT 360.) THEN message,'ra > 360.'
  IF (dec LT -90.) THEN message,'dec < -90.'
  IF (dec GT 90.) THEN message,'dec > 90.'

  depth = long(depth)           ;Make sure long for C program
  d = cos( angle )              ;C-program takes in cos(angle)
  IF n_elements(numdefault) EQ 0 THEN BEGIN 
      numdefault = 10000L        ;default # of values in leaflist
  ENDIF 

  leaflist = replicate( ulong(0), numdefault )

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; define shared object file and entry point
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  sofile = getenv('SDSSIDL_DIR')+'/src/lib/htmIntersectRadec.so'
  entry = 'main'

  tmp = call_external(value=[0B,0B,0B,0B,0B], sofile, entry,$
                      depth, ra, dec, d, leaflist)

  w=where(leaflist NE 0,nw)
  IF nw NE 0 THEN leaflist=leaflist[w] ELSE leaflist = -1


END 
