
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
;
; NAME:
;    htmLookupRadec
;       
; PURPOSE:
;    THIS PROGRAM IS OBSOLETE.  For IDL v5.5 or later, use id=htm_index()
;
;    find htm triangles associated with the input (ra,dec) vectors
;
; CALLING SEQUENCE:
;    htmLookupRadec, ra, dec, depth, leafid
;
; INPUTS: 
;    ra,dec: double precision positions, either scalar or array.
;
; OPTIONAL INPUTS:
;    NONE
;
;       
; OUTPUTS: 
;    leafid: the leaf id's for the (ra,dec) pairs.
;
; OPTIONAL OUTPUTS:
;    NONE
;
; CALLED ROUTINES:
;    DATATYPE
;    ISARRAY
;    htmLookupRadec.so (shared object library)
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

PRO htmLookupRadec, ra, dec, depth, leafid

  IF n_params() LT 3 THEN BEGIN 
      print,'THIS PROGRAM IS OBSOLETE.  For IDL v5.5 or later, use result=htm_index()'
      print,'-Syntax: htmLookupRadec, ra, dec, depth, leafid'
      print
      print,'Use doc_library,"htmLookupRadec"'
      return
  ENDIF 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; do some type checking
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  nra =  n_elements(ra)
  ndec = n_elements(dec)
  IF nra NE ndec THEN message,'ra and dec must be same size'

  IF (datatype(ra) NE 'DOU') OR (datatype(dec) NE 'DOU') THEN $
    message,'ra and dec must be of type double'
  IF NOT isarray(ra) THEN ra = [ra]
  IF NOT isarray(dec) THEN dec = [dec]
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; check for bad values
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  wbad = where(ra LT 0.,nbad)
  IF nbad NE 0 THEN message,'ra contains values < 0'
  wbad = where(ra GT 360., nbad)
  IF nbad NE 0 THEN message,'ra contains values > 360'

  wbad = where(dec LT -90.,nbad)
  IF nbad NE 0 THEN message,'dec contains values < -90'
  wbad = where(dec GT 90., nbad)
  IF nbad NE 0 THEN message,'dec contains values > 90'

  depth = ulong(depth)          ;convert to uint32
  n = ulong( n_elements(ra) )   ;;convert to uint32
  IF n_elements(leafid) NE 0 THEN delvarx,leafid
  leafid = replicate(ulong(0), n) ;leafid will be array of uint32

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; define shared object file and entry point
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  sofile = getenv('SDSSIDL_DIR')+'/src/lib/htmLookupRadec.so'
  entry = 'main'

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; call the shared object library, get leafid
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  time=systime(1)
  tmp = call_external(value=[0B,0B,0B,0B,0B], sofile,entry,$
                      ra, dec, depth, n, leafid)

END 
