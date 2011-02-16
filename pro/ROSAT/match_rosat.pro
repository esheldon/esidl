pro  match_rosat, photostr, photo_match, radius=radius, rosat=rosat,allow=allow

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
;
; NAME:  match_rosat
;       
; PURPOSE: match rosat catalog to a photo structure.
;	
;
; CALLING SEQUENCE: match_rosat, photostr, photo_match, radius=radius, $
;      rosat=rosat,allow=allow
;      
;                 
;
; INPUTS: photostr: photo structure.
;
; INPUT KEYWORD PARAMETERS:
;         radius: radius in arcseconds to match with. (default 20)
;         rosat:  the rosat catalog.
;         allow:  the number of allowed matches (default 5)
;       
; OUTPUTS: photo_match
;
; CALLED ROUTINES:
;                    CLOSE_MATCH
; 
;
; REVISION HISTORY:
;	Erin Scott Sheldon  Umich 5/25/99
;       
;                                      
;-                                       
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


  if N_params() eq 0 then begin
     print,'-Syntax: match_rosat, photostr, photo_match, radius=radius, rosat=rosat, allow=allow'
     print,''
     print,'Use doc_library,"match_rosat"  for more help.'  
     return
  endif

  IF (n_elements(rosat) EQ 0) THEN BEGIN 
    rosat_file = '/sdss3/data2/rosat/rass-bsc-1.2rxs.fit'
    rosat = mrdfits(rosat_file, 1, hdr)
  ENDIF 

;;;;; Find the indices of the phstruct that match the rosat struct
;;;;; in ra and dec.  Since there is a nominal 6" systematic error in the
;;;;; rosat positions, we will make a close match to be something bigger
;;;;; that 6".  i.e. use factor >= 1.0


  IF n_elements(radius) EQ 0 THEN radius = 20.0  ; radius in arcsec
  print,'Using ',strtrim(string(radius),2),' arcsecond tolerance'
  tolerance = radius*2.8e-4

;;;; Remember that rosat is the known catalog, so we input it as first
;;;; arrays in close_match.  This lets us return many matches from 
;;;; photostr.  Set allow to an appropriate value.

  IF (NOT keyword_set(allow)) THEN allow = 5
  close_match, rosat.rosat_ra, rosat.rosat_dec, photostr.ra, photostr.dec, $
             rosat_match, photo_match, tolerance, allow
 
  run = photostr[0].run
  runstr = run2string(run)
  dir = run_dir(run)
  IF (rosat_match[0] NE -1) THEN BEGIN
    nn = n_elements(rosat_match)
    ss = create_struct('photo_id', 0L, 'photo_field',0L, $
                       'photo_camcol', 0L, 'rosat_index', 0L, $
                       'rosat_ra', 0.0, 'rosat_dec', 0.0 )
    matches = replicate(ss, nn )
    FOR i=0, nn-1 DO BEGIN
      matches[i].photo_id = photostr[ photo_match[i] ].id
      matches[i].photo_field = photostr[ photo_match[i] ].field
      matches[i].photo_camcol = photostr[ photo_match[i] ].camcol
      matches[i].rosat_index = rosat_match[i]
    ENDFOR
    ;;;; this will create or append the file
    mwrfits, matches, dir+'tsObj-'+runstr+'-*-0-0011-rosat_match_key.fit'
  ENDIF
  

return
end













