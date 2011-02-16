PRO rosat_circ, matchstr, index, radius, clr, rerun=rerun,fieldzero=fieldzero, rosat=rosat

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
;
; NAME: rosat_circ
;       
; PURPOSE: circle object from matched structure
;	
;
; CALLING SEQUENCE: 
;      rosat_circ, matchstr, index, radius, clr,  rerun=rerun,
;                  fieldzero=fieldzero, rosat=rosat
;                 
;
; INPUTS: matchstr:  structure containing photo and rosat info for matched
;                    objects.
;         index:     index of matchstr
;
; INPUT KEYWORD PARAMETERS:
;         rerun:  the rerun number in integer form
;         fielszero: field zero in integer form
;         rosat: the rosat catalog
;
; CALLED ROUTINES:
; 
; PROCEDURE: 
;	            TSOBJ_NAME
;                   READ_PHOTO_COL
;                   FCHART_CIRC_RADEC
;	
;
; REVISION HISTORY:
;	Erin Scott Sheldon  Umich  5/25/99
;       
;                                      
;-                                       
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


  if N_params() eq 0 then begin
     print,'-Syntax: rosat_circ, matchstr, index, radius, clr, rerun=rerun,fieldzero=fieldzero, rosat=rosat'
     print,''
     print,'Use doc_library,""  for more help.'  
     return
  endif

  IF (n_elements(rosat) EQ 0) THEN BEGIN 
    rosatname = '/sdss3/data2/rosat/rass-bsc-1.2rxs.fit'
    rosat = mrdfits(rosatname, 1, hdr)
  ENDIF 

  ;;; use .4'' per pixel
  old_circ_rad = rosat[ matchstr[index].rosat_index ].pos_err
  circ_rad = old_circ_rad/.4
  print,'Rosat position error: ',strtrim(string(old_circ_rad),2),' arcsec  ',strtrim(string(circ_rad),2),' pixels'

  run = matchstr[index].run
  field = matchstr[index].field
  camcol=matchstr[index].camcol
  id = matchstr[index].id
  ra = matchstr[index].rosat_ra
  dec = matchstr[index].rosat_dec

  IF (n_elements(rerun) EQ 0) THEN rerun = 0
  IF (n_elements(fieldzero) EQ 0) THEN fieldzero = 11

  name=tsobj_name(run, camcol, fieldzero=fieldzero, rerun=rerun)
  dir = run_dir(run)
  filename = dir + name

  ;;; WARNING:  THIS ASSUMES OFFSET BETWEEN FIELD AND HDU IS fieldzero-1!
  ;;; READ 3 FIELDS AROUND THE OBJECT
  read_photo_col, filename, str, start = field - fieldzero, nframes=3

  w=where(str.field EQ field AND str.id EQ id)
  IF (w[0] EQ -1 ) THEN BEGIN
    print,'Something is terribly wrong'
    return
  ENDIF
  
  fchart_circ_radec, str, ra, dec, radius, clr, dir=dir, $
    image=image, photoid=w[0], /circ_obj, circ_rad=circ_rad, $
    objx=objx,objy=objy

;;; if you wanted to circle (or box) another thing on the image:
;  circ_radec, str, w[0], objx, objy, other_ra, other_dec, image,$
;                /nodisplay, color=other_color, radius=other_radius, box=box


return
end








