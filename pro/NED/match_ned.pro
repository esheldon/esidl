pro  match_ned, pstruct, type, tol, ph_id, ned_id, $
                ned_cat=ned_cat, nedfile=nedfile, $
                plot=plot, atlas=atlas, dir=dir, $
                silent=silent

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
;
; NAME:   
;    MATCH_NED
;       
; PURPOSE:  
;    Match the ned catalog with a photo struct 
;	
;
; CALLING SEQUENCE:match_ned, pstruct, type, tol, ph_id, ned_cat=ned_cat,
;        ned_id=ned_id, plot=plot, atlas=atlas,dirdir=dirdir
;      
;                 
;
; INPUTS: pstruct: photo struct 
;         type: what ned type to use
;         tol: tolerance for matching by ra and dec
;       
; OPTIONAL INPUTS:
;         dir: directory to look for atlas images
;         nedfile: the fits file holding the ned catalog.
;         WARNING: the default value will have to be changed on other
;                  machines.
;
; KEYWORD PARAMETERS:
;         /plot: make color-color plots and size plot
;         /atlas: view atlas images
;         /silent: Don't print messages except errors
;
; OUTPUTS: ph_id the id's of the photo struct.
;
; OPTIONAL OUTPUTS: ned_id: the id's of the ned structure
; 
; PROCEDURE: 
;	
;	
;
; REVISION HISTORY:  
;    Author: Erin Scott Sheldon UM  2/9/98
;	
;       
;                                      
;-                                       
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


  IF N_params() LT 3 THEN BEGIN 
	print,'Syntax: match_ned, pstruct, type, tol [, ph_id, ned_id, ned_cat=ned_cat, nedfile=nedfile, plot=plot, atlas=atlas, dir=dir, silent=silent]'
	print,'1 arcsecond = 2.8e-4 degrees'
        print
        print,'Use doc_library,"match_ned" for more help'
	return
  ENDIF
  
  IF NOT keyword_set(silent) THEN silent = 0
  IF NOT keyword_set(atlas) THEN atlas = 0
  IF NOT keyword_set(plot) THEN plot=0
  IF NOT keyword_set(dir) THEN dir = ''
  
  
  IF NOT keyword_set(ned_cat) THEN BEGIN 
    IF NOT keyword_set(nedfile) THEN BEGIN 
      nedfile = '/sdss3/data2/Ned/sdss_equatorial_with_z.fit'
    ENDIF 
    get_ned, ned_cat, file=nedfile
  ENDIF 

  ned_types = ned_cat[ rem_dup(ned_cat.type) ].type

  wtype = where( ned_types EQ type, nwt)
  IF nwt EQ 0 THEN BEGIN
    print,'Type "',type,'" not found in ned'
    print,'Try one of these:'
    print,ned_types
    return
  ENDIF 


  w = where(ned_cat.type EQ type, nw)
  IF NOT silent THEN BEGIN
    print,ntostr(nw),' ',type,' found in ned'
    print,'Matching by ra and dec'
  ENDIF 
  
  allow=1
  close_match_radec, ned_cat[w].ra,ned_cat[w].dec, $
                     pstruct.ra, pstruct.dec, $
                     ned_id, ph_id, $
                     tol, allow

  IF ph_id[0] EQ -1 THEN return

  IF plot THEN BEGIN
    extract_stars, pstruct, 2, stars, silent=1
    plot_colors, stars, oplot_str=pstruct[ph_id], /size
    IF atlas THEN BEGIN 
      print,'Hit space to see atlas images.'
      key=get_kbrd(1)
    ENDIF
  ENDIF  
  IF atlas THEN get_atlas,pstruct,ph_id,dir=dir

  return
END 



