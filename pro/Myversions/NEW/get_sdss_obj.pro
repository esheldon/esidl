PRO get_sdss_obj, run, rerun, camcol, field, id, objstr, dir, atldir,all=all

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
;
; NAME:
;    GET_SDSS_OBJ
;       
; PURPOSE:
;    Get the photo structure for the specified object.
;
; CALLING SEQUENCE:
;    get_sdss_obj, run, rerun, camcol, field, id [, objstr, dir, $
;                  atldir, all=all]
;
; INPUTS: 
;    run,rerun,camcol,field,id
;
; OPTIONAL INPUTS:
;    
;
; KEYWORD PARAMETERS:
;    /all: get all objects from the field containing the input object
;       
; OUTPUTS: 
;    objstr: PHOTO structure of the object or the whole field if /all
;
; OPTIONAL OUTPUTS:
;    dir: directory containing the tsObj file
;    atldir: directory containing the atlas files 
;
; CALLED ROUTINES:
;    FETCH_DIR
;    READ_TSOBJ
; 
; PROCEDURE: 
;    
;	
;
; REVISION HISTORY:
;    Created ??-??-2001: Erin S. Sheldon UofMich
;       
;                                      
;-                                       
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


  IF N_params() EQ 0 THEN BEGIN 
     print,'-Syntax: get_sdss_obj, run, rerun, camcol, field, id [, objstr, dir, atldir, all=all]'
     print,''
     print,'Use doc_library,""  for more help.'  
     return
  ENDIF 

  fetch_dir, run, camcol, rerun, dir, atldir, /check
  IF dir[0] EQ '' THEN BEGIN 
      print,'No such run/rerun'
      delvarx, objstr
      return
  ENDIF 
  read_tsobj, dir, objstr, start = field
  w=where(objstr.id EQ id, nw)  

  IF nw NE 0 THEN BEGIN 
      IF NOT keyword_set(all) THEN BEGIN 
          objstr = objstr[w]
      ENDIF
  ENDIF ELSE BEGIN 
      print,'No such object in the field!'
      delvarx, objstr
  ENDELSE 
  
      
  return 
END 
