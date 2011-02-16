pro  combine_rosat, run, fieldzero=fieldzero, rerun=rerun

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
;
; NAME: combine_rosat
;       
; PURPOSE: combine info from rosat matched photo structures with rosat data
;	
;
; CALLING SEQUENCE: combine_rosat, run, fieldzero=fieldzero, rerun=rerun
;      
;
; INPUTS: run: run in integer form
;
; INPUT KEYWORD PARAMETERS:
;         fieldzero:  field zero in integer form
;         rurun:      rerun number in integer form
;       
; OUTPUTS: combined structures (to files)
;
;
; CALLED ROUTINES:
;                   RUN2STRING
;                   RUN_DIR
;                   COMBINE_STRUCTS
; 
; PROCEDURE: 
;	
;	
;
; REVISION HISTORY:
;	Erin Scott Sheldon Umich 5/25/99
;       
;                                      
;-                                       
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


  if N_params() eq 0 then begin
     print,'-Syntax: combine_rosat, run, fieldzero=fieldzero'
     print,''
     print,'Use doc_library,""  for more help.'  
     return
  endif


dir = run_dir(run)
IF (dir EQ '') THEN BEGIN
  print,'Bad run number'
  return
ENDIF 

runstr = run2string(run)
IF n_elements(fieldzero) eq 0 THEN fieldzerostr = '0011' ELSE $
                fieldzerostr=field2string(fieldzero)
IF n_elements(rerun) EQ 0 THEN rerunstr = '0' ELSE $
                               rerunstr = strtrim(string(rerun),2)
keyname=dir + 'tsObj-'+runstr+'-*-'+'0-0011-rosat_match_key.fit'
key = mrdfits(keyname, 1, hdr)

rosatname = '/sdss3/data2/rosat/rass-bsc-1.2rxs.fit'
rosat = mrdfits(rosatname, 1, hdr)

r = create_struct('rosat_index', 0L)
rr = replicate(r, n_elements(key) )
rr.rosat_index = key.rosat_index

FOR i=1,6 DO BEGIN
  name = dir + 'tsObj-'+runstr+'-'+strtrim(string(i),2) + $
           '-'+rerunstr+'-'+fieldzerostr+'-rosat_match.fit'
  cat = mrdfits(name,1,hdr)
  ;;; note everything should already be in the right order
  w = where( key.photo_camcol EQ i)
  ;;; this code is found in daves stuff
  combine_structs, cat, rosat[ key[w].rosat_index ], outcat_tmp
  combine_structs, outcat_tmp, rr[w], outcat

  ;;; now sort by the rosat index
  s=sort(outcat.rosat_index)
  outcat = outcat[s]
  mwrfits, outcat, name,/create
ENDFOR 


return
end





