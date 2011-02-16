pro rosat_extract,run,fieldzero=fieldzero,rerun=rerun,addon=addon,groupn=groupn,taglist=taglist, startcol=startcol, endcol=endcol, command_name=command_name

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
;
; NAME: rosat_extract
;       
; PURPOSE: create catalogues of matches between rosat and photo.
;	
;
; CALLING SEQUENCE: rosat_extract, run, fieldzero=fieldzero, rerun=rerun,
;     addon=addon, groupn=groupn, taglist=taglist, startcol=startcol, $
;     endcol=endcol, command_name=command_name
;      
;                 
;
; INPUTS: run:  the run number in integer form
;
; INPUT KEYWORD PARAMETERS:
;         fieldzero:  field zero in integer form
;         rerun:      rerun number in integer form
;         addon:      file name formatting
;         groupn:     number read in by read_photo_col at a time.
;         taglist:    the taglist for read_photo_col
;         startcol:   which column to do first
;         command_name:  string containing the command of the matching routine
;         
; OUTPUTS: reformats the match_key file.  The programs it calls have outputs.
;     in particular, sdss_extract ouputs the matched files for each column
;     and combine_rosat combines the photo info with these files.
;
; CALLED ROUTINES:
;                  RUN2STRING
;                  FIELD2STRING
;                  TSOBJ_NAME
;                  SDSS_EXTRACT
;                  COMBINE_ROSAT
;                  
; PROCEDURE: 
;	
;	
;
; REVISION HISTORY:
;	Erin Scott Sheldon  Umich  5/25/99
;       
;                                      
;-                                       
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;



  if n_params() eq 0 then begin
    print,'rosat_extract,run, fieldzero=fieldzero, rerun=rerun, addon=addon, groupn=groupn, taglist=taglist, startcol=startcol, endcol=endcol, command_name=command_name'
    return
  endif


  IF (n_elements(groupn) EQ 0) then groupn = 50
  IF (n_elements(fieldzero) EQ 0)  then fieldzero = 11
  IF (n_elements(rerun) EQ 0) THEN rerun=0
  IF (n_elements(addon) EQ 0) THEN addon = '-rosat_match.fit'
  IF (n_elements(startcol) EQ 0) THEN startcol=1
  IF (n_elements(endcol) EQ 0) THEN endcol=6
  IF (n_elements(command_name) EQ 0) THEN command_name = 'match_rosat'

  dir = run_dir(run)
  IF (dir EQ '') THEN BEGIN
    print,'Bad run number'
    return
  ENDIF 

  runstr = run2string(run)
  fieldzerostr = field2string(fieldzero)
  rerunstr = strtrim(string(rerun), 2)

  for i=startcol, endcol do begin

    print,''
    print,'Processing column:',i
    fname=tsobj_name(run, i, fieldzero=fieldzero, rerun=rerun,/full)
    oname = dir + 'tsObj-'+ runstr +'-'+$
      strtrim(string(i),2)+'-'+rerunstr+'-'+fieldzerostr+addon

    sdss_extract,fname,oname,command_name,groupn=groupn,taglist=taglist

  endfor


;;; put match key into one big file
  name = dir + 'tsObj-'+runstr+'-*-0-0011-rosat_match_key.fit'
  fits_info, name , /silent, n_ext=n_ext

  FOR i=1, n_ext DO BEGIN
    s = mrdfits(name, i, hdr, structyp='match_key')
    IF (i EQ 1) THEN cs = s ELSE cs = [cs,s]
  ENDFOR

  mwrfits, cs, name, /create

;;; now combine all the rosat stuff

  combine_rosat, run, fieldzero=fieldzero, rerun=rerun

  return
end









