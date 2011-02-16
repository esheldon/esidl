
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
;
; NAME:
;    GET_GOOD_LENSTRIPE
;       
; PURPOSE:
;    Get all the runs ready for lensing in a particular stripe, 
;    those with tsObj files, 
;    fpAtlas files, an asTrans file, and psField files. Since psField
;    files can be used for more than one rerun, we also return the
;    variable psfieldrerun, which can point at a rerun that does have
;    the psfield files.  This is a wrapper for GET_GOOD_LENSRUNS which
;    also checks the stripe number
;
;
; CALLING SEQUENCE:
;    get_good_lenstripe, stripe, runs, reruns, stripes, strips, 
;          indices, /silent
;
; INPUTS: 
;    stripe: stripe number in integer form
;
; KEYWORD PARAMETERS:
;    /silent:
;       
;
; OUTPUTS: 
;    runs,reruns,psfieldreruns,stripes,strips,indices (in !run_status struct)
;
; OPTIONAL OUTPUTS:
;    none
;
; CALLED ROUTINES:
;    GET_GOOD_LENSRUNS
;
; PROCEDURE: 
;    call GET_GOOD_LENSRUNS and then choose those with the input stripe
;	
;
; REVISION HISTORY:
;    Created: Feb-03-2003  Erin Scott Sheldon UChicago
;                            
;-                                       
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PRO get_good_lenstripe, stripe, runs, reruns, stripes, strips, indices, ngood,$
                        silent=silent

  delvarx,runs,reruns,stripes,strips

  IF N_params() LT 1 THEN BEGIN 
     print,'-Syntax: get_good_lenstripe, stripe, runs, reruns, stripes, strips, indices, ngood, /silent'
     print,''
     print,'Use doc_library,"get_good_lenstripe"  for more help.'  
     return
  ENDIF 

  rs = sdss_runstatus()
  wstripe=where(rs.stripe eq stripe, nstripe)
  if nstripe eq 0 then begin
      message,'No runs found for stripe: '+ntostr(stripe)
  endif


  fs = {astrans_exist:'y',tsobj_exist:'y'}
  good=sdss_flag_select(rs.flags, 'runstatus', fs, ngood, input_index=wstripe)

  IF ngood eq 0 THEN BEGIN 
      print
      print,'No good runs/reruns found'
      return
  ENDIF

  runs = rs[good].run
  reruns = rs[good].rerun
  stripes = rs[good].stripe
  strips = rs[good].strip

  IF NOT keyword_set(silent) THEN BEGIN 
      print
      print,'        runs  reruns strips'
      colprint,runs,reruns,'         '+strips
  ENDIF 

  return

END 
