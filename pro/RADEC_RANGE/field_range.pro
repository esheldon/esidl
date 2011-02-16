;+
;
; NAME:
;    FIELD_RANGE
;       
; PURPOSE:
;    For a given run, create a fits file with lambda and
;    eta range info for each camcol/field.  It is a fits file with an array of
;    structures. Fields are not exactly rectangular in this system so there
;    may be error.  An alternative to find_radec.pro (which usees these
;    outputs) is find_radec_astrans.pro, which is more accurate but much
;    slower.
;
; CALLING SEQUENCE:
;    field_range, run, bad=
;
; INPUTS: 
;    run: The SDSS run
;
; OPTIONAL INPUTS:
;    None.
;
; KEYWORD PARAMETERS:
;    bad=: if returns 0 then ok, if returns 1 then bad
;       
; OUTPUTS: 
;    Writes out files to the directory "base" specified in the code. 
;    NOTE: This should be configured for your machine.
;
; OPTIONAL OUTPUTS:
;    None.
;
; CALLED ROUTINES:
;   sdss_read()
; 
; REVISION HISTORY:
;    Erin Scott Sheldon 10/27/99  UMich
;    Switched to survey coords 9-OCT-2000 E.S.S.
;    Switched to corrected survey coords 2007-05-21 E.S.S.
;       
;                                      
;-                                       

PRO field_range, run, bad=bad, astrans_dir=astrans_dir

  IF N_params() LT 1 THEN BEGIN 
     print,'-Syntax: field_range, run, bad=bad'
     print,''
     print,'Use doc_library,"field_range"  for more help.'  
     return
  ENDIF 

; Base directory in which to put the radec-range files.  

  base = sdssidl_config('radec_search_dir')

  bad=0

  colmin = 1
  colmax = 6

  typ = 'tmp'+ntostr(  long(systime(1)) )
  s=create_struct(name=typ, $
                  'field', 0L, $
                  'cetamin', 0d, $
                  'cetamax', 0d, $
                  'clambdamin', 0d, $
                  'clambdamax', 0d )  

; Find approximage ranges for each field

  clr=2                         ;read red astrans file
  mincol = 0
  maxcol = 2048-1

  file_base = 'clamceta-range-'+run2string(run)+'-'

  run_status = sdss_runstatus()
  w=where(run_status.run EQ run, nmatch)
  IF nmatch EQ 0 THEN message,'No such run in run_status struct: '+ntostr(run)

  maxrerun = long( max(run_status[w].rerun) )

  FOR camcol=colmin, colmax DO BEGIN 

      ;; output file
      addstr = ntostr(camcol)+'.fit'
      writefile = path_join(base,file_base + addstr)

      trans=sdss_read('astrans', run, camcol, band=clr, node=node, inc=inc, $
                      status=astatus, indir=astrans_dir, /silent)

      IF astatus ne 0 THEN BEGIN
          print,'Could not read astrans'
          bad=1
          return
      ENDIF 

      field = trans.field
      nfield = n_elements(field)

      struct = replicate(s, nfield)
      fieldmin = min(field)
      fieldmax = max(field)
      
      FOR fi=0L, nfield-1 DO BEGIN 
          ; don't include overlap region except first and last frames
          minrow = (field[fi] GT fieldmin)*64 + 1 
          maxrow = 1489 - (field[fi] LT fieldmax)*64 -1 
          ;; assume square
          rowcol2munu, trans, field[fi], minrow, mincol, mu1, nu1
          gc2csurvey, mu1, nu1, node, inc, lambda1, eta1
          
          rowcol2munu, trans, field[fi], maxrow, maxcol, mu2, nu2
          gc2csurvey, mu2, nu2, node, inc, lambda2, eta2

          rowcol2munu, trans, field[fi], minrow, maxcol, mu3, nu3
          gc2csurvey, mu2, nu2, node, inc, lambda3, eta3

          rowcol2munu, trans, field[fi], maxrow, mincol, mu4, nu4
          gc2csurvey, mu2, nu2, node, inc, lambda4, eta4

          etas = [eta1, eta2, eta3, eta4]
          lambdas = [lambda1, lambda2, lambda3, lambda4]
          cetamin = min(etas)
          cetamax = max(etas)
          clambdamin = min(lambdas)
          clambdamax = max(lambdas)

          struct[fi].field = field[fi]
          struct[fi].cetamin = cetamin
          struct[fi].cetamax = cetamax
          struct[fi].clambdamin = clambdamin
          struct[fi].clambdamax = clambdamax
      ENDFOR
      mwrfits, struct, writefile, /create
      delvarx, struct
  ENDFOR 

  return 
END 
