FUNCTION adatc_dbid, RUN_OR_PHOTOID, camcol, field, id, rerun=rerun

  ;; Get the database ids that correspond to the input photo identification
  ;; This is very fast because the database is sorted by photoid

  np = n_params()
  IF np NE 1 AND np NE 4 THEN BEGIN 
      print,'-Syntax: list = adatc_dbid(photoids)'
      print,'              -- OR --'
      print,'         list = adatc_dbid(run,camcol,field,id, rerun=)'
      return,-1
  ENDIF 

  IF np EQ 1 THEN BEGIN 
      photoid_extract, RUN_OR_PHOTOID, run, rerun, camcol, field, id
  ENDIF ELSE BEGIN 
      run = RUN_OR_PHOTOID
  ENDELSE 

  IF n_elements(rerun) EQ 0 THEN BEGIN 
      rerun = sdss_rerun(run)
  ENDIF 

  adatc_dbopen, run[0], camcol[0], rerun=rerun[0], status=status
  IF status NE 0 THEN return,-1

  nobj = n_elements(run)

  IF np EQ 4 THEN BEGIN 
      IF n_elements(rerun) NE nobj THEN BEGIN 
          rerun = replicate(rerun, nobj)
      ENDIF 

      photoid = photoid(run,rerun,camcol,field,id)

      list = dbget('photoid', photoid)
  ENDIF ELSE BEGIN 
      ;; run is actually the photoid
      list = dbget('photoid', RUN_OR_PHOTOID)
  ENDELSE 

  dbclose

  return,list

END 
