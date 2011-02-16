
PRO runstatus_new, run_status

  ;; A replacement for good_reruns

  sdss_data_dir = sdssidl_config('data_dir')
  dir =  sdss_data_dir + "dbm/collate/"

  file = dir + "runstatus.dat"

  format = "L,L,L,L,A,L,L,L,L,L,L,F,F,F,F,F,A,A,A"

  readcol, file, $
    run, rerun, camcol, stripe, strip, $
    n_tsObj, n_tsField, n_fpAtlas, n_fpM, n_psField, n_adatc, $
    tsObj_photo_v, fpAtlas_photo_v, adatc_photo_v, baye_ver, phtz_ver, $
    imagingDir, adatcDir, host, $
    format=format

  ntot = n_elements(camcol)

  s=create_struct("stripe", -1, "strip", "?", "run", 0L, "rerun", -2, $
                  "tsobj_photo_v", -1.0, "fpatlas_photo_v", -1.0, $
                  "adatc_photo_v", -1.0, 'baye_ver', -1.0, 'phtz_ver', -1.0, $
                  "bad", 0L, "bad2", 0L)

  ;; Get unique identifire for run-rerun
  ten = ulong64(10)
  super = ulong64(run)*ten^4 + ulong64(rerun)

  rmd = rem_dup(super)
  nst = n_elements(rmd)

  run_status = replicate(s, nst)

  ;; bad flags
  tsObj_offset = 0
  fpAtlas_offset = 6
  psField_offset = 12
  fpM_offset = 18
  adatc_offset = 24
  
  ;; bad2 flags
  tsField_offset = 0
  photoz_offset = 6

  ;; Data must be sorted by run/rerun/camcol

  ;; i is the run/rerun index
  i = 0L
  FOR ii=0L, ntot-1 DO BEGIN 

      col = camcol[ii]

      IF col EQ 1 THEN BEGIN 
          ;; initialize the bitmasks
          bad = 0L
          bad2 = 0L
      ENDIF 

      IF n_tsObj[ii]   EQ 0 THEN bad  = bad  + 2L^(col+tsObj_offset)
      IF n_tsField[ii] EQ 0 THEN bad2 = bad2 + 2L^(col+tsField_offset)
      
      IF n_fpAtlas[ii] EQ 0 THEN bad  = bad  + 2L^(col+fpAtlas_offset)
      IF n_psField[ii] EQ 0 THEN bad  = bad  + 2L^(col+psField_offset)
      IF n_fpM[ii]     EQ 0 THEN bad  = bad  + 2L^(col+fpM_offset)

      IF n_adatc[ii]   EQ 0 THEN bad  = bad  + 2L^(col+adatc_offset)

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; Look for the photoz file
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      photoz_file = sdssidl_config('photoz_dir')+$
        'tsObj_ascii_'+ntostr(run[ii])+$
        '_'+ntostr(rerun[ii])+'_'+ntostr(col)+'.dat'
      IF NOT fexist(photoz_file) THEN bad2 = bad2 + 2L^(col+photoz_offset)

      ;; Only need to copy this in once. Should be camcol 6 since
      ;; we want to copy bitmasks
      IF col EQ 6 THEN BEGIN 

          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
          ;; Look for the asTrans file
          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

          transdir = SDSS_DATA_DIR+ $
            ntostr(run[ii])+'/'+$
            ntostr(rerun[ii])+'/astrom/'
          transfile = transdir + 'asTrans-'+run2string(run[ii])+'.fit'
          IF NOT fexist(transfile) THEN bad = bad + 2L^0

          run_status[i].stripe = stripe[ii]
          run_status[i].strip = strip[ii]

          run_status[i].run = run[ii]
          run_status[i].rerun = rerun[ii]

          run_status[i].tsObj_Photo_v = tsObj_Photo_v[ii]
          run_status[i].fpAtlas_Photo_v = fpAtlas_Photo_v[ii]
          run_status[i].adatc_Photo_v = adatc_Photo_v[ii]
          run_status[i].baye_ver = baye_ver[ii]
          run_status[i].phtz_ver = phtz_ver[ii]

          run_status[i].bad = bad
          run_status[i].bad2 = bad2
          i = i+1
      ENDIF 
      
  ENDFOR 
return

  run_status_file = sdssidl_config('RUN_STATUS_FILE')

  print
  print,'RUN_STATUS_FILE: ',run_status_file
  print
  mwrfits, run_status, run_status_file, /create

return
END 
