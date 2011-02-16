PRO good_corr_reruns

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
;
; NAME:
;    GOOD_CORR_RERUNS
;       
; PURPOSE:
;    Find the rerun for runs in the !SDSS_SHAPECORR_DIR, and see if there
;    are any adatc files there.
;
; CALLING SEQUENCE:
;    good_corr_reruns
;
; INPUTS: 
;    NONE: gets its info from system variables defined in SDSSIDL_SETUP
;
; OPTIONAL INPUTS:
;    NONE:
;
; KEYWORD PARAMETERS:
;    NONE
;       
; OUTPUTS: 
;    Outputs !RUNCORR_STATUS_FILE, a fits file containing info about
;    known runs.
;
; OPTIONAL OUTPUTS:
;    NONE
;
; CALLED ROUTINES:
;    SDSSIDL_SETUP
;    ADD_ARRVAL
;    FETCH_RERUN
;    FETCH_DIR
;    FETCH_FILE_LIST
;    MWRFITS
;
; PROCEDURE: 
;    
;	
;
; REVISION HISTORY:
;    14-NOV-2000 Erin Scott Sheldon UofMich
;       
;                                      
;-                                       
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  sdssidl_setup
  IF NOT !sdssidl_def.sdss_shapecorr_dir_defined THEN BEGIN 
      message,'!SDSS_SHAPECORR_DIR must be defined'
  ENDIF 
  IF NOT !sdssidl_def.runcorr_status_file_defined THEN BEGIN 
      message,'!RUNCORR_STATUS_FILE must be defined'
  ENDIF 
  
  spawn, 'ls '+!SDSS_SHAPECORR_DIR,answer
  IF answer[0] EQ '' THEN BEGIN
      message,'Nothing Found!'
  ENDIF 
  
  print
  print,'-------------------------------'
  print,'Ignore type conversion errors'
  print,'-------------------------------'
  print
  nans=n_elements(answer)
  FOR i=0L, nans-1 DO IF strnumber(answer[i]) THEN add_arrval, long(answer[i]), runs

  nrun = n_elements(runs)

  s=create_struct("stripe", -1, "run", 0L, "rerun", -2, "bad", 0L)
  run_status = replicate(s, nrun)
  run_status.run = runs

  ;; link known to point at non-existing directory
;  w=where(runs NE 752 AND runs NE 745, ngood)
  ngood = n_elements(runs)
  w=lindgen(ngood)

  ;; check if link properly set
  FOR i=0L, ngood-1 DO BEGIN 
      fetch_rerun, run_status[w[i]].run, rer
      run_status[w[i]].rerun = rer
  ENDFOR 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; check if directories are OK
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  w=where(run_status.rerun GE 0, ngood)
  FOR i=0L, ngood-1 DO BEGIN 
      ind=w[i]

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; check for asTrans file
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      transdir = !SDSS_SHAPECORR_DIR+ ntostr(run_status[ind].run)+'/astrom/'
      spawn,'ls '+transdir+' | grep asTrans | wc -w',nastrans
      nastrans = long(nastrans[0])
      IF nastrans EQ 0 THEN run_status[ind].bad = run_status[ind].bad + 2L^0
      FOR camcol=1,6 DO BEGIN 
          fetch_dir, run_status[ind].run, camcol, run_status[ind].rerun, dir, atldir

          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
          ;; check for tsObj,fpAtlas,psField files
          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

          fetch_file_list, dir, files, fnums
          ;; If there are tsObj files, check for fpAtlas or psField files.
          spawn,'ls '+atldir+' | grep fpAtlas | wc -w',natlas
          spawn,'ls '+atldir+' | grep psField | wc -w',npsfield
          natlas = long(natlas[0])
          npsfield = long(npsfield[0])
          IF (files[0] EQ '') THEN BEGIN 
              run_status[ind].bad = run_status[ind].bad + 2L^(camcol)
          ENDIF ELSE BEGIN
              IF run_status[ind].stripe EQ -1 THEN BEGIN
                  hdr=headfits(files[0])
                  run_status[ind].stripe = sxpar(hdr, 'stripe')
              ENDIF 
          ENDELSE 
          IF natlas EQ 0 THEN run_status[ind].bad = run_status[ind].bad + 2L^(camcol+6)
          IF npsfield EQ 0 THEN run_status[ind].bad = run_status[ind].bad + 2L^(camcol+12)
      ENDFOR 
      jump:
  ENDFOR 

  print
  print,'!RUNCORR_STATUS_FILE: ',!RUNCORR_STATUS_FILE
  print
  mwrfits, run_status, !RUN_STATUS_FILE, /create
  
  return
END 
