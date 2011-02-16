PRO comb_scripts, runs, reruns, email=email

  IF n_params() LT 2 THEN BEGIN
      print,'-Syntax: comb_scripts, runs, reruns, email=email'
      return
  ENDIF 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Parameters
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  sdssidl_setup,/silent
  setup_mystuff
  IF n_elements(email) EQ 0 THEN email = 'esheldon@umich.edu'

  nrun = n_elements(runs)
  nrerun = n_elements(reruns)

  IF nrun NE nrerun THEN message,'#runs ne #reruns'
  rstr = strarr(nrun)
  rrstr = strarr(nrun)
  combdir = strarr(nrun)

  ;; string versions of runs,reruns
  home = !SDSS_SHAPECORR_DIR
  FOR i=0L, nrun-1 DO BEGIN 
      rstr[i] = ntostr(runs[i])
      rrstr[i] = ntostr(reruns[i])
      combdir[i] = home+'corr'+rstr[i]+'/'+rrstr[i]+'/combined/'
      IF i EQ 0 THEN allruns = rstr[i] $
      ELSE allruns = allruns + '_' + rstr[i]
  ENDFOR 

  ;; run,rerun arrays in string form
  runarrstr = '['
  rerunarrstr = '['
  FOR i=0L, nrun-1 DO BEGIN
      runarrstr = runarrstr + rstr[i]
      rerunarrstr = rerunarrstr + rrstr[i]
      IF i NE nrun-1 THEN BEGIN
          runarrstr = runarrstr + ', '
          rerunarrstr = rerunarrstr + ', '
      ENDIF 
  ENDFOR 
  runarrstr = runarrstr + ']'
  rerunarrstr = rerunarrstr + ']'

  colors=['u','g','r','i','z']
  clr = [1,2,3]                 ;g,r,i
  nclr=n_elements(clr)

  head = '#! /bin/tcsh'

  ;; Put scripts in run1 directory
  rundir = home+'fscripts/'+rstr[0]
  script_dir = rundir+'/fscripts/'
  outdir = rundir+'/outfiles/'

  print
  print,' COMB_SCRIPTS: Creating shape-correction and comb. scripts for'
  colprint,'               Run: '+rstr+' Rerun: '+rrstr
  print,' Putting these scripts in: ',script_dir

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; First check if this is a valid run/rerun
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  print
  print,' Checking if these are valid runs/reruns: ', format='(a,$)'

  FOR i=0L, nrun-1 DO BEGIN 
      fetch_dir, runs[i], 1, reruns[i], dir,atldir,/check
      IF (dir EQ '') OR (atldir EQ '') THEN $
        message,' ERROR: Not valid run/rerun: '+rstr[i]+'/'+rrstr[i]

      fetch_file_list,dir,files,fnums
      IF files[0] EQ '' THEN $
        message,'ERROR: Not valid run/rerun: '+rstr+'/'+rrstr
  ENDFOR 
  print,'OK'

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Make sure corrected directories exist
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  FOR i=0L, nrun-1 DO BEGIN 
      print
      print,' Checking for '+rstr[i]+' combined file directory: ',format='(a,$)'
      spawn,'ls -l '+combdir[i], answer
      IF answer[0] EQ '' THEN BEGIN
          print
          print,' Making combined directory ',combdir1
          spawn,'mkdir '+combdir1
      ENDIF ELSE print,'Ok'
  ENDFOR 
  print

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; callscript: sends scripts to the appropriate queues with
  ;; the right dependencies.
  ;; The other scripts run the IDL routines
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  callscript = script_dir+'callcomb'+allruns+'.sh'
  
  combinebase    = 'combine'+allruns
  combine_file = script_dir+combinebase+'.sh'

  combspectrabase = 'combspec'+allruns
  combspectrafile = script_dir+combspectrabase+'.sh'

  print,' Creating call script: '
  print,'   '+callscript

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; The queues for each script
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  combineq = '30min'

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Open main call script
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  
  openw, call_lun, callscript, /get_lun
  printf, call_lun, head
  printf, call_lun
  printf, call_lun, '. $FBATCH_DIR/bin/fbatch_setpgp.sh'

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; combine script
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  print,' Creating combine script: '
  print,'   '+combine_file

  openw, lun, combine_file, /get_lun
  printf, lun, head
  printf, lun

  FOR i=0, nclr-1 DO BEGIN 

      printf, lun, 'idl<<EOF'
      printf, lun, 'clr='+ntostr(clr[i])
      printf, lun, 'runs = '+runarrstr
      printf, lun, 'reruns = '+rerunarrstr
      printf, lun, 'combine_stripe, runs, reruns, clr'
      printf, lun, 'EOF'
      printf, lun, 'idl<<EOF'
      printf, lun, 'clr='+ntostr(clr[i])
      printf, lun, 'runs = '+runarrstr
      printf, lun, 'reruns = '+rerunarrstr
      printf, lun, 'combine_stripe, runs, reruns, clr, /lcat'
      printf, lun, 'EOF'

  ENDFOR 
  free_lun, lun
  spawn,'chmod 755 '+combine_file

  ;; write to call script

  printf, call_lun
  printf, call_lun, 'echo'
  printf, call_lun,'echo Sending Script '+combinebase+'.sh to '+combineq+' queue'
  printf, call_lun
  printf, call_lun, 'fbatch_sub -q '+combineq+ ' -R "fsgi03"'    $
    + ' -J '+combinebase                              $
    + ' -o "' +outdir+combinebase+'.out"'             $
    + ' -e "' +outdir+combinebase+'.err"'             $
    + ' -N -u '+email  $
    + ' "' +combine_file+ '"'


  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; combine spectra script
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  print,' Creating combine spectra script: '
  print,'   '+combspectrafile

  openw, lun, combspectrafile, /get_lun
  printf, lun, head
  printf, lun

  printf, lun, 'idl<<EOF'
  printf, lun, 'runs = '+runarrstr
  printf, lun, 'reruns = '+rerunarrstr
  printf, lun, 'combine_spectra_lcat, runs, reruns'
  printf, lun, 'EOF'

  free_lun, lun
  spawn,'chmod 755 '+combspectrafile

  ;; write to call script

  printf, call_lun
  printf, call_lun, 'echo'
  printf, call_lun,'echo Sending Script '+combspectrabase+'.sh to '+combineq+' queue'
  printf, call_lun
  printf, call_lun, 'fbatch_sub -q '+combineq+ ' -R "fsgi03"'    $
    + ' -J '+combspectrabase                              $
    + ' -o "' +outdir+combspectrabase+'.out"'             $
    + ' -e "' +outdir+combspectrabase+'.err"'             $
    + ' -N -u '+email  $
    + ' "' +combspectrafile+ '"'

  ;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; free call script lun
  ;;;;;;;;;;;;;;;;;;;;;;;;;;

  free_lun, call_lun
  spawn, 'chmod 755 '+callscript
  
  print
  print,' COMB_SCRIPTS: Excecution Successful'
  print
  print,' To combine Runs ',allruns,$
    ' on fsgi03, use these commands: '
  print,' setup fbatch'
  print,' '+callscript
  print

  return
END 
