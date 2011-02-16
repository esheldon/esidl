PRO makecat_scripts, runs, reruns, email=email

  IF n_params() LT 2 THEN BEGIN
      print,'-Syntax: makecat_scripts, runs, reruns, email=email'
      return
  ENDIF 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Parameters
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  sdssidl_setup,/silent
  setup_mystuff
  IF n_elements(email) EQ 0 THEN email = 'esheldon@sdss4.physics.lsa.umich.edu'

  nrun = n_elements(runs)
  nrerun = n_elements(reruns)

  IF nrun NE nrerun THEN message,'#runs ne #reruns'
  rstr = strarr(nrun)
  rrstr = strarr(nrun)
  combdir = strarr(nrun)

  ;; set up run strings and directories
  home = !SDSS_SHAPECORR_DIR
  FOR i=0L, nrun-1 DO BEGIN 
      rstr[i] = ntostr(runs[i])
      rrstr[i] = ntostr(reruns[i])
      combdir[i] = home+'corr'+rstr[i]+'/'+rrstr[i]+'/combined/'
      IF i EQ 0 THEN allruns = rstr[i] $
      ELSE allruns = allruns+'_' + rstr[i]
  ENDFOR 

  colors=['u','g','r','i','z']
  clr = [1,2,3]                 ;g,r,i
  nclr=n_elements(clr)

  head = '#! /bin/tcsh'

  ;; Put scripts in run1 directory
  rundir = home+'fscripts/'+rstr[0]
  script_dir = rundir+'/fscripts/'
  outdir = rundir+'/outfiles/'

  print
  print,' MAKECAT_SCRIPTS: Creating shape-correction scripts for'
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

  callscript = script_dir+'callmkcat'+allruns+'.sh'
  
  mkscatbase = strarr(nclr)
  mkscat_file = strarr(nclr)

  FOR i=0, nclr-1 DO BEGIN 

      mkscatbase[i] = 'mkscat'+allruns+'_'+colors[ clr[i] ]
      mkscat_file[i] = script_dir+mkscatbase[i] + '.sh'

  ENDFOR 

  mkspectrabase = 'mkspectra'+allruns
  mkspectrafile = script_dir+mkspectrabase+'.sh'

  print,' Creating call script: '
  print,'   '+callscript

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; The queues for each script
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  mkcatq = '4hr'

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Open main call script
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  
  openw, call_lun, callscript, /get_lun
  printf, call_lun, head
  printf, call_lun
  printf, call_lun, '. $FBATCH_DIR/bin/fbatch_setpgp.sh'


  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; make_scat scripts
  ;; want to run all colors at once
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  print,' Creating mkscat scripts: '
  colprint,'   '+mkscat_file

  FOR i=0, nclr-1 DO BEGIN 

      
      openw, lun, mkscat_file[i], /get_lun
      printf, lun, head

      FOR j=0L, nrun-1 DO BEGIN 
          printf, lun
          printf, lun, 'idl<<EOF'
          printf, lun, 'clr = '+ntostr(clr[i])
          printf, lun, 'run='+rstr[j]
          printf, lun, 'rerun='+rrstr[j]
          printf, lun, 'make_survey_scat, run, rerun, clr'
          printf, lun, 'EOF'
          printf, lun
      ENDFOR 

      free_lun, lun
      spawn,'chmod 755 '+mkscat_file[i]

  ENDFOR 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Print to call script for each bandpass
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  FOR i=0, nclr-1 DO BEGIN 
      
      printf, call_lun
      printf, call_lun, 'echo'
      printf, call_lun, 'echo Sending Script '+mkscatbase[i]+'.sh to '+$
                                                       mkcatq+' queue'
      printf, call_lun
      printf, call_lun, 'fbatch_sub -q '+mkcatq+ ' -R "fsgi03"'      $
        + ' -J '+mkscatbase[i]                                $
        + ' -o "' +outdir+mkscatbase[i]+'.out"'               $
        + ' -e "' +outdir+mkscatbase[i]+'.err"'               $
        + ' -N -u '+email     $
        + ' "' +mkscat_file[i]+ '"'
  ENDFOR 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; make spectra cats
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  print,' Creating mkspectra script: '
  print,'   '+mkspectrafile

  openw, lun, mkspectrafile, /get_lun
  printf, lun, head
  FOR j=0L, nrun-1 DO BEGIN 
      printf, lun
      printf, lun, 'idl<<EOF'
      printf, lun, 'run='+rstr[j]
      printf, lun, 'rerun='+rrstr[j]
      printf, lun, 'make_spectra_lcat, run, rerun'
      printf, lun, 'EOF'
      printf, lun
  ENDFOR 
  
  free_lun, lun
  spawn, 'chmod 755 '+mkspectrafile

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; write to call script
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  printf, call_lun
  printf, call_lun, 'echo'
  printf, call_lun, 'echo Sending Script '+mkspectrabase+'.sh to '+$
                                                       mkcatq+' queue'

  printf, call_lun
  printf, call_lun, 'fbatch_sub -q '+mkcatq+ ' -R "fsgi03"'      $
    + ' -J '+mkspectrabase                                       $
    + ' -o "' +outdir+mkspectrabase+'.out"'               $
    + ' -e "' +outdir+mkspectrabase+'.err"'               $
    + ' -N -u '+email     $
    + ' "' +mkspectrafile+ '"'

  ;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; free call script lun
  ;;;;;;;;;;;;;;;;;;;;;;;;;;

  free_lun, call_lun
  spawn, 'chmod 755 '+callscript
  
  print
  print,' MAKECAT_SCRIPTS: Excecution Successful'
  print
  print,' To make shape cats for Runs ',allruns,$
    ' on fsgi03, use these commands: '
  print,' setup fbatch'
  print,' '+callscript
  print

  return
END 
