PRO print_status_check, lun

  printf, lun, 'status=$?'
  printf, lun, 'if [ $status -ne 0 ]'
  printf, lun, 'then'
  printf, lun, '    echo Error in `basename $0`: Halting Execution'
  printf, lun, '    exit $status'
  printf, lun, 'fi'
  printf, lun

END 

PRO make_corrscripts, run, rerun, maxmag, $
                      scriptDir=scriptDir, $
                      outFileDir=outFileDir, $
                      disk=disk, $
                      addphotoz=addphotoz, addbayes=addbayes, $
                      nonewcrun=nonewcrun

  IF n_params() LT 3 THEN BEGIN
      print,'-Syntax: make_corrscripts, run, rerun, maxmag, $'
      print,'           scriptDir=scriptDir, $'
      print,'           outFileDir=outFileDir, $'
      print,'           disk=disk, /addphotoz, /addbayes, /nonewcrun'
      print
      print,'Use /nonewcrun if the corrected directory already'
      print,'exists on another machine'
      print,'Note that he "disk" parameter is irrelevant in that case'
      return
  ENDIF 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Parameters
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF n_elements(disk) EQ 0 THEN disk='data0'

  setup_mystuff
  binDir = !MYBINDIR

  rstr = ntostr(run)
  rrstr =  ntostr(rerun)
  maxmagstr = ntostr(maxmag)

  head = '#! /bin/sh'
;  head = '#!/bin/tcsh'
  camcol = ['1','2','3','4','5','6']
  ncol = n_elements(camcol)

  home = sdssidl_config('SHAPECORR_DIR')
  runDir     = home+'scripts/'+rstr + '/'
  rerunDir   = runDir + rrstr +'/'
  scriptDir = rerunDir+'scripts/'
  outFileDir   = rerunDir+'outfiles/'

  IF keyword_set(addphotoz) THEN photozstr = ', /addphotoz' ELSE photozstr=''
  IF keyword_set(addbayes) THEN bayesstr = ', /addbayes' ELSE bayesstr=''

  print
  print,' CORR_SCRIPTS: Creating shape-correction scripts for'
  print,'               Run: ',rstr,' Rerun: ',rrstr

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; First check if this is a valid run/rerun
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  print
  print,' Checking if this is a valid run/rerun: ', format='(a,$)'

  input_index = where(!run_status.tsObj_photo_v GE 5.4)
  get_goodruns, gruns, greruns, gstripes, gstrips, gindices, $
    /silent, /tsObj, /asTrans, input_index=input_index

  w=where((gruns EQ run) AND (greruns EQ rerun), nw)
  IF nw EQ 0 THEN BEGIN 
      message,'ERROR: Invalid run/rerun: '+rstr+'/'+rrstr
  ENDIF 
  print,' OK'

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Make sure corrected directories exist
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF NOT keyword_set(nonewcrun) THEN BEGIN 
      print,' Setting up corrected file directories'
      print,' --------------------------------------------------------'
      command = binDir+'newcrun_home '+rstr+' '+rrstr+' '+disk
      spawn,command
      print,' --------------------------------------------------------'
      print
  ENDIF ELSE BEGIN 
      print,' --------------------------------------------------------'
      print,' Not doing newcrun_home'
      print,' --------------------------------------------------------'
      print
  ENDELSE 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Make script directories (if they don't already exist)
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  print,' Setting up script directories: ', format='(a,$)'
  print
  IF NOT file_test(runDir, /directory) THEN BEGIN 
      print,' Making directory ',runDir
      file_mkdir, runDir
  ENDIF 
  IF NOT file_test(rerunDir, /directory) THEN BEGIN 
      print,' Making directory ',rerunDir
      file_mkdir, rerunDir
  ENDIF 
  IF NOT file_test(scriptDir, /directory) THEN BEGIN 
      print,' Making directory ',scriptDir
      file_mkdir, scriptDir
  ENDIF 
  IF NOT file_test(outFileDir, /directory) THEN BEGIN 
      print,' Making directory ',outFileDir
      file_mkdir, outFileDir
  ENDIF 


  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; callscript: sends scripts to the appropriate queues with
  ;; the right dependencies.
  ;; The other scripts run the routines
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;  rotationbase = 'rotate'+rstr
  rsmearbase = 'rsmear_cuts'+rstr
  makecorbase = 'make_corrected'+rstr
  corrtestbase = 'corrtest'+rstr

  call_123script = scriptDir+'correct'+rstr+'-col1-2-3.sh'
  call_456script = scriptDir+'correct'+rstr+'-col4-5-6.sh'

;  rotation_123file = scriptDir + rotationbase+'-col1-2-3.sh'
;  rotation_456file = scriptDir + rotationbase+'-col4-5-6.sh'

  makecor_123file = scriptDir+makecorbase+'-col1-2-3.sh'
  makecor_456file = scriptDir+makecorbase+'-col4-5-6.sh'

  rsmear_file = scriptDir+rsmearbase+'.sh'

  corrtest_file = scriptDir+corrtestbase+'.sh'
  corrtest_123file = scriptDir+corrtestbase+'-col1-2-3.sh'
  corrtest_456file = scriptDir+corrtestbase+'-col4-5-6.sh'

  FOR icol=1, 6 DO BEGIN 

      camstr = ntostr(icol)

      IF (icol EQ 1) OR (icol EQ 4) THEN BEGIN 

          IF icol EQ 1 THEN BEGIN 
              callscript = call_123script
;              rotation_file = rotation_123file
              makecor_file = makecor_123file
              tcorrtest_file = corrtest_123file
          ENDIF ELSE BEGIN 
              callscript = call_456script
;              rotation_file = rotation_456file
              makecor_file = makecor_456file
              tcorrtest_file = corrtest_456file
          ENDELSE 
          
          print
          print
          print,' Creating call script: '
          print,'    '+callscript

          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
          ;; Open scripts
          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  
          openw, call_lun, callscript, /get_lun
          printf, call_lun, head
          printf, call_lun

;          openw, rot_lun, rotation_file, /get_lun 
;          printf, rot_lun, head
          openw, makecor_lun, makecor_file, /get_lun 
          printf, makecor_lun, head

      ENDIF 

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; rotation script
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;      IF (icol EQ 1) OR (icol EQ 4) THEN BEGIN 
;          print,' Creating rotation script: '
;          print,'    '+rotation_file
;      ENDIF 
;      printf, rot_lun
;      printf, rot_lun, 'idl<<EOF'
;      printf, rot_lun, 'setup_mystuff'
;      printf, rot_lun, 'run=',rstr
;      printf, rot_lun, 'rerun=',rrstr
;      printf, rot_lun, 'camcol='+camstr
;      printf, rot_lun, 'sdss_survey_rot,run,rerun,camcol, status=status'
;      printf, rot_lun, 'IF status ne 0 THEN exit,status=45'
;      printf, rot_lun, 'EOF'
;      print_status_check, rot_lun

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; makecor script
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      IF (icol EQ 1) OR (icol EQ 4) THEN BEGIN 
          print,' Creating makecor script: '
          print,'    '+makecor_file
      ENDIF 

      printf, makecor_lun
      printf, makecor_lun, 'idl<<EOF'
      printf, makecor_lun, 'setup_mystuff'
      printf, makecor_lun, 'maxmag='+maxmagstr
      printf, makecor_lun, 'run=',rstr
      printf, makecor_lun, 'rerun=',rrstr
      printf, makecor_lun, 'camcol=',camstr
      printf, makecor_lun
      printf, makecor_lun, 'make_corrected_files,maxmag,run,rerun,camcol,status=status'+photozstr+bayesstr
      printf, makecor_lun, 'IF status ne 0 THEN exit, status=45'
      printf, makecor_lun, 'EOF'
      print_status_check, makecor_lun

      IF (icol EQ 3) OR (icol EQ 6) THEN BEGIN 

          ;; now print to call script
;          printf, call_lun, rotation_file
;          print_status_check, call_lun
          printf, call_lun, makecor_file
          print_status_check, call_lun
          printf, call_lun, tcorrtest_file
          print_status_check, call_lun

          printf, call_lun
          printf, call_lun, 'echo Done'

;          free_lun, rot_lun & spawn,'chmod 755 '+rotation_file
          free_lun, makecor_lun & spawn,'chmod 755 '+makecor_file

          free_lun, call_lun & spawn,'chmod 755 '+callscript

      ENDIF 

  ENDFOR 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; corrtest scripts
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  
;  print,' Creating corrtest script: '
;  print,'    '+corrtest_file

;  openw, lun, corrtest_file, /get_lun
;  printf, lun, head
;  printf, lun
 ; printf, lun, 'idl<<EOF'
;  printf, lun, 'setup_mystuff'
;  printf, lun, 'run=',rstr
;  printf, lun, 'rerun=',rrstr
;  printf, lun, 'clr=1'
;  printf, lun, 'corrtest_plots_allcols, run, rerun, clr, status=status'
;  printf, lun, 'IF status ne 0 THEN exit, status=45'
;  printf, lun, 'clr=2'
;  printf, lun, 'corrtest_plots_allcols, run, rerun, clr, status=status'
;  printf, lun, 'IF status ne 0 THEN exit, status=45'
;  printf, lun, 'clr=3'
;  printf, lun, 'corrtest_plots_allcols, run, rerun, clr, status=status'
;  printf, lun, 'IF status ne 0 THEN exit, status=45'
;  printf, lun, 'EOF'
;  print_status_check, lun
;  free_lun, lun
;  spawn,'chmod 755 '+corrtest_file

  openw, lun, corrtest_123file, /get_lun
  printf, lun, head
  printf, lun
  printf, lun, 'idl<<EOF'
  printf, lun, 'setup_mystuff'
  printf, lun, 'run=',rstr
  printf, lun, 'rerun=',rrstr
  printf, lun, 'camcol=1'
  printf, lun, 'clr=1'
  printf, lun, 'corrtest_plots, run, rerun, camcol, clr, status=status'
  printf, lun, 'IF status ne 0 THEN exit, status=45'
  printf, lun, 'clr=2'
  printf, lun, 'corrtest_plots, run, rerun, camcol, clr, status=status'
  printf, lun, 'IF status ne 0 THEN exit, status=45'
  printf, lun, 'clr=3'
  printf, lun, 'corrtest_plots, run, rerun, camcol, clr, status=status'
  printf, lun, 'IF status ne 0 THEN exit, status=45'
  printf, lun, 'camcol=2'
  printf, lun, 'clr=1'
  printf, lun, 'corrtest_plots, run, rerun, camcol, clr, status=status'
  printf, lun, 'IF status ne 0 THEN exit, status=45'
  printf, lun, 'clr=2'
  printf, lun, 'corrtest_plots, run, rerun, camcol, clr, status=status'
  printf, lun, 'IF status ne 0 THEN exit, status=45'
  printf, lun, 'clr=3'
  printf, lun, 'corrtest_plots, run, rerun, camcol, clr, status=status'
  printf, lun, 'IF status ne 0 THEN exit, status=45'
  printf, lun, 'camcol=3'
  printf, lun, 'clr=1'
  printf, lun, 'corrtest_plots, run, rerun, camcol, clr, status=status'
  printf, lun, 'IF status ne 0 THEN exit, status=45'
  printf, lun, 'clr=2'
  printf, lun, 'corrtest_plots, run, rerun, camcol, clr, status=status'
  printf, lun, 'IF status ne 0 THEN exit, status=45'
  printf, lun, 'clr=3'
  printf, lun, 'corrtest_plots, run, rerun, camcol, clr, status=status'
  printf, lun, 'IF status ne 0 THEN exit, status=45'
  printf, lun, 'EOF'
  print_status_check, lun
  free_lun, lun
  spawn,'chmod 755 '+corrtest_123file
  
  openw, lun, corrtest_456file, /get_lun
  printf, lun, head
  printf, lun
  printf, lun, 'idl<<EOF'
  printf, lun, 'setup_mystuff'
  printf, lun, 'run=',rstr
  printf, lun, 'rerun=',rrstr
  printf, lun, 'camcol=4'
  printf, lun, 'clr=1'
  printf, lun, 'corrtest_plots, run, rerun, camcol, clr, status=status'
  printf, lun, 'IF status ne 0 THEN exit, status=45'
  printf, lun, 'clr=2'
  printf, lun, 'corrtest_plots, run, rerun, camcol, clr, status=status'
  printf, lun, 'IF status ne 0 THEN exit, status=45'
  printf, lun, 'clr=3'
  printf, lun, 'corrtest_plots, run, rerun, camcol, clr, status=status'
  printf, lun, 'IF status ne 0 THEN exit, status=45'
  printf, lun, 'camcol=5'
  printf, lun, 'clr=1'
  printf, lun, 'corrtest_plots, run, rerun, camcol, clr, status=status'
  printf, lun, 'IF status ne 0 THEN exit, status=45'
  printf, lun, 'clr=2'
  printf, lun, 'corrtest_plots, run, rerun, camcol, clr, status=status'
  printf, lun, 'IF status ne 0 THEN exit, status=45'
  printf, lun, 'clr=3'
  printf, lun, 'corrtest_plots, run, rerun, camcol, clr, status=status'
  printf, lun, 'IF status ne 0 THEN exit, status=45'
  printf, lun, 'camcol=6'
  printf, lun, 'clr=1'
  printf, lun, 'corrtest_plots, run, rerun, camcol, clr, status=status'
  printf, lun, 'IF status ne 0 THEN exit, status=45'
  printf, lun, 'clr=2'
  printf, lun, 'corrtest_plots, run, rerun, camcol, clr, status=status'
  printf, lun, 'IF status ne 0 THEN exit, status=45'
  printf, lun, 'clr=3'
  printf, lun, 'corrtest_plots, run, rerun, camcol, clr, status=status'
  printf, lun, 'IF status ne 0 THEN exit, status=45'
  printf, lun, 'EOF'
  print_status_check, lun
  free_lun, lun
  spawn,'chmod 755 '+corrtest_456file

  print
  print,' CORR_SCRIPTS: Excecution Successful'
  print
  print,' To correct Run ',rstr,' use these commands: '
  print,'   '+call_123script
  print,'   '+call_456script
  print
  print,' When those are done, use these commands to make QA plots: '
  ;;print,'   '+corrtest_file
  print,'   '+corrtest_123file
  print,'   '+corrtest_456file
  print

  return
END 
