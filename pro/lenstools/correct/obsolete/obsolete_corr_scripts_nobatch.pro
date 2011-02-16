PRO print_status_check, lun

  printf, lun, 'status=$?'
  printf, lun, 'if [ $status -ne 0 ]'
  printf, lun, 'then'
  printf, lun, '    echo Error in `basename $0`: Halting Execution'
  printf, lun, '    exit $status'
  printf, lun, 'fi'
  printf, lun

END 

PRO corr_scripts_nobatch, run, rerun, color_index, email=email, nocheck=nocheck, $
                        psFieldrerun=psFieldrerun

  IF n_params() LT 3 THEN BEGIN
      print,'-Syntax: corr_scripts_sdss5, run, rerun, color_index, email=email, nocheck=nocheck, psFieldrerun=psFieldrerun'
      print,'/nocheck does not check if #atlas files = #tsObj files'
      ;print,' If you do not give the email address, then email="esheldon@umich.edu"'
      print,'We use color_index = 2'
      return
  ENDIF 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Parameters
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF n_elements(psFieldrerun) EQ 0 THEN psFieldrerun=rerun

  IF n_elements(email) EQ 0 THEN email = 'esheldon@umich.edu'

  setup_mystuff
  bindir = !MYBINDIR

  rstr = ntostr(run)
  rrstr =  ntostr(rerun)
  psrrstr = ntostr(psFieldrerun)
  head = '#! /bin/sh'
;  head = '#!/bin/tcsh'
  camcol = ['1','2','3','4','5','6']
  clr = ntostr(color_index)
  ncol = n_elements(camcol)

  sdssidl_setup, /silent
  setup_mystuff

  home = !SDSS_SHAPECORR_DIR
  rundir = home+'scripts/'+rstr
  script_dir = rundir+'/scripts/'
  outdir = rundir+'/outfiles/'

  print
  print,' CORR_SCRIPTS: Creating shape-correction scripts for'
  print,'               Run: ',rstr,' Rerun: ',rrstr

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; First check if this is a valid run/rerun
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  print
  print,' Checking if this is a valid run/rerun: ', format='(a,$)'
  fetch_dir, run, 1, rerun, dir,atldir,/check
  IF (dir EQ '') OR (atldir EQ '') THEN BEGIN
      print,' ERROR: Not valid run/rerun: '+rstr+'/'+rrstr
      return
  ENDIF 
  fetch_file_list,dir,files,fnums,start=start,nframes=nframes
  IF files[0] EQ '' THEN BEGIN
      print,'ERROR: Not valid run/rerun: '+rstr+'/'+rrstr
      return
  ENDIF 
  print,' OK'
  ntot = n_elements(fnums)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; See if target has been run
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  notarget = 'NOCVS:ts'
  chhdr = headfits(files[0])
  chtarget = sxpar(chhdr, 'TARG_VER')
  IF datatype(chtarget) EQ 'STR' THEN BEGIN 
      IF ntostr(chtarget) EQ notarget THEN BEGIN 
          print,' Target not run on these files. Want to continue (y/n)?'
          key=get_kbrd(1)
          CASE key OF
              'n': return
              'N': return
              ELSE: 
          ENDCASE 
      ENDIF ELSE BEGIN
          print
          print,' Version of target: ',chtarget
          print
      ENDELSE 
  ENDIF ELSE BEGIN 
      print,' Header value "TARG_VER" not present. Want to continue (y/n)?'
      key=get_kbrd(1)
      CASE key OF
          'n': return
          'N': return
          ELSE:
      ENDCASE 
  ENDELSE 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Check if number of tsObj files equals number of fpAtlas files
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  IF NOT keyword_set(nocheck) THEN BEGIN 
      print,' Checking if #fpAtlas >= #tsObj: ', format='(a,$)'
      spawn,'ls '+atldir+' | grep fpAtlas | wc -w', natlas
      natlas = long(natlas[0])
      IF natlas LT ntot THEN BEGIN
          print,' ERROR: #fpAtlas < #tsobj '+ntostr(ntot)+'/'+ntostr(natlas)
          return
      ENDIF 
      print,' OK  '+ntostr(natlas)+'/'+ntostr(ntot)
      print
  ENDIF ELSE BEGIN 
      print,' Not checking #fpAtlas >= #tsObj'
      print
  ENDELSE 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; make sure photo version is the same for tsObj and fpAtlas
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  print, ' Checking if PHOTO versions same for tsObj and fpAtlas: ',format='(a,$)'
  q=where( (!run_status.run EQ run) AND $
           (!run_status.rerun EQ rerun), nq)
  IF nq EQ 0 THEN message,'WHAT!!!'
  
  IF (!run_status[q[0]].tsobj_photo_v NE $
      !run_status[q[0]].fpatlas_photo_v) THEN BEGIN 
      message,'Photo versions differ!!!',/inf
      message,'  tsObj PHOTO VERSION = '+!run_status[q[0]].tsobj_photo_v,/inf
      message,'  fpAtlas PHOTO VERSION = '+!run_status[q[0]].fpatlas_photo_v
  ENDIF 
  print,' OK  '+!run_status[q[0]].tsobj_photo_v+'  '+$
        !run_status[q[0]].fpatlas_photo_v
  print

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Make sure corrected directories exist
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  print,' Setting up corrected file directories'
  print,' --------------------------------------------------------'
  command = bindir+'newcrun_home '+rstr+' '+rrstr
  spawn,command
  print,' --------------------------------------------------------'
  print

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Make script directories (if they don't already exist)
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  print,' Setting up script directories: ', format='(a,$)'
  spawn,['ls','-l',rundir], answer1,count=count1,/noshell
  pnext = 'Directories already exist'
  IF count1 EQ 0 THEN BEGIN
      pnext=''
      print
      print,' Making directory ',rundir
      spawn,'mkdir '+rundir
      print,' Making directory ',script_dir
      spawn,'mkdir '+script_dir
      print,' Making directory ',outdir
      spawn,'mkdir '+outdir
  ENDIF ELSE BEGIN 
      spawn,['ls','-l',script_dir],answer2,count=count2,/noshell
      spawn,['ls','-l',outdir],answer3,count=count3,/noshell
      IF strlen(answer2[0]) EQ 0 THEN BEGIN
          print
          pnext=''
          print,' Making directory ',script_dir
          spawn,'mkdir '+script_dir
      ENDIF 
      IF strlen(answer3[0]) EQ 0 THEN BEGIN
          IF pnext EQ '' THEN print
          pnext=''
          print,' Making directory ',outdir
          spawn,'mkdir '+outdir
      ENDIF
  ENDELSE 
  IF pnext EQ '' THEN print ELSE BEGIN
      print,pnext
      print
  ENDELSE 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; callscript: sends scripts to the appropriate queues with
  ;; the right dependencies.
  ;; The other scripts run the routines
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  rotationbase = 'rotate'+rstr
  admomIDLbase = 'admomIDL'+rstr
  admomCbase = 'admomC'+rstr
  convertbase = 'convert'+rstr
  fitmombase = 'fitmom'+rstr
  makecorbase = 'makecor'+rstr
  corshapebase = 'corshape'+rstr

  call_123script = script_dir+'correct'+rstr+'-col1-2-3.sh'
  call_456script = script_dir+'correct'+rstr+'-col4-5-6.sh'

  rotation_123file = script_dir + rotationbase+'-col1-2-3.sh'
  rotation_456file = script_dir + rotationbase+'-col4-5-6.sh'

  admomIDL_123file = script_dir + admomIDLbase + '-col1-2-3.sh'
  admomIDL_456file = script_dir + admomIDLbase + '-col4-5-6.sh'

  admomC_123file = script_dir + admomCbase+'-col1-2-3.sh'
  admomC_456file = script_dir + admomCbase+'-col4-5-6.sh'

  convert_123file = script_dir + convertbase+'-col1-2-3.sh'
  convert_456file = script_dir + convertbase+'-col4-5-6.sh'

  fitmom_123file = script_dir + fitmombase+'-col1-2-3.sh'
  fitmom_456file = script_dir + fitmombase+'-col4-5-6.sh'

  makecor_123file = script_dir+makecorbase+'-col1-2-3.sh'
  makecor_456file = script_dir+makecorbase+'-col4-5-6.sh'

  corshape_file = script_dir+corshapebase+'.sh'
  corshape_123file = script_dir+corshapebase+'-col1-2-3.sh'
  corshape_456file = script_dir+corshapebase+'-col4-5-6.sh'

  FOR icol=1, 6 DO BEGIN 

      camstr = ntostr(icol)

      IF (icol EQ 1) OR (icol EQ 4) THEN BEGIN 

          IF icol EQ 1 THEN BEGIN 
              callscript = call_123script
              rotation_file = rotation_123file
              admomIDL_file = admomIDL_123file
              admomC_file = admomC_123file
              convert_file = convert_123file
              fitmom_file = fitmom_123file
              makecor_file = makecor_123file
              tcorshape_file = corshape_123file
          ENDIF ELSE BEGIN 
              callscript = call_456script
              rotation_file = rotation_456file
              admomIDL_file = admomIDL_456file
              admomC_file = admomC_456file
              convert_file = convert_456file
              fitmom_file = fitmom_456file
              makecor_file = makecor_456file
              tcorshape_file = corshape_456file
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

          openw, rot_lun, rotation_file, /get_lun      & printf, rot_lun, head
          openw, admomIDL_lun, admomIDL_file, /get_lun & printf, admomIDL_lun, head
          openw, admomC_lun, admomC_file, /get_lun     & printf, admomC_lun, head
          openw, convert_lun, convert_file, /get_lun   & printf, convert_lun, head
          openw, fitmom_lun, fitmom_file, /get_lun     & printf, fitmom_lun, head
          openw, makecor_lun, makecor_file, /get_lun   & printf, makecor_lun, head

      ENDIF 

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; rotation script
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


      IF (icol EQ 1) OR (icol EQ 4) THEN BEGIN 
          print,' Creating rotation script: '
          print,'    '+rotation_file
      ENDIF 
      printf, rot_lun
      printf, rot_lun, 'idl<<EOF'
      printf, rot_lun, 'setup_mystuff'
      printf, rot_lun, 'run=',rstr
      printf, rot_lun, 'rerun=',rrstr
      printf, rot_lun, 'camcol='+camstr
      printf, rot_lun, 'sdss_survey_rot,run,rerun,camcol, status=status'
      printf, rot_lun, 'IF status ne 0 THEN exit,status=45'
      printf, rot_lun, 'EOF'
      print_status_check, rot_lun

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; admom IDL script
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      IF (icol EQ 1) OR (icol EQ 4) THEN BEGIN 
          print,' Creating admom IDL script: '
          print,'    '+admomIDL_file
      ENDIF 
      printf, admomIDL_lun
      printf, admomIDL_lun, 'idl<<EOF'
      printf, admomIDL_lun, 'setup_mystuff'
      printf, admomIDL_lun, 'run=',rstr
      printf, admomIDL_lun, 'rerun=',rrstr
      printf, admomIDL_lun, 'psFieldrerun='+psrrstr
;      printf, admomIDL_lun, 'start=',ntostr(start)
;      printf, admomIDL_lun, 'nf=',ntostr(nframes)
      printf, admomIDL_lun
      printf, admomIDL_lun, 'color_index=',clr
      printf, admomIDL_lun
      printf, admomIDL_lun, 'camcol='+camstr
;      printf, admomIDL_lun, 'admomatlas, color_index, run, rerun, camcol, start=start,nf=nf,psFieldrerun=psFieldrerun'
      printf, admomIDL_lun, 'admomatlas, color_index, run, rerun, camcol, psFieldrerun=psFieldrerun, status=status'
      printf, admomIDL_lun, 'IF status ne 0 THEN exit,status=45'
      printf, admomIDL_lun, 'EOF'
      print_status_check, admomIDL_lun

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; admomatlas C scripts
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      IF (icol EQ 1) OR (icol EQ 4) THEN BEGIN 
          print,' Creating admom C scripts: '
          colprint,'    '+admomC_file
      ENDIF 
      
      printf, admomC_lun
      fetch_dir, run, icol, rerun, tdir, atldir, corrdir=corrdir
      admomatlas_infile, run, rerun, icol, adinfile, adoutfile
      adinfile = corrdir+adinfile & adoutfile = corrdir+adoutfile
      
      printf, admomC_lun, 'admomatlas '+adinfile+' '+adoutfile+' '+atldir
      print_status_check, admomC_lun

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; convert admomout files from ascii to fits
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      IF (icol EQ 1) OR (icol EQ 4) THEN BEGIN 
          print,' Creating convert script: '
          print,'    '+convert_file
      ENDIF 
  
      printf, convert_lun
      printf, convert_lun, 'idl<<EOF'
      printf, convert_lun, 'setup_mystuff'
      printf, convert_lun

      fetch_dir, run, icol, rerun, tdir, atldir, corrdir=corrdir
      admomatlas_infile, run, rerun, icol, adinfile, adoutfile
      adoutfile = corrdir+adoutfile

      printf, convert_lun, 'infile = "'+adoutfile+'"'
      printf, convert_lun, 'admomascii2fits, infile,status=status'
      printf, convert_lun, 'IF status ne 0 THEN exit, status=45'
      printf, convert_lun, 'EOF'
      print_status_check, convert_lun

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; fitmom script. New way is slow, break it up to 6 scripts
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      IF (icol EQ 1) OR (icol EQ 4) THEN BEGIN 
          print,' Creating fitmom scripts: '
          colprint,'    '+fitmom_file
      ENDIF 

      printf, fitmom_lun
      printf, fitmom_lun, 'idl<<EOF'
      printf, fitmom_lun, 'setup_mystuff'
      printf, fitmom_lun, 'run=',rstr
      printf, fitmom_lun, 'camcol='+camstr
      printf, fitmom_lun, 'rerun=',rrstr
      printf, fitmom_lun, 'psFieldrerun='+psrrstr
      printf, fitmom_lun
      printf, fitmom_lun, 'fitmom, run, rerun, camcol,psFieldrerun=psFieldrerun, status=status'
      printf, fitmom_lun, 'IF status ne 0 THEN exit, status=45'
      printf, fitmom_lun, 'EOF'
      print_status_check, fitmom_lun

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
      printf, makecor_lun, 'run=',rstr
      printf, makecor_lun, 'rerun=',rrstr
      printf, makecor_lun, 'camcol=',camstr
;      printf, makecor_lun, 'start=',ntostr(start)
;      printf, makecor_lun, 'nf=',ntostr(ntot)
      printf, makecor_lun
;      printf, makecor_lun, 'makecor,run,rerun,camcol,start=start,nframes=nf'
      printf, makecor_lun, 'makecor,run,rerun,camcol,status=status'
      printf, makecor_lun, 'IF status ne 0 THEN exit, status=45'
      printf, makecor_lun, 'EOF'
      print_status_check, makecor_lun

      IF (icol EQ 3) OR (icol EQ 6) THEN BEGIN 

          ;; now print to call script
          printf, call_lun, rotation_file
          print_status_check, call_lun
          printf, call_lun, admomIDL_file
          print_status_check, call_lun
          printf, call_lun, admomC_file
          print_status_check, call_lun
          printf, call_lun, convert_file
          print_status_check, call_lun
          printf, call_lun, fitmom_file
          print_status_check, call_lun
          printf, call_lun, makecor_file
          print_status_check, call_lun
          printf, call_lun, tcorshape_file
          print_status_check, call_lun

          free_lun, rot_lun & spawn,'chmod 755 '+rotation_file
          free_lun, admomIDL_lun & spawn,'chmod 755 '+admomIDL_file
          free_lun, admomC_lun & spawn,'chmod 755 '+admomC_file
          free_lun, convert_lun & spawn,'chmod 755 '+convert_file
          free_lun, fitmom_lun & spawn,'chmod 755 '+fitmom_file
          free_lun, makecor_lun & spawn,'chmod 755 '+makecor_file

          free_lun, call_lun & spawn,'chmod 755 '+callscript

      ENDIF 

  ENDFOR 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; corshape scripts
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  
  print,' Creating corshape script: '
  print,'    '+corshape_file

  openw, lun, corshape_file, /get_lun
  printf, lun, head
  printf, lun
  printf, lun, 'idl<<EOF'
  printf, lun, 'setup_mystuff'
  printf, lun, 'run=',rstr
  printf, lun, 'rerun=',rrstr
  printf, lun, 'clr=1'
  printf, lun, 'corshape_allcols, run, rerun, clr, status=status'
  printf, lun, 'IF status ne 0 THEN exit, status=45'
  printf, lun, 'clr=2'
  printf, lun, 'corshape_allcols, run, rerun, clr, status=status'
  printf, lun, 'IF status ne 0 THEN exit, status=45'
  printf, lun, 'clr=3'
  printf, lun, 'corshape_allcols, run, rerun, clr, status=status'
  printf, lun, 'IF status ne 0 THEN exit, status=45'
  printf, lun, 'EOF'
  print_status_check, lun
  free_lun, lun
  spawn,'chmod 755 '+corshape_file

  openw, lun, corshape_123file, /get_lun
  printf, lun, head
  printf, lun
  printf, lun, 'idl<<EOF'
  printf, lun, 'setup_mystuff'
  printf, lun, 'run=',rstr
  printf, lun, 'rerun=',rrstr
  printf, lun, 'camcol=1'
  printf, lun, 'clr=1'
  printf, lun, 'corshape, run, rerun, camcol, clr, status=status'
  printf, lun, 'IF status ne 0 THEN exit, status=45'
  printf, lun, 'clr=2'
  printf, lun, 'corshape, run, rerun, camcol, clr, status=status'
  printf, lun, 'IF status ne 0 THEN exit, status=45'
  printf, lun, 'clr=3'
  printf, lun, 'corshape, run, rerun, camcol, clr, status=status'
  printf, lun, 'IF status ne 0 THEN exit, status=45'
  printf, lun, 'camcol=2'
  printf, lun, 'clr=1'
  printf, lun, 'corshape, run, rerun, camcol, clr, status=status'
  printf, lun, 'IF status ne 0 THEN exit, status=45'
  printf, lun, 'clr=2'
  printf, lun, 'corshape, run, rerun, camcol, clr, status=status'
  printf, lun, 'IF status ne 0 THEN exit, status=45'
  printf, lun, 'clr=3'
  printf, lun, 'corshape, run, rerun, camcol, clr, status=status'
  printf, lun, 'IF status ne 0 THEN exit, status=45'
  printf, lun, 'camcol=3'
  printf, lun, 'clr=1'
  printf, lun, 'corshape, run, rerun, camcol, clr, status=status'
  printf, lun, 'IF status ne 0 THEN exit, status=45'
  printf, lun, 'clr=2'
  printf, lun, 'corshape, run, rerun, camcol, clr, status=status'
  printf, lun, 'IF status ne 0 THEN exit, status=45'
  printf, lun, 'clr=3'
  printf, lun, 'corshape, run, rerun, camcol, clr, status=status'
  printf, lun, 'IF status ne 0 THEN exit, status=45'
  printf, lun, 'EOF'
  print_status_check, lun
  free_lun, lun
  spawn,'chmod 755 '+corshape_123file
  
  openw, lun, corshape_456file, /get_lun
  printf, lun, head
  printf, lun
  printf, lun, 'idl<<EOF'
  printf, lun, 'setup_mystuff'
  printf, lun, 'run=',rstr
  printf, lun, 'rerun=',rrstr
  printf, lun, 'camcol=4'
  printf, lun, 'clr=1'
  printf, lun, 'corshape, run, rerun, camcol, clr, status=status'
  printf, lun, 'IF status ne 0 THEN exit, status=45'
  printf, lun, 'clr=2'
  printf, lun, 'corshape, run, rerun, camcol, clr, status=status'
  printf, lun, 'IF status ne 0 THEN exit, status=45'
  printf, lun, 'clr=3'
  printf, lun, 'corshape, run, rerun, camcol, clr, status=status'
  printf, lun, 'IF status ne 0 THEN exit, status=45'
  printf, lun, 'camcol=5'
  printf, lun, 'clr=1'
  printf, lun, 'corshape, run, rerun, camcol, clr, status=status'
  printf, lun, 'IF status ne 0 THEN exit, status=45'
  printf, lun, 'clr=2'
  printf, lun, 'corshape, run, rerun, camcol, clr, status=status'
  printf, lun, 'IF status ne 0 THEN exit, status=45'
  printf, lun, 'clr=3'
  printf, lun, 'corshape, run, rerun, camcol, clr, status=status'
  printf, lun, 'IF status ne 0 THEN exit, status=45'
  printf, lun, 'camcol=6'
  printf, lun, 'clr=1'
  printf, lun, 'corshape, run, rerun, camcol, clr, status=status'
  printf, lun, 'IF status ne 0 THEN exit, status=45'
  printf, lun, 'clr=2'
  printf, lun, 'corshape, run, rerun, camcol, clr, status=status'
  printf, lun, 'IF status ne 0 THEN exit, status=45'
  printf, lun, 'clr=3'
  printf, lun, 'corshape, run, rerun, camcol, clr, status=status'
  printf, lun, 'IF status ne 0 THEN exit, status=45'
  printf, lun, 'EOF'
  print_status_check, lun
  free_lun, lun
  spawn,'chmod 755 '+corshape_456file

  print
  print,' CORR_SCRIPTS: Excecution Successful'
  print
  print,' To correct Run ',rstr,' use these commands: '
  print,'   '+call_123script
  print,'   '+call_456script
  print
  print,' When those are done, use these commands to make allcol QA plots: '
  print,'   '+corshape_file
  ;;print,'   '+corshape_123file
  ;;print,'   '+corshape_456file
  print

  return
END 
