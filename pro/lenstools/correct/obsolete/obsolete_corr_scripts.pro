PRO corr_scripts, run, rerun, color_index, email=email, nocheck=nocheck, $
                  psFieldrerun=psFieldrerun, $
                  sdss_data_dir=sdss_data_dir, sdss_shapecorr_dir=sdss_shapecorr_dir

  IF n_params() LT 3 THEN BEGIN
      print,'-Syntax: corr_scripts, run, rerun, color_index, email=email, nocheck=nocheck, psFieldrerun=psFieldrerun, sdss_data_dir=sdss_data_dir, sdss_shapecorr_dir=sdss_shapecorr_dir'
      print,'/nocheck does not check if #atlas files = #tsObj files'
      print,' If you do not give the email address, then email="esheldon@sdss4.physics.lsa.umich.edu"'
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
  camcol = ['1','2','3','4','5','6']
  clr = ntostr(color_index)
  ncol = n_elements(camcol)

  sdssidl_setup, /silent
  IF (n_elements(sdss_data_dir) NE 0) $
    OR ( n_elements(sdss_shapecorr_dir) NE 0) THEN redefine=1 $
  ELSE redefine=0
  IF n_elements(sdss_data_dir) NE 0 THEN !sdss_data_dir = sdss_data_dir
  IF n_elements(sdss_shapecorr_dir) NE 0 THEN !sdss_shapecorr_dir=sdss_shapecorr_dir

  home = !SDSS_SHAPECORR_DIR
  rundir = home+'fscripts/'+rstr
  script_dir = rundir+'/fscripts/'
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
      print,' OK'
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
  print,' OK  '+ntostr(natlas)+'/'+ntostr(ntot)
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
  pnext = 'Directories aleady exist'
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
      IF count2 EQ 0 THEN BEGIN
          print
          pnext=''
          print,' Making directory ',script_dir
          spawn,'mkdir '+script_dir
      ENDIF 
      IF count3 EQ 0 THEN BEGIN
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

  callscript = script_dir+'correct'+rstr+'.sh'

  rotationbase = 'rotate'+rstr
  admomIDLbase = 'admomIDL'+rstr
  admomCbase = 'admomC'+rstr
  convertbase = 'convert'+rstr
  fitmombase = 'fitmom'+rstr
  makecorbase = 'makecor'+rstr
  corshapebase = 'corshape'+rstr

  rotation_file = script_dir + rotationbase+'.sh'

  admomIDL_file = script_dir + admomIDLbase + '.sh'

  admomC_file1 = admomCbase+'_col1.sh'
  admomC_file2 = admomCbase+'_col2.sh'
  admomC_file3 = admomCbase+'_col3.sh'
  admomC_file4 = admomCbase+'_col4.sh'
  admomC_file5 = admomCbase+'_col5.sh'
  admomC_file6 = admomCbase+'_col6.sh'

  admomCfiles = script_dir + [admomC_file1, admomC_file2, admomC_file3, $
                              admomC_file4, admomC_file5, admomC_file6]
  ncolumns = n_elements(admomCfiles)

  convertfile = script_dir + convertbase+'.sh'

  fitmom_file1 = fitmombase+'_col1.sh'
  fitmom_file2 = fitmombase+'_col2.sh'
  fitmom_file3 = fitmombase+'_col3.sh'
  fitmom_file4 = fitmombase+'_col4.sh'
  fitmom_file5 = fitmombase+'_col5.sh'
  fitmom_file6 = fitmombase+'_col6.sh'
  fitmomfiles = script_dir + [fitmom_file1, fitmom_file2, fitmom_file3, $
                              fitmom_file4, fitmom_file5, fitmom_file6]

  makecor_file = script_dir+makecorbase+'.sh'
  corshape_file = script_dir+corshapebase+'.sh'

  print,' Creating call script: '
  print,'    '+callscript

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; The queues for each script
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  rotationq = '30min'
  makeskyq = '4hr'
  admomIDLq = '4hr'               ;admomq may change below
  admomCq = '4hr'
  convertq = '30min'
  fitmomq = '4hr'
  makecorq = '4hr'
  corshapeq = '4hr'

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Open main call script
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  
  openw, call_lun, callscript, /get_lun
  printf, call_lun, head
  printf, call_lun
  printf, call_lun, '. $FBATCH_DIR/bin/fbatch_setpgp.sh'

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; rotation script
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  print,' Creating rotation script: '
  print,'    '+rotation_file
  openw, lun, rotation_file, /get_lun
  printf, lun, head
  printf, lun
  printf, lun, 'idl<<EOF'
  IF redefine THEN BEGIN 
      printf, lun, 'sdssidl_setup, /silent'
      printf, lun, '!SDSS_DATA_DIR = '+sdss_data_dir
      printf, lun, '!SDSS_SHAPECORR_DIR = '+sdss_shapecorr_dir
  ENDIF 
  printf, lun, 'run=',rstr
  printf, lun, 'rerun=',rrstr
  FOR icam=1, 6 DO BEGIN 
      camstr = ntostr(icam)
      printf, lun, 'camcol='+camstr
      printf, lun, 'sdss_survey_rot,run,rerun,camcol'
  ENDFOR 
  printf, lun, 'EOF'
  free_lun, lun
  spawn,'chmod 755 '+rotation_file

  ;; Write to call script

  printf, call_lun, ''
  printf, call_lun, 'echo'
  printf, call_lun, 'echo Sending Script '+rotationbase+'.sh to '+rotationq+' queue'
  printf, call_lun
  printf, call_lun, 'fbatch_sub -q '+rotationq+' -R "fsgi03"' $
               + ' -J '+rotationbase                              $
               + ' -o "' +outdir+rotationbase+'.out"'             $
               + ' -e "' +outdir+rotationbase+'.err"'             $
               + ' -N -u '+email $
               + ' "' +rotation_file+ '"'

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; admom IDL script
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  print,' Creating admom IDL script: '
  print,'    '+admomIDL_file

  openw, lun, admomIDL_file, /get_lun
  printf, lun, head
  printf, lun
  printf, lun, 'idl<<EOF'
  IF redefine THEN BEGIN 
      printf, lun, 'sdssidl_setup, /silent'
      printf, lun, '!SDSS_DATA_DIR = '+sdss_data_dir
      printf, lun, '!SDSS_SHAPECORR_DIR = '+sdss_shapecorr_dir
  ENDIF 
  printf, lun, 'run=',rstr
  printf, lun, 'rerun=',rrstr
  printf, lun, 'psFieldrerun='+psrrstr
  printf, lun, 'start=',ntostr(start)
  printf, lun, 'nf=',ntostr(nframes)
  printf, lun
  printf, lun, 'color_index=',clr
  printf, lun
  printf, lun, 'camcol=1'
  printf, lun, 'admomatlas2, color_index, run, rerun, camcol, start=start,nf=nf,psFieldrerun=psFieldrerun'
  printf, lun, 'camcol=2'
  printf, lun, 'admomatlas2, color_index, run, rerun, camcol, start=start,nf=nf,psFieldrerun=psFieldrerun'
  printf, lun, 'camcol=3'
  printf, lun, 'admomatlas2, color_index, run, rerun, camcol, start=start,nf=nf,psFieldrerun=psFieldrerun'
  printf, lun, 'camcol=4'
  printf, lun, 'admomatlas2, color_index, run, rerun, camcol, start=start,nf=nf,psFieldrerun=psFieldrerun'
  printf, lun, 'camcol=5'
  printf, lun, 'admomatlas2, color_index, run, rerun, camcol, start=start,nf=nf,psFieldrerun=psFieldrerun'
  printf, lun, 'camcol=6'
  printf, lun, 'admomatlas2, color_index, run, rerun, camcol, start=start,nf=nf,psFieldrerun=psFieldrerun'
  printf, lun, 'EOF'

  free_lun, lun
  spawn, 'chmod 755 '+admomIDL_file

  ;; write to call script

  printf, call_lun
  printf, call_lun, 'echo'
  printf, call_lun, 'echo Sending Script '+admomIDLbase+' to '+admomIDLq+' queue'
  printf, call_lun
  printf, call_lun, 'fbatch_sub -q '+admomIDLq+' -R "fsgi03"' $
               + ' -J '+admomIDLbase                              $
               + ' -w "done('+rotationbase+')"'                   $
               + ' -o "' +outdir+admomIDLbase+'.out"'             $
               + ' -e "' +outdir+admomIDLbase+'.err"'             $
               + ' -N -u '+email $
               + ' "' +admomIDL_file+ '"'

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; admomatlas C scripts
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  print,' Creating admom C scripts: '
  colprint,'    '+admomCfiles

  col = [1, 2, 3, 4, 5, 6]
  FOR ifile=0, ncolumns-1 DO BEGIN 
      openw, lun, admomCfiles[ifile], /get_lun
      
      printf, lun, head
      printf, lun
      fetch_dir, run, col[ifile], rerun, tdir, atldir, corrdir=corrdir
      admomatlas_infile, run, rerun, col[ifile], adinfile, adoutfile
      adinfile = corrdir+adinfile & adoutfile = corrdir+adoutfile

      printf, lun, 'admomatlas '+adinfile+' '+adoutfile+' '+atldir

      free_lun, lun
      spawn,'chmod 755 '+admomCfiles[ifile]
  ENDFOR 

  ;; Write to call script

  FOR ifile=0, ncolumns-1 DO BEGIN 

      IF ifile LT 3 THEN BEGIN 
          donestring = ' -w "done('+admomIDLbase+')"'
      ENDIF ELSE BEGIN 
          donestring = ' -w "done('+admomCbase+'_col'+camcol[ifile-3]+')"' 
      ENDELSE 
      printf, call_lun, ''
      printf, call_lun, 'echo'
      printf, call_lun, 'echo Sending Script '+admomCbase+'_col'+camcol[ifile]+'.sh to ' $
                                                                 +admomCq+' queue'
      printf, call_lun
      printf, call_lun, 'fbatch_sub -q '+admomCq+ ' -R "fsgi03"'    $
               + ' -J '+admomCbase+'_col'+camcol[ifile]                    $
               + donestring                 $
               + ' -o "' +outdir+admomCbase+'_col'+camcol[ifile]+'.out"'             $
               + ' -e "' +outdir+admomCbase+'_col'+camcol[ifile]+'.err"'             $
               + ' -N -u '+email $
               + ' "' +admomCfiles[ifile]+ '"'
  ENDFOR 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; convert admomout files from ascii to fits
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  print,' Creating convert script: '
  print,'    '+convertfile

  openw, lun, convertfile, /get_lun
  
  printf, lun, head
  printf, lun
  printf, lun, 'idl<<EOF'
  IF redefine THEN BEGIN 
      printf, lun, 'sdssidl_setup, /silent'
      printf, lun, '!SDSS_DATA_DIR = '+sdss_data_dir
      printf, lun, '!SDSS_SHAPECORR_DIR = '+sdss_shapecorr_dir
  ENDIF 
  FOR ifile=0, ncolumns-1 DO BEGIN
      printf, lun

      fetch_dir, run, col[ifile], rerun, tdir, atldir, corrdir=corrdir
      admomatlas_infile, run, rerun, col[ifile], adinfile, adoutfile
      adoutfile = corrdir+adoutfile

      printf, lun, 'infile = "'+adoutfile+'"'
      printf, lun, 'admomascii2fits, infile'
  ENDFOR 
  printf, lun, 'EOF'
  free_lun, lun
  spawn,'chmod 755 '+convertfile

  ;; write to call script

  printf, call_lun, ''
  printf, call_lun, 'echo'
  printf, call_lun, 'echo Sending Script '+convertbase+' to '+convertq+' queue'
  printf, call_lun
  printf, call_lun, 'fbatch_sub -q '+convertq $
                     + ' -J '+convertbase $
                     + ' -w "done('+admomCbase+'_col1)' $
                     + '  && done('+admomCbase+'_col2)' $
                     + '  && done('+admomCbase+'_col3)' $
                     + '  && done('+admomCbase+'_col4)' $
                     + '  && done('+admomCbase+'_col5)' $
                     + '  && done('+admomCbase+'_col6)"' $
                     + ' -o "' +outdir+convertbase+'.out"' $
                     + ' -e "' +outdir+convertbase+'.err"' $
                     + ' -N -u '+email $
                     + ' "' +convertfile+ '"'

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; fitmom script. New way is slow, break it up to 6 scripts
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  print,' Creating fitmom scripts: '
  colprint,'    '+fitmomfiles

  FOR ifile=0, ncolumns-1 DO BEGIN 
      openw, lun, fitmomfiles[ifile], /get_lun

      printf, lun, head
      printf, lun
      printf, lun, 'idl<<EOF'
      IF redefine THEN BEGIN 
          printf, lun, 'sdssidl_setup, /silent'
          printf, lun, '!SDSS_DATA_DIR = '+sdss_data_dir
          printf, lun, '!SDSS_SHAPECORR_DIR = '+sdss_shapecorr_dir
      ENDIF 
      printf, lun, 'run=',rstr
      printf, lun, 'camcol=',ntostr(col[ifile])
      printf, lun, 'rerun=',rrstr
      printf, lun, 'psFieldrerun='+psrrstr
      printf, lun
      printf, lun, 'fitmom2, run, rerun, camcol,psFieldrerun=psFieldrerun'
      printf, lun, 'EOF'

      free_lun, lun
      spawn,'chmod 755 '+fitmomfiles[ifile]
  ENDFOR 

  ;; write to call script

  FOR ifile=0, ncolumns-1 DO BEGIN 
      
      IF ifile LT 3 THEN BEGIN 
          donestring = ' -w "done('+convertbase+')"'
      ENDIF ELSE BEGIN 
          donestring = ' -w "done('+fitmombase+'_col'+camcol[ifile-3]+')"' 
      ENDELSE 

      printf, call_lun
      printf, call_lun, 'echo'
      printf, call_lun,'echo Sending Script '+fitmombase+'_col'+camcol[ifile]+'.sh to '+$
                                                                 fitmomq+' queue'
      printf, call_lun
      printf, call_lun, 'fbatch_sub -q '+fitmomq+ ' -R "fsgi03"'  $
               + ' -J '+fitmombase+'_col'+camcol[ifile]           $
               + donestring               $
               + ' -o "' +outdir+fitmombase+'_col'+camcol[ifile]+'.out"'  $
               + ' -e "' +outdir+fitmombase+'_col'+camcol[ifile]+'.err"'  $
               + ' -N -u '+email $
               + ' "' +fitmomfiles[ifile]+ '"'
  ENDFOR 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; makecor script
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  print,' Creating makecor script: '
  print,'    '+makecor_file

  openw, lun, makecor_file, /get_lun
  printf, lun, head
  printf, lun
  FOR icol=0, ncol-1 DO BEGIN 
      printf, lun, 'idl<<EOF'
      IF redefine THEN BEGIN 
          printf, lun, 'sdssidl_setup, /silent'
          printf, lun, '!SDSS_DATA_DIR = '+sdss_data_dir
          printf, lun, '!SDSS_SHAPECORR_DIR = '+sdss_shapecorr_dir
      ENDIF 
      printf, lun, 'run=',rstr
      printf, lun, 'rerun=',rrstr
      printf, lun, 'camcol=',camcol[icol]
      printf, lun, 'start=',ntostr(start)
      printf, lun, 'nf=',ntostr(ntot)
      printf, lun
      printf, lun, 'makecor2,run,rerun,camcol,start=start,nframes=nf'
      printf, lun, 'EOF'
  ENDFOR 
  free_lun, lun
  spawn,'chmod 755 '+makecor_file

  ;; write to call script

  donestring = '"done('+fitmombase+'_col1)'                $
         + ' &&  done('+fitmombase+'_col2)'                $
         + ' &&  done('+fitmombase+'_col3)'                $
         + ' &&  done('+fitmombase+'_col4)'                $
         + ' &&  done('+fitmombase+'_col5)'                $
         + ' &&  done('+fitmombase+'_col6)"'

  printf, call_lun
  printf, call_lun, 'echo'
  printf, call_lun, 'echo Sending Script '+makecorbase+'.sh to ' $
                    +makecorq+' queue'
  printf, call_lun
  printf, call_lun, 'fbatch_sub -q '+makecorq+ ' -R "fsgi03"'   $
               + ' -J '+makecorbase                             $
               + ' -w '+donestring                              $  
               + ' -o "' +outdir+makecorbase+'.out"'            $
               + ' -e "' +outdir+makecorbase+'.err"'            $
               + ' -N -u '+email  $
               + ' "' +makecor_file+ '"'

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; corshape script
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  
  print,' Creating corshape script: '
  print,'    '+corshape_file

  openw, lun, corshape_file, /get_lun
  printf, lun, head
  printf, lun
  printf, lun, 'idl<<EOF'
  IF redefine THEN BEGIN 
      printf, lun, 'sdssidl_setup, /silent'
      printf, lun, '!SDSS_DATA_DIR = '+sdss_data_dir
      printf, lun, '!SDSS_SHAPECORR_DIR = '+sdss_shapecorr_dir
  ENDIF 
  printf, lun, 'run=',rstr
  printf, lun, 'rerun=',rrstr
  printf, lun, 'clr=1'
  printf, lun, 'corshape_allcols2, run, rerun, clr'
  printf, lun, 'EOF'
  printf, lun, 'idl<<EOF'
  IF redefine THEN BEGIN 
      printf, lun, 'sdssidl_setup, /silent'
      printf, lun, '!SDSS_DATA_DIR = '+sdss_data_dir
      printf, lun, '!SDSS_SHAPECORR_DIR = '+sdss_shapecorr_dir
  ENDIF 
  printf, lun, 'run=',rstr
  printf, lun, 'rerun=',rrstr
  printf, lun, 'clr=2'
  printf, lun, 'corshape_allcols2, run, rerun, clr'
  printf, lun, 'EOF'
  printf, lun, 'idl<<EOF'
  IF redefine THEN BEGIN 
      printf, lun, 'sdssidl_setup, /silent'
      printf, lun, '!SDSS_DATA_DIR = '+sdss_data_dir
      printf, lun, '!SDSS_SHAPECORR_DIR = '+sdss_shapecorr_dir
  ENDIF 
  printf, lun, 'run=',rstr
  printf, lun, 'rerun=',rrstr
  printf, lun, 'clr=3'
  printf, lun, 'corshape_allcols2, run, rerun, clr'
  printf, lun, 'EOF'
  free_lun, lun
  spawn,'chmod 755 '+corshape_file

  ;; write to call script
  printf, call_lun
  printf, call_lun, 'echo'
  printf, call_lun, 'echo Sending Script '+corshapebase+'.sh to ' $
                    +corshapeq+' queue'
  printf, call_lun
  printf, call_lun, 'fbatch_sub -q '+corshapeq+ ' -R "fsgi03"'  $
               + ' -J '+corshapebase                            $
               + ' -w "done('+makecorbase+')"'                  $
               + ' -o "' +outdir+corshapebase+'.out"'           $
               + ' -e "' +outdir+corshapebase+'.err"'           $
               + ' -N -u esheldon@sdss4.physics.lsa.umich.edu ' $
               + ' "' +corshape_file+ '"'
  free_lun, call_lun
  spawn,'chmod 755 '+callscript
  
  print
  print,' CORR_SCRIPTS: Excecution Successful'
  print
  print,' To correct Run ',rstr,' use these commands: '
  print,' setup fbatch'
  print,' '+callscript
  print

  return
END 
