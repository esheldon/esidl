PRO mysql_stuff_scripts, scriptnum, adatc=adatc

  w=where(!run_status.tsobj_photo_v GT 5.4, nw)
  eachdef = nw/4
  remainder = nw - (eachdef*4)
  CASE scriptnum OF
      1: w = w[0:eachdef-1]
      2: w = w[eachdef:2*eachdef-1]
      3: w = w[2*eachdef:3*eachdef-1]
      4: w = w[3*eachdef:nw-1]
      ELSE: message,'unknown script number: '+ntostr(scriptnum)
  ENDCASE 

  nw = n_elements(w)
  dir = '/net/cheops2/home/esheldon/idlscripts/stuff/'

  IF keyword_set(adatc) THEN BEGIN 
      dir = dir + 'adatc/'
      type = 'adatc'
  ENDIF ELSE BEGIN 
      dir = dir + 'tsObj/'
      type = 'tsobj'
  ENDELSE 
  file = dir + type+'_stuff'+ntostr(scriptnum)+'.sh'

  openw, lun, file, /get_lun

  printf,lun,'#!/bin/bash'
  printf,lun

  FOR i=0L, nw-1 DO BEGIN 

      run = !run_status[w[i]].run
      rerun = !run_status[w[i]].rerun

      FOR camcol=1,6 DO BEGIN 

          rcstr = 'run '+ntostr(run)+' camcol '+ntostr(camcol)

          printf,lun,'idl<<EOF'
          printf,lun,'  run='+ntostr(run)
          printf,lun,'  rerun='+ntostr(rerun)
          printf,lun,'  camcol='+ntostr(camcol)
          printf,lun
          printf,lun,'  mysql_'+type+'_stuff, run, camcol, rerun=rerun, status=status'
          printf,lun,'  if status ne 0 then exit,status=45'
          printf,lun,'EOF'
          printf,lun,'status=$?'
          printf,lun,'if [ $status -ne 0 ]'
          printf,lun,'then'
          printf,lun,'    echo Error in `basename $0` for '+rcstr
          printf,lun,'    err="Stuff error for '+rcstr+'"'
          printf,lun,'    echo "$err" | mail esheldon@cfcp.uchicago.edu -s "$err"'
          printf,lun,'fi'
          printf,lun
          printf,lun

      ENDFOR 

  ENDFOR 

  printf,lun,'dt=`date`'
  printf,lun,'message="Finished stuffing database: '+ntostr(scriptnum)+'"'
  printf,lun,'echo "$message:  $dt" | mail esheldon@cfcp.uchicago.edu -s "$message"'


  free_lun, lun

  spawn,['chmod','755',file],/noshell

END 
