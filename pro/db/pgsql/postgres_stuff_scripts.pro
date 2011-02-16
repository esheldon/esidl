PRO postgres_stuff_scripts_getruns, scriptnum, runs, reruns

  ;; 5042 just has empty files, how stupid can you get
  w=where(!run_status.tsobj_photo_v GT 5.4 AND $
          !run_status.run NE 5042, nw)
  eachdef = nw/4
  remainder = nw - (eachdef*4)
  CASE scriptnum OF
      1: w = w[0:eachdef-1]
      2: w = w[eachdef:2*eachdef-1]
      3: w = w[2*eachdef:3*eachdef-1]
      4: w = w[3*eachdef:nw-1]
      ELSE: message,'unknown script number: '+ntostr(scriptnum)
  ENDCASE 

  runs = !run_status[w].run
  reruns = !run_status[w].rerun

END 



PRO postgres_stuff_scripts, tablename, scriptnum

  IF n_params() LT 2 THEN BEGIN 
      print,'-Syntax: postgres_stuff_scripts, tablename, scriptnum'
      return
  ENDIF 

  postgres_getruns, scriptnum, runs, reruns, /new

  nrun = n_elements(runs)

  dir = '/net/cheops2/home/esheldon/idlscripts/stuff/'+tablename+'/'

  IF NOT fexist(dir) THEN BEGIN 
      print,'Making script directory'
      spawn,'mkdir '+dir
  ENDIF 

  file = dir + tablename + '_stuff'+ntostr(scriptnum)+'.sh'

  objtype = 'postgres_'+tablename

  openw, lun, file, /get_lun

  printf,lun,'#!/bin/bash'
  printf,lun

  FOR i=0L, nrun-1 DO BEGIN 

      run = runs[i]
      rerun = reruns[i]

      FOR camcol=1,6 DO BEGIN 

          rcstr = tablename + ' run '+ntostr(run)+' camcol '+ntostr(camcol)

          printf,lun,'idl<<EOF'
          printf,lun,'  run='+ntostr(run)
          printf,lun,'  rerun='+ntostr(rerun)
          printf,lun,'  camcol='+ntostr(camcol)
          printf,lun
          printf,lun,'  pga=obj_new("'+objtype+'", run, camcol, rerun=rerun)'
          printf,lun,'  pga->write_input'
          printf,lun,'  if pga->write_status() ne 0 then exit,status=45'
          printf,lun,'  pga->stuff'
          printf,lun,'  if pga->stuff_status() ne 0 then exit,status=46'
          printf,lun,'  pga->delete_input'
          printf,lun,'EOF'
          printf,lun,'status=$?'
          printf,lun,'if [ $status -eq 45 ]'
          printf,lun,'then'
          printf,lun,'    echo Error in `basename $0` for '+rcstr
          printf,lun,'    err="Write error for '+rcstr+'"'
          printf,lun,'    echo "$err" | mail esheldon@cfcp.uchicago.edu -s "$err"'
          printf,lun,'fi'
          printf,lun,'if [ $status -eq 46 ]'
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
  printf,lun,'message="Finished stuffing '+tablename+' database: '+$
    ntostr(scriptnum)+'"'
  printf,lun,'echo "$message:  $dt" | mail esheldon@cfcp.uchicago.edu -s "$message"'


  free_lun, lun

  spawn,['chmod','755',file],/noshell

END 
