FUNCTION postgres_zphot::init, run, rerun=rerun, type=type, princeton=princeton

  nrun = n_elements(run)
  IF nrun EQ 0 THEN BEGIN 
      message,'You must initialize the run/stripe',/inf
      message,"-Syntax: pgz = obj_new('postgres_zphot',run/stripe,rerun=,/princeton)",/inf
      return,0
  ENDIF 
  IF NOT (self->postgres::init()) THEN BEGIN 
      message,'Could not initialize postgres',/inf
      return,0
  ENDIF 
  IF NOT (self->photoz_uchicago::init(type=type)) THEN BEGIN 
      message,'Could not initialize photoz_uchicago',/inf
      return,0
  ENDIF 

  IF keyword_set(princeton) THEN BEGIN 
      self.princeton = 1
      self->reset, stripe=run
  ENDIF ELSE BEGIN 
      self.princeton = 0
      self->reset, run=run, rerun=rerun
  ENDELSE 

  return,1

END 

PRO postgres_zphot::reset, run=run, rerun=rerun, stripe=stripe

  IF n_elements(stripe) NE 0 THEN BEGIN 
      self.stripe = stripe
      self.run = -1
      self.rerun = -1
  ENDIF ELSE IF n_elements(run) NE 0 THEN BEGIN 

      self.stripe = -1
      self.run = run
          
      IF n_elements(rerun) EQ 0 THEN BEGIN 
          sf = obj_new('sdss_files')
          rerun=sf->rerun(run)
      ENDIF 
      
      self.rerun = rerun
  ENDIF 

  self.write_status = -1
  self.stuff_status = -1

END 

FUNCTION postgres_zphot::run
  return,self.run
END 
FUNCTION postgres_zphot::rerun
  return,self.rerun
END 
FUNCTION postgres_zphot::stripe
  return,self.stripe
END 

FUNCTION postgres_zphot::write_status
  return,self.write_status
END 

FUNCTION postgres_zphot::stuff_status
  return,self.stuff_status
END 

PRO postgres_zphot::tabledef, sqlfile

  IF n_elements(sqlfile) EQ 0 THEN BEGIN 
      lun = -1
  ENDIF ELSE BEGIN 
      openw,lun,sqlfile,/get_lun
  ENDELSE 

  struct = self->postgres_struct()
  tags = tag_names(struct)

  ;; inherited function
  coldefs = self->struct2coldefs(struct)

  ;; Add the primary key
  coldefs = [coldefs, 'PRIMARY KEY (photoid)']

  ncoldefs = n_elements(coldefs)
  printf,lun,'CREATE TABLE zphot'
  printf,lun,'('
  FOR i=0L, ncoldefs-2 DO BEGIN 
      printf,lun,coldefs[i]+', '
  ENDFOR 
  printf,lun,coldefs[i]
  printf,lun,');'

  IF n_elements(sqlfile) NE 0 THEN free_lun,lun

END 

FUNCTION postgres_zphot::photoz_file_read, count=count
  IF self.princeton THEN BEGIN 
      zphot = $
        self->photoz_uchicago::princeton_output_read(self.stripe, count=count)
  ENDIF ELSE BEGIN 
      zphot = $
        self->photoz_uchicago::output_read(self.run, $
                                           rerun=self.rerun, count=count)
  ENDELSE 
  return,zphot
END 
FUNCTION postgres_zphot::postgres_file

  IF self.princeton THEN BEGIN 
      input_file = $
        self->photoz_uchicago::princeton_output_file(self.stripe)
      output_file = repstr(input_file, 'etbl.gz', 'pgsql')
  ENDIF ELSE BEGIN 
      ;; This only works for nn, with its xst.gz end
      input_file = $
        self->photoz_uchicago::output_file(self.run, rerun=self.rerun)
      output_file = repstr(input_file, 'xst.gz', 'pgsql')
  ENDELSE 
  return,output_file
END 

FUNCTION postgres_zphot::postgres_struct

  ;; Just the things we need to stuff
  struct = {$
             photoid: 0ULL, $
             photoz_z:-1.0,     $
             photoz_zerr1:-1.0, $
             photoz_zerr2:-1.0, $
             photoz_zerr3:-1.0, $
             photoz_zerr4:-1.0, $
             photoz_zwarning:0b $
           }
  
  return,struct

END 

PRO postgres_zphot::write_postgres_input, status=status

  status = 1

  outfile = self->postgres_file()

  ;; Read in outputs from the photoz code
  zphot = self->photoz_file_read(count=count)

  IF count EQ 0 THEN return

  pgstruct = self->postgres_struct()

  pgstruct = replicate(pgstruct, count)

  pgstruct.photoz_z = zphot.photoz_z
  pgstruct.photoz_zerr4 = zphot.photoz_zerr4
  pgstruct.photoz_zwarning = zphot.photoz_zwarning

  IF NOT self.princeton THEN BEGIN 
      photoid = photoid(zphot)
      pgstruct.photoid = photoid

      pgstruct.photoz_zerr1 = zphot.photoz_zerr1
      pgstruct.photoz_zerr2 = zphot.photoz_zerr2
      pgstruct.photoz_zerr3 = zphot.photoz_zerr3
  ENDIF ELSE BEGIN 
      pgstruct.photoid = zphot.photoid
  ENDELSE 


  print,'Writing to file: ',outfile  
  ;; This is a uproc
  ascii_write, pgstruct, outfile, /bracket_arrays
  ;;write_postgres_input, pgstruct, outfile, error=status

  self.write_status = 0

END 

PRO postgres_zphot::stuff

  pgfile = self->postgres_file()

  IF self.write_status NE 0 THEN BEGIN 
      ;; See if the file exists
      IF NOT fexist(pgfile) THEN BEGIN 
          message,$
            'self.write_status indicates an error ('+$
            ntostr(self.write_status)+') and input file does not exist'+$
            pgfile,/inf
      return
      ENDIF 
  ENDIF 

  IF NOT fexist(pgfile) THEN BEGIN 
      message,'postgres input file does not exist: '+pgfile,/inf
      return
  ENDIF 
  
  self.connect_info = 'user=postgres'

  query = "copy zphot FROM '"+pgfile+"'"
  print
  print,query

  result = self->query(query)

  IF self.query_status NE self->status_val('no_result') THEN BEGIN
      self.stuff_status = 1
  ENDIF ELSE BEGIN 
      self.stuff_status = 0
  ENDELSE 

END 

PRO postgres_zphot::delete_postgres_input
  print
  print,'Deleting file: ',self->postgres_file()
  file_delete, self->postgres_file(), /quiet
END 

PRO postgres_zphot::stuff_scripts, status=status

  status=1

  w=where(!run_status.tsobj_photo_v GT 5.4 AND $
          !run_status.run NE 5042, nw)

  dir = '/net/cheops2/home/esheldon/idlscripts/stuff/zphot/'

  file = dir + 'zphot_stuff.sh'

  openw, lun, file, /get_lun

  printf,lun,'#!/bin/bash'
  printf,lun

  FOR i=0L, nw-1 DO BEGIN 

      run = !run_status[w[i]].run
      rerun = !run_status[w[i]].rerun
      rcstr = 'run '+ntostr(run)

      printf,lun,'idl<<EOF'
      printf,lun,'  run='+ntostr(run)
      printf,lun,'  rerun='+ntostr(rerun)

      printf,lun,'  pgz = obj_new("postgres_zphot", run, rerun=rerun)'
      printf,lun,'  pgz->write_postgres_input'
      printf,lun,'  if pgz->write_status() ne 0 then exit,status=45'
      printf,lun,'  pgz->stuff'
      printf,lun,'  if pgz->stuff_status() ne 0 then exit,status=46'
      printf,lun,'  pgz->delete_postgres_input'
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

  printf,lun,'dt=`date`'
  printf,lun,'message="Finished stuffing database"'
  printf,lun,'echo "$message:  $dt" | mail esheldon@cfcp.uchicago.edu -s "$message"'

  free_lun, lun
  spawn,['chmod','755',file],/noshell

END 


PRO postgres_zphot::princeton_stuff_scripts, status=status

  status=1

  ;; will just skip those that don't exist
  nstripes = 86
  stripes = 1+lindgen(nstripes)

  dir = '/net/cheops2/home/esheldon/idlscripts/stuff/zphot/'

  file = dir + 'princeton_zphot_stuff.sh'

  openw, lun, file, /get_lun

  printf,lun,'#!/bin/bash'
  printf,lun

  FOR i=0L, nstripes-1 DO BEGIN 

      stripe = stripes[i]

      tfile = self->photoz_uchicago::princeton_output_file(stripe)
      IF fexist(tfile) THEN BEGIN 
          sstr = ntostr(stripe)

          rcstr = 'stripe '+sstr
          
          printf,lun,'idl<<EOF'
          printf,lun,'  stripe='+sstr
          printf,lun
          printf,lun,'  pgz = obj_new("postgres_zphot", stripe, /princeton)'
          printf,lun,'  pgz->write_postgres_input'
          printf,lun,'  if pgz->write_status() ne 0 then exit,status=45'
          printf,lun,'  pgz->stuff'
          printf,lun,'  if pgz->stuff_status() ne 0 then exit,status=46'
          printf,lun,'  pgz->delete_postgres_input'
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
      ENDIF 
  ENDFOR 

  printf,lun,'dt=`date`'
  printf,lun,'message="Finished stuffing database"'
  printf,lun,'echo "$message:  $dt" | mail esheldon@cfcp.uchicago.edu -s "$message"'

  free_lun, lun
  file_chmod, file, /u_execute
  status = 0

END 




PRO postgres_zphot__define

  struct = { postgres_zphot, $
             princeton:0, $
             run:0L, $
             rerun:0, $
             stripe: 0, $
             write_status: 0, $
             stuff_status: 0, $
             INHERITS postgres, $
             INHERITS photoz_uchicago $
           }

END 
