FUNCTION postgres_stuff_camcol::init, run, camcol, tablename, rerun=rerun

  nrun = n_elements(run)
  ncamcol = n_elements(camcol)
  ntable = n_elements(tablename)

  IF nrun EQ 0 OR ncamcol EQ 0 OR ntable EQ 0 THEN BEGIN 
      message,'You must initialize the run and camcol and tablename',/inf
      message,"-Syntax: pga = obj_new('--', run, camcol, tablename, rerun=)",/inf
      return,0
  ENDIF 
  IF NOT (self->postgres::init()) THEN BEGIN 
      message,'Could not initialize postgres',/inf
      return,0
  ENDIF 

  ;; To avoid conflict with reset in inherited classes
  self->postgres_stuff_camcol::reset, $
    run=run, camcol=camcol, tablename=tablename, rerun=rerun

  return,1
END 

PRO postgres_stuff_camcol::reset, $
                         run=run, camcol=camcol, tablename=tablename, $
                         rerun=rerun

  nrun = n_elements(run)
  ncamcol = n_elements(camcol)
  ntable = n_elements(tablename)

  IF nrun NE 0 THEN BEGIN 
      self.run = run
      
      IF n_elements(rerun) EQ 0 THEN BEGIN 
          sf = obj_new('sdss_files')
          rerun=sf->rerun(run)
      ENDIF 

      self.rerun = rerun
  ENDIF 

  IF ncamcol NE 0 THEN BEGIN 
      self.camcol = camcol
  ENDIF 

  IF n_elements(tablename) NE 0 THEN BEGIN 
      self.tablename = tablename
  ENDIF 
  
  self.stuff_status = -1

END 

FUNCTION postgres_stuff_camcol::run
  return,self.run
END 
FUNCTION postgres_stuff_camcol::rerun
  return,self.rerun
END 
FUNCTION postgres_stuff_camcol::camcol
  return,self.camcol
END 
FUNCTION postgres_stuff_camcol::tablename
  return,self.tablename
END 

FUNCTION postgres_stuff_camcol::write_status
  return,self.write_status
END 

FUNCTION postgres_stuff_camcol::stuff_status
  return,self.stuff_status
END 


FUNCTION postgres_stuff_camcol::input_file

  host = getenv('HOST')

;  dir = '/net/cheops1/data6/db/tmp/'
  dir = '/net/'+host+'/data0/esheldon/tmp/'
  sf = obj_new('sdss_files')
  runstr = sf->run2string(self.run)
  cstr = ntostr(self.camcol)
  file = dir + $
    self.tablename+'-stuff-'+runstr+'-'+cstr+'.pgsql'

  return,file

END 

PRO postgres_stuff_camcol::stuff

  ;; note, -1 is allowed since we may have created the file before
  IF self.write_status EQ 1 THEN BEGIN 
      message,$
        'self.write_status indicates an error: '+ntostr(self.write_status),/inf
      return
  ENDIF 

  tm = systime(1)

  ;; Need to be superuser to stuff into database. This assumes that
  ;; your ~/.pgpass file has the password
  self.connect_info = 'user=postgres'

  file = self->input_file()

  IF NOT fexist(file) THEN BEGIN 
      message,'File not found: '+file,/inf
      return
  ENDIF 

  print
  query = "COPY "+self.tablename+" FROM '"+file+"'"
  print,query

  self->query,query

  IF self.query_status NE self->status_val('no_result') THEN BEGIN
      self.stuff_status = 1
  ENDIF ELSE BEGIN 
      self.stuff_status = 0
  ENDELSE 

  ptime,systime(1)-tm

END 

PRO postgres_stuff_camcol::delete_input
  print
  print,'Deleting file: ',self->input_file()
  file_delete, self->input_file(), /quiet
END 


PRO postgres_stuff_camcol__define

  struct = {$
             postgres_stuff_camcol, $
             run:0L, $
             camcol: 0, $
             rerun:0, $
             tablename:'', $
             write_status: 0, $
             stuff_status: 0, $
             INHERITS postgres $
           }

END 
