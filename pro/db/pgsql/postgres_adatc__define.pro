FUNCTION postgres_adatc::init, run, camcol, rerun=rerun

  nrun = n_elements(run)
  ncamcol = n_elements(camcol)

  IF nrun EQ 0 OR ncamcol EQ 0 THEN BEGIN 
      message,'You must initialize the run and camcol',/inf
      message,$
        "-Syntax: obj=obj_new('postgres_adatc',run,camcol,rerun=')",/inf
      return,0
  ENDIF 

  tablename = 'adatc'

  test = self->postgres_stuff_camcol::init(run,camcol,tablename, rerun=rerun)
  IF NOT test THEN BEGIN 
      message,'Could not initialize postgres_stuff_camcol',/inf
      return,0
  ENDIF 
  return,1

END 

PRO postgres_adatc::reset, $
                  run=run, camcol=camcol, rerun=rerun

  self->postgres_stuff_camcol::reset, run=run, camcol=camcol, rerun=rerun

END 

PRO postgres_adatc::tabledef, sqlfile

  IF n_elements(sqlfile) EQ 0 THEN BEGIN 
      lun = -1
  ENDIF ELSE BEGIN 
      openw,lun,sqlfile,/get_lun
  ENDELSE 

  struct = make_corrected_struct()

  coldefs = self->struct2coldefs(struct)

  ;; Add the primary key
  coldefs = [coldefs, 'PRIMARY KEY (photoid)']

  ncoldefs = n_elements(coldefs)
  printf,lun,'CREATE TABLE adatc'
  printf,lun,'('
  FOR i=0L, ncoldefs-2 DO BEGIN 
      printf,lun,coldefs[i]+', '
  ENDFOR 
  printf,lun,coldefs[i]
  printf,lun,');'

  IF n_elements(sqlfile) NE 0 THEN free_lun,lun
  
END 

PRO postgres_adatc::write_input

  ;; status is 1 unless we reach end
  self.write_status = 1

  time=systime(1)

  ;; HTM depth
  htm_depth = 17
  
  sf = obj_new('sdss_files')

  run = self.run
  rerun = self.rerun
  camcol = self.camcol

  runstr = sf->run2string(run)
  cstr = ntostr(camcol)

  files = sf->filelist('tsObj',run,camcol,rerun=rerun,fields=fields)
  nfields = n_elements(fields)

  tmpfile = self->input_file()

  file_delete, tmpfile, /quiet

  print
  print,'Writing to temporary file: ',tmpfile


  FOR fi = 0L, nfields-1 DO BEGIN 

      field = fields[fi]
      fstr=ntostr(field)

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; Read tsObj file
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      tsobj_struct = sdss_read('tsobj',run, camcol, rerun=rerun,$
                                     field=field, $
                                     taglist=taglist, $
                                     ex_struct=newstruct, verbose=0, $
                                     status=ts_status)


      IF ts_status EQ 0 THEN BEGIN 

          nps = n_elements(tsobj_struct)
          print
          print,'Run: ',runstr,' Camcol: ',cstr,' Field: ',fstr,' Nobj: ',$
            ntostr(nps)

          adatc = make_corrected(tsobj_struct, used, $
                                 rotstruct=rotstruct, status=cstatus)

          IF cstatus EQ 0 THEN BEGIN 

              nused = n_elements(used)
              print,'Writing '+ntostr(nused)

              adatc.htm_index = htm_index(tsobj_struct[used].ra, $
                                          tsobj_struct[used].dec, $
                                          htm_depth)
              
              ;; write file
              IF fexist(tmpfile) THEN append=1 ELSE append=0

              ascii_write, adatc,  tmpfile, /bracket_arrays, append=append,$
                status=status

              adatc = 0

              IF status NE 0 THEN BEGIN 
                  message,'Error writing input file',/inf
                  return
              ENDIF 

;              postgres_write_input, adatc, tmpfile, $
;                format=format, append=append, $
;                error=error
          ENDIF 

      ENDIF 

  ENDFOR 

  ptime,systime(1)-time
  self.write_status = 0
  return


END 

PRO postgres_adatc__define

  struct = {$
             postgres_adatc, $
             INHERITS postgres_stuff_camcol $
           }

END 
