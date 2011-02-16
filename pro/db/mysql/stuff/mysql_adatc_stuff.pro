
PRO mysql_adatc_stuff, run, camcol, rerun=rerun, $
                       start=start, nframes=nframes,$
                       status=status, arrays=arrays


  ;; status is 1 unless we reach end
  status =1

  IF n_params() LT 2 THEN BEGIN ;Help message
      print,'-syntax mysql_adatc_stuff, run, camcol, rerun=rerun, $'
      print,'               start=start, nframes=nframes,$'
      print,'               status=status'
      return 
  ENDIF 

  time=systime(1)

  ;; HTM depth
  htm_depth = 17
  
  sf = obj_new('sdss_files')
  tmpdir = sf->filedir('calibChunks',run, camcol=camcol, /corrected, exists=exists)
  IF NOT exists THEN BEGIN 
      print,'The corrected directory seems to be missing: ',tmpdir
      return
  ENDIF 


  files = sf->filelist('tsObj',run,camcol,fields=fields)
  nfields = n_elements(fields)

  runstr = sf->run2string(run)
  cstr = ntostr(camcol)

  tmpfile = tmpdir + 'adatc-stuff-'+runstr+'-'+cstr
  IF keyword_set(arrays) THEN tmpfile = tmpfile + '.pgsql_arrays'$
  ELSE tmpfile = tmpfile + '.pgsql'

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

              write_pgsql_input, adatc, tmpfile, $
                format=format, append=append, $
                error=error, arrays=arrays
          ENDIF 

      ENDIF 

  ENDFOR 

  ptime,systime(1)-time

  return

  ;; Skip the header. This is necessary because the blank line
  ;; after the header gets read as a record.
  print
  print,'Loading file: ',tmpfile,' into adatc database'
  ignorelines = '30'
  query = $
    "LOAD DATA LOCAL INFILE '"+tmpfile+"' "+$
    "INTO TABLE adatc IGNORE "+ignorelines+" LINES"
  out = mysql_query(query, user='esheldon', pwd='cheopsDB', status=mstatus)
 
  ;; Do not delete file if status problems.  It may just have timed out, and
  ;; we can try to re-stuff the file
  ;; Note: for these queries we expect status=2, NO_RESULT
  IF mstatus EQ 1 THEN BEGIN 
      print,'Error'
      return
  ENDIF 

  file_delete, tmpfile
  ptime,systime(1)-time  


  status=0

  return

END 
