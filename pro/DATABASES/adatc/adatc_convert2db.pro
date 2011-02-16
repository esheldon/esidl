PRO adatc_convert2db, run, camcol, rerun=rerun

    IF n_params() LT 2 THEN BEGIN 
      print,'-Syntax: adatc_convert2db, run, camcol, rerun=rerun'
      return
  ENDIF 

  rstr = run2string(run)
  cstr = ntostr(camcol)

  depth = 10

  adatc_create_dbdfile, run, camcol, rerun=rerun, database=database

  adatc_files = sdss_filelist('adatc', run, camcol=camcol, rerun=rerun, $
                              fields=fields)

  nfiles = n_elements(adatc_files)

  ;; Create the database
  !priv=2
  dbcreate, database, 1, 1

  tm=systime(1)
  FOR i=0L, nfiles-1 DO BEGIN 

      adatc = sdss_read('adatc',run, camcol, field=fields[i], rerun=rerun, $
                        ex_struct={htmindex:0L}, verbose=0, $
                              status=status)
                              
      IF status EQ 0 THEN BEGIN 
          adatc.htmindex = htm_index(adatc.ra, adatc.dec, depth)
          
          print,adatc_files[i]
          sdss_dbbuild_struct, adatc, database, $
            /noindex, /silent, $
            /sort_photoid, /noprompt
      ENDIF 

  ENDFOR 

  ;; Now index the database
  !priv=2
  dbopen, database, 1
  dbindex
  dbclose

  ptime,systime(1)-tm

END 
