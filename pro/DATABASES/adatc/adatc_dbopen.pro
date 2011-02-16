PRO adatc_dbopen, run, camcol, rerun=rerun, status=status

  IF n_params() LT 2 THEN BEGIN 
      print,'-Syntax: adatc_dbopen, run, camcol, rerun=, status='
      return
  ENDIF 

  database = adatc_dbfile(run, camcol, rerun=rerun, dir=dir)
  database = dir + database
  dbopen, database, unavail=status

END 
