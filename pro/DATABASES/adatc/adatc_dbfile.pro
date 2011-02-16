FUNCTION adatc_dbfile, run, camcol, rerun=rerun, status=status, dir=dir

  status = 1
  IF n_params() LT 2 THEN BEGIN 
      print,'-Syntax: database = adatc_dbfile(run, camcol, rerun=, dir=, status=)'
      return,-1
  ENDIF 

  dbfile = 'adatc-'+run2string(run)+'-'+ntostr(camcol)

  IF arg_present(dir) THEN BEGIN 
      dir = sdss_filedir('adatc_db', run, $
                         camcol=camcol, rerun=rerun, /corr, /silent)
  ENDIF 
  
  status=0
  return,dbfile

END 
