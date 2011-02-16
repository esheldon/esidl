PRO skyname, run, camcol, skyfile

  IF n_params() LT 2 THEN BEGIN 
      print,'-Syntax: skyname, run, camcol, skyfile'
      return
  ENDIF 

  rstr = run2string(run)
  cstr = ntostr(camcol)
  skyfile = 'Sky-'+rstr+'-'+cstr+'.fit'

  return
END 
