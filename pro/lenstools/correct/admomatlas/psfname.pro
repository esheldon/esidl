PRO psfname, run, camcol, fname

  IF n_params() LT 2 THEN BEGIN 
      print,'-Syntax: psfname, run, camcol, fname'
      return
  ENDIF 

  fname = 'psffit-'+run2string(run)+'-'+ntostr(camcol)+'.fit'

  return
END 
