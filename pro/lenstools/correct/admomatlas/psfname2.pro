PRO psfname2, run, rerun, camcol, field, fname

  IF n_params() LT 2 THEN BEGIN 
      print,'-Syntax: psfname2, run, rerun, camcol, field, fname'
      return
  ENDIF 

  runstr = run2string(run)
  camcolstr = ntostr(long(camcol))
  rerunstr = ntostr(long(rerun))
  fieldstr = field2string(field)
  fname = 'psffit-'+runstr+'-'+camcolstr+'-'+rerunstr+'-'+fieldstr+'.fit'

  return
END 
