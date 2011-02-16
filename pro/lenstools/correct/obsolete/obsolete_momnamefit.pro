PRO psfname, run, camcol, fname

  IF n_params() LT 2 THEN BEGIN
      print,'-Syntax: psfname, run, camcol, fname'
      return
  ENDIF 
  colors = ['u','g','r','i','z']
  rstr = run2string(run)
  cstr = ntostr(camcol)

  fname = 'psffit'+'-'+rstr+'-col'+cstr+'.fit'

return
END 
