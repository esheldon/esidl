FUNCTION shapelet_idstring, strOrRun, rerun, camcol, field, id
  
  IF size(strOrRun,/tname) EQ 'STRUCT' THEN BEGIN 

      run    = strOrRun.run
      rerun  = strOrRun.rerun
      camcol = strOrRun.camcol
      field  = strOrRun.field
      id     = strOrRun.id

  ENDIF ELSE IF n_params() EQ 5 THEN BEGIN 

      run = strOrRun

  ENDIF ELSE BEGIN
      print,'-Syntax: idString = shapelet_idstring(struct)'
      print,'                -- OR --'
      print,'         idString = shapelet_idstring(run,rerun,camcol,field,id)'
      return,'ERROR'
  ENDELSE 

  runstr = run2string(run)
  camcolstr = ntostr(long(camcol))
  rerunstr = ntostr(long(rerun))
  fieldstr = field2string(field)
  idstr = id2string(id)

  idstring = runstr+'-'+rerunstr+'-'+camcolstr+'-'+fieldstr+'-'+idstr

  return,idstring

END
