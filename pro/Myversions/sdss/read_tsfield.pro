FUNCTION read_tsfield, run, rerun, camcol, field, _extra=_extra, status=status

  status = 1

  IF n_params() LT 4 THEN BEGIN 
      print,'-Syntax: struct = read_tsField(run, rerun, camcol, field, _extra=_extra, status=)'
      print,'status=0 is success'
      return,-1
  ENDIF 

  rstr = run2string(run)
  rstr2 = ntostr(run)

  rrstr = ntostr(rerun)
  cstr = ntostr(camcol)
  fstr = field2string(field)

  file = $
    'tsField-'+rstr+'-'+cstr+'-'+rrstr+'-'+fstr+'.fit'

  dir = sdssidl_config('data_dir') + $
    rstr2+'/'+rrstr+'/calibChunks/'+cstr+'/'
  
  file = dir + file
  tsFieldStruct = mrdfits(file, 1, silent=silent, _extra=_extra)

  IF size(tsFieldStruct, /tname) NE 'STRUCT' THEN return,-1

  status = 0
  return,tsFieldStruct

END 
